# De novo sequencing for ion-trap HCD spectra (Casanovo fine-tuning)

Stage 3 (`peaksparity/denovo.py`). PEAKS's de novo result feeds its `denovoOnly` export and the
Tag Length / CAA (%) columns of the DB tables (docs/PEAKS_ALGORITHMS.md section 1). Our open
substitute is Casanovo 5.2.1. Its released model (`casanovo_orbitrap_v5-2-0.ckpt`) was trained on
Orbitrap MS2 and does poorly on the low-resolution ion-trap (ITMS) HCD spectra of the Jurkat
reference run. This file records how we fine-tuned it on ion-trap HCD PSMs from our own pipeline
and what it reaches. The PEAKS binary was never used; PEAKS's de novo calls were only read
from its exported CSVs.

Scripts: `scripts/denovo_make_mgf.py` (training/test MGFs), `scripts/denovo_finetune.py`
(training under a wall-clock budget), `scripts/denovo_eval.py` (evaluation). Working data:
`/drive_a/drew/peaks-parity-work/denovo_finetune/`.

## API (`peaksparity/denovo.py`)

| function | purpose |
|---|---|
| `to_casanovo(peptide)` / `from_casanovo(seq)` | PEAKS `C(+57.02)`, `M(+15.99)`, `N(+0.98)`, `Q(+0.98)`, `(+42.01)` <-> Casanovo tokens `C[Carbamidomethyl]`, `M[Oxidation]`, `N[Deamidated]`, `Q[Deamidated]`, `[Acetyl]-`; other mods -> None |
| `training_psms(psms_csv, min_score=30)` | labels from a PEAKS-format `db.psms.csv`: -10LgP >= 30, in-vocabulary mods, I/L variant rows merged, best PSM per (peptide, z) |
| `write_annotated_mgf(entries, mzml, out)` | MGF with `SEQ=` lines, PEPMASS/CHARGE from the PSM (feature-corrected m/z), SCANS = raw scan, peaks from the mzML |
| `finetune_checkpoint(src, dst)` | makes the released checkpoint trainable (see "Fine-tuning" below) |
| `sequence_mzml(mzml, model_ckpt, out_mztab, scans=None, precursors=None, n_beams=1, top_match=1)` | wraps `casanovo sequence`; returns DataFrame(scan, sequence (ProForma), score, per_residue_scores, peptide (PEAKS notation), charge, exp_mz, calc_mz, precursor_ok, ppm). `scans`/`precursors` write a temporary MGF (subset of scans, or corrected precursor m/z and z per row; the precursor table's other columns, e.g. feature_id, rt_min, area, are carried to the output rows). With `top_match > 1` the best-scoring candidate within 15 ppm (isotope 0/1) is kept |
| `read_mztab(path)` | Casanovo mzTab -> same DataFrame (raw scan from the mzML native id or the MGF `SCANS=`) |
| `evaluate(pred, truth, key)` | exact modified sequence (I != L), exact with I = L, Casanovo peptide/AA precision and AA recall, accuracy-vs-score curve |
| `denovo_only_table(denovo_df, accepted_psm_scans, source_file, precursors=None, min_alc=50, precursor_tol_ppm=15, fasta=None)` | rows in the PEAKS `<Sample>.db.denovoOnly.csv` layout, ready for `peaks_export.write_table(..., table="denovo_only")` |

Score conventions. Casanovo's peptide score is the product of per-residue softmax
probabilities (in 5.2.1 there is no precursor-mass penalty in de novo mode; the
`precursor_mass_tol` config key only affects DB-search mode). `denovo_only_table` uses
ALC (%) = mean(per-residue score) x 100 and local confidence (%) = per-residue score x 100
rounded to integers, as requested; this is a stand-in for PEAKS's candidate-agreement LC
(PEAKS_ALGORITHMS 1.2), not a reimplementation of it.

`denovo_only_table` details, from the reference export (Jurkat, 101,551 rows):
- a spectrum is excluded if its scan (or (scan, Feature Id) when feature ids are given) has an
  accepted DB PSM;
- rows with ALC < 50 are dropped; rows are sorted by ALC descending (true in the reference);
- every PEAKS row has |ppm| <= 15.1 (the precursor tolerance), and ppm is uniformly spread over
  [-15, 15], i.e. PEAKS only reports mass-consistent sequences. Casanovo does not enforce this,
  so rows outside `precursor_tol_ppm` are dropped; ppm = (observed - theoretical mass) /
  theoretical x 1e6, which reproduces the reference values;
- `Tag length` = `Length`, `tag(>=0%)` = `Peptide` (threshold 0 %), `mode` = `HCD`;
- `PTM`: one name per modified site, sorted alphabetically and joined with `; `
  (`Carbamidomethylation; Oxidation (M)` even when M precedes C);
- `Accession`: filled only when a FASTA is given (exact I=L substring match). PEAKS fills it for
  2.5 % of rows; its rule is not established;
- `Area`, `RT`, `Feature Id` come from `precursors` (feature table) when given, otherwise
  empty / 0 (PEAKS: 67 % of reference rows have Feature Id 0 and no Area).
- PEAKS reports lengths 4-57; Casanovo's default `min_peptide_len` is 6 (1.1 % of PEAKS rows
  are shorter).

## Training data

Pipeline 1 %-FDR PSMs (Comet + rescoring, semi-tryptic, CAM fixed) of the PXD004940 HeLa
ion-trap HCD runs (`151118_Q_Hela_*_HCD_IT_*`, Orbitrap Fusion Lumos per the mzML, ITMS rapid scan, HCD 28),
filtered by `training_psms` (-10LgP >= 30, unique (peptide, z) per run):

| run | role | spectra |
|---|---|---|
| 03_01 (500 ng) | train | 37,241 |
| 03_02 (500 ng) | train | 36,385 |
| 05_01 (1000 ng) | train | 41,777 |
| 03_03 (500 ng) | validation (whole run held out) | 37,161 |
| 05_02, 05_03 (1000 ng) | train, second model only | 40,434 + 41,258 |

When the MGFs were built, four runs were finished; 05_02 and 05_03 finished while the first
model was training and were used only for the second model (`ft5runs`).

Test set (never trained on, never used for model selection): the Jurkat run
DrewR_Adam_91526_DDA_OT_4 (also Orbitrap Fusion Lumos; ITMS turbo scan, HCD 30), PEAKS DB PSMs from
`Analysis 1 - Sample 4/db.psms.csv`: 7,378 spectra, one per (scan, peptide, z) with I/L
variants merged, using PEAKS's corrected precursor m/z. No score filter.

## Fine-tuning

`casanovo train --model casanovo_orbitrap_v5-2-0.ckpt` fails in 5.2.1 with
`TypeError: Adam.__init__() got an unexpected keyword argument 'precursor_mass_tol'`: the
released checkpoint stores two DB-search hyperparameters that `Spec2Pep` forwards to the
optimizer. `finetune_checkpoint` writes a copy without them (weights untouched) and training
then runs normally from the released weights. The vocabulary of the checkpoint and of the
default config are identical (28 tokens including plain `C`), so no token remapping occurs.
Casanovo replaces I by L in training labels and never predicts I.

Config (Casanovo defaults except): learning_rate 1e-4 (the released checkpoint records 2.5e-4;
Casanovo's from-scratch default is 5e-4), warmup_iters 1000, cosine_schedule_period_iters = planned steps, train_batch_size 128,
val_check_interval 3000 (2500), calculate_precision true, precision 32-true, devices 1.
Preprocessing: min_mz 50, max_peaks 150, max_charge 4, remove_precursor_tol 2 Da (all Casanovo
defaults, which already suit low-resolution data). `dropout` and other architecture keys cannot
be changed when fine-tuning (Casanovo keeps the checkpoint's). One RTX A6000, ~2.4 steps/s
(~300 spectra/s) at batch 128, 8 GB GPU memory; validation with decoding of the 37k held-out
spectra takes ~3 min.

Runs (all from the released weights, same data pipeline; wall clock 03:22-07:15, 3.9 h;
7.2 GPU-hours on two A6000s run in parallel):

| run | train runs | lr | steps | wall time | best valid CE loss (step) | valid peptide precision (Casanovo) |
|---|---|---|---|---|---|---|
| ft_lr1e-4 | 3 | 1e-4 | ~15,200 | 1.8 h | 0.472 (6,000) | 37.6 % |
| ft_lr3e-4 | 3 | 3e-4 | ~7,300 (stopped early) | 0.9 h | 0.469 (6,000) | 37.1 % |
| ft_lr3e-5 | 3 | 3e-5 | ~17,300 | 2.4 h | 0.503 (15,000) | 36.4 % |
| **ft5runs** | 5 | 1e-4 | ~15,500 | 2.05 h | **0.424 (15,000)** | **41.4 %** |

With 3 training runs the validation loss is lowest after 6 epochs (6,000 steps) and then rises
(0.529 at 12k, 0.568 at 15k) while training loss keeps falling: the model memorises the
~115k training spectra. Two more runs (197k spectra) delay this. The model selected by
validation loss (never by the test set) is ft5runs, step 15,000:

`/drive_a/drew/peaks-parity-work/denovo_finetune/models/casanovo_iontrap_hcd_ft5runs_step15000.ckpt`

(also kept: `casanovo_iontrap_hcd_ft3runs_step6000.ckpt`, the best 3-run model).
The start checkpoint is `.../denovo_finetune/casanovo_orbitrap_v5-2-0.finetune_start.ckpt`.

## Results

Metrics. `exact`: predicted modified sequence identical to the reference (I != L; Casanovo
never predicts I, so every I-containing reference counts as wrong). `exact I=L`: the same with
I and L identified; this is the number to use. `Casanovo pep`, `AA prec`, `AA rec`: Casanovo's
own mass-based definitions (`casanovo.denovo.evaluate`; residues match within 0.1 Da, so e.g.
N(+0.98) = D). Spectra without a prediction count as wrong. AUC = area under the
accuracy (I=L)-vs-coverage curve when spectra are ranked by Casanovo score.
"beam 5 + precursor" = `n_beams 5, top_match 5` and keep the best-scoring candidate within
15 ppm (isotope 0/1) of the precursor (`pick_precursor_consistent`); Casanovo 5.2.1 itself does
not use the precursor mass to choose among beams.

### Jurkat test run, 7,378 PEAKS PSM spectra (PEAKS precursor m/z and z)

| model | decoding | exact | exact I=L | Casanovo pep | AA prec | AA rec | AUC |
|---|---|---|---|---|---|---|---|
| released | beam 1 | 7.0 % | 12.5 % | 13.2 % | 40.6 % | 39.6 % | 0.35 |
| released | beam 5 + precursor | 10.8 % | 18.4 % | 18.9 % | 46.7 % | 45.6 % | 0.43 |
| ft5runs (selected) | beam 1 | 7.5 % | 14.2 % | 15.6 % | 42.0 % | 40.6 % | 0.38 |
| ft5runs (selected) | beam 5 + precursor | 11.6 % | 20.4 % | 21.1 % | 46.4 % | 45.0 % | 0.45 |
| ft3runs step 6000 | beam 1 | 8.1 % | 15.4 % | 16.7 % | 44.2 % | 42.3 % | |
| ft3runs step 6000 | beam 5 + precursor | 12.9 % | 22.5 % | 23.1 % | 48.8 % | 47.0 % | 0.48 |

Accuracy (I=L) at a score threshold keeping the top x % of spectra:

| model, decoding | 10 % | 25 % | 50 % | 75 % | 100 % |
|---|---|---|---|---|---|
| released, beam 1 | 74.8 % | 46.2 % | 24.9 % | 16.6 % | 12.5 % |
| released, beam 5 + precursor | 80.4 % | 56.2 % | 35.8 % | 24.6 % | 18.4 % |
| ft5runs, beam 1 | 79.1 % | 50.9 % | 28.1 % | 18.9 % | 14.2 % |
| ft5runs, beam 5 + precursor | 85.8 % | 57.9 % | 37.9 % | 27.1 % | 20.4 % |
| ft3runs step 6000, beam 5 + precursor | 88.3 % | 62.1 % | 41.1 % | 30.0 % | 22.5 % |

Plot: `/drive_a/drew/peaks-parity-work/denovo_finetune/eval/accuracy_vs_score_jurkat.png`;
curves: `eval/<model>.peaks[.b5t5].curve.tsv`.

Sequencing the whole mzML with its own precursor records (189k MS2 spectra; evaluated on the
7,322 scans with a PEAKS PSM, best PSM per scan), beam 1: released 12.2 % exact I=L (AA prec
40.1 %), ft5runs 13.2 % (AA prec 41.2 %).

Every checkpoint saved during training (every 30 min, best-so-far by validation loss) was
evaluated on Jurkat (beam 1, exact I=L): 14.2-16.3 % for all fine-tuned models, against
12.5 % for the released model, with no trend towards higher values as validation loss falls;
the earliest checkpoints (step 2,500-3,000) are among the best (16.2 %).

### Held-out HeLa run 03_03 (37,161 spectra), beam 1

| model | exact | exact I=L | Casanovo pep | AA prec | AA rec |
|---|---|---|---|---|---|
| released | 9.2 % | 16.8 % | 17.6 % | 45.0 % | 45.3 % |
| ft3runs step 6000 | 20.5 % | 37.6 % | 37.8 % | 61.2 % | 61.3 % |
| ft5runs | 22.9 % | 42.1 % | 42.3 % | 64.0 % | 64.2 % |

### Why the gain does not transfer to Jurkat

Split by whether the reference peptide occurs in the training labels (exact I=L, beam 1; the
first percentage / released value refers to the 3-run labels, the second to the 5-run labels):

| test set | peptides | released | ft3runs step 6000 | ft5runs |
|---|---|---|---|---|
| HeLa 03_03 | seen in training (86 % / 91 %) | 18.7 % / 18.1 % | 41.7 % | 45.3 % |
| HeLa 03_03 | novel (14 % / 9 %) | 5.4 % / 3.7 % | 13.1 % | 10.5 % |
| Jurkat | seen in training (63 % / 66 %) | 13.7 % | 18.7 % | 17.4 % |
| Jurkat | novel | 10.4 % | 9.9 % | 7.8 % |

Within the HeLa series fine-tuning also helps on unseen peptides (5.4 -> 13.1 %), so the model
learns the acquisition, not only the peptides. On Jurkat it helps only for peptides it has seen,
and slightly hurts unseen ones. The two acquisitions differ: HeLa MS2 is ITMS *rapid* scan,
HCD 28, first mass 140, fragment m/z offset -0.065 Da, median |error| 0.035 Da; Jurkat is ITMS
*turbo* scan, HCD 30, first mass 113 (7.7 % of intensity below m/z 140), offset +0.005 Da,
median |error| 0.059 Da, ~880 vs ~520 centroids per spectrum. Two checks did not explain the gap:
shifting Jurkat fragments by -0.07 Da (15.6 % vs 15.4 %) and cutting them at m/z 140
(15.2 % vs 15.4 %). The composition also differs (Jurkat: 18 % Met-oxidised PSMs, almost no
Cys peptides; HeLa: ~2 % oxidised, carbamidomethyl Cys common): on Jurkat, oxidised peptides are
predicted correctly half as often (8.8 % vs 17.9 %, ft3runs step 3000). Charge 3+ is nearly
never right (1.9 %), charge 2+ 18 %.

### Comparison with PEAKS's own de novo (denovoOnly spectra)

PEAKS exports de novo sequences only for spectra without a DB hit (`denovoOnly`, 101,551 rows,
ALC >= 50); for DB-identified scans only Tag Length / CAA survive, so a direct comparison on
those is not possible. Agreement with PEAKS's de novo peptide on the denovoOnly spectra
(PEAKS precursor m/z and z per row; I=L; PEAKS de novo never reports I):

| model, decoding | all rows | PEAKS ALC 80-90 (5,577) | ALC >= 90 (499) | AA prec vs PEAKS |
|---|---|---|---|---|
| released, beam 1 | 0.19 % | 1.5 % | 11.8 % | 16.6 % |
| ft5runs, beam 1 | 0.25 % | 1.6 % | 11.6 % | 17.2 % |
| released, beam 5 + precursor | 0.37 % | 3.0 % | 16.6 % | 18.0 % |
| ft5runs, beam 5 + precursor | 0.45 % | 3.2 % | 16.2 % | 17.1 % |

Agreement is near zero because neither side is usually right on these spectra. For the 5,135
denovoOnly rows that our own pipeline identifies by database search (same scan, z and m/z;
1 % FDR; `eval/denovo_only_vs_pipeline_db.tsv`), the fraction whose de novo sequence equals the
pipeline's DB peptide (I=L) is: PEAKS de novo 4.2 %, released Casanovo 4.9 %, ft5runs 6.6 %
(beam 1); in PEAKS ALC >= 90 (106 rows): PEAKS 38.7 %, released 35.8 %, ft5runs 49.1 %.
With beam 5 + precursor selection: released 8.3 %, ft5runs 10.8 % overall; PEAKS ALC >= 80
(882 rows): PEAKS 16.0 %, released 24.5 %, ft5runs 30.2 %; ALC >= 90: 38.7 / 48.1 / 57.5 %.
By this measure (pipeline IDs as reference, which carry their own ~1 % error) the fine-tuned
model with precursor-aware beam search is right about 2.5 times as often as PEAKS de novo on
these spectra, and PEAKS's ALC is a useful confidence for both.

## denovoOnly export from our model

```python
from peaksparity import denovo
from peaksparity.peaks_export import write_table
# prec: one row per precursor (chimeric scans may have several): scan, charge, mz, rt_min, area, feature_id
pred = denovo.sequence_mzml(mzml, ckpt, "out/run.denovo.mztab", precursors=prec, n_beams=5, top_match=5)
# accepted: (scan, feature_id) pairs (or plain scans) with an accepted DB PSM
rows = denovo.denovo_only_table(pred, accepted, "Run.raw")
write_table(rows, "Sample 1.db.denovoOnly.csv", table="denovo_only")
```

Rows produced for the 101,551 PEAKS denovoOnly spectra of Jurkat (PEAKS precursors), versus
101,551 in PEAKS's export by construction:

| model, decoding | ALC >= 50 | ALC >= 50 and \|ppm\| <= 15 |
|---|---|---|
| released, beam 1 | 52,144 | 2,717 |
| ft5runs, beam 1 | 85,697 | 3,511 |
| released, beam 5 + precursor | 70,795 | 8,052 |
| ft5runs, beam 5 + precursor | 99,056 | 9,417 |

The mass filter is what separates us from PEAKS: every PEAKS de novo row fits the precursor
within 15 ppm, but even with 5 beams Casanovo finds a mass-consistent sequence for only ~9 %
of these spectra. Reproducing PEAKS's row count would need mass-constrained decoding (e.g.
Casanovo's DB-search-style scoring of mass-matched candidates, or a larger beam with a faster
implementation) and an LC calibrated to PEAKS's ALC distribution (PEAKS median 63 %; ours is
higher: 98 % of ft5runs beam-5 predictions have ALC >= 50).

## Open issues

- Transfer to other ion-trap acquisitions is weak (see above). Training data with the target's
  acquisition (ITMS turbo, HCD 30, same instrument) is what is missing; PXD004940 is rapid scan.
  Mixing Orbitrap MassIVE-KB spectra into fine-tuning (to limit forgetting) was not tried.
- Model selection used a held-out run of the same HeLa series; its peptides overlap the
  training set by 86 %, so validation loss rewards memorisation. A peptide-disjoint validation
  set would select for generalisation.
- ALC and local confidence here are mean / per-residue softmax probabilities, not PEAKS's
  candidate-agreement LC; the ALC >= 50 threshold therefore selects a different number of rows.
- `Accession` rule of PEAKS denovoOnly rows (2.5 % filled) is unknown; ours is an exact I=L
  substring match when a FASTA is given.
- Casanovo's `min_peptide_len` 6 excludes the 1.1 % of PEAKS de novo rows with length 4-5.
- Beam 10 was too slow in Casanovo 5.2.1 (> 30 min for 7k spectra on a shared GPU) and was not
  evaluated; beam 5 takes ~2.5 min per 7k spectra on a free A6000 (~40 min for 100k).
