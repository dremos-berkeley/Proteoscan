# Status and known limitations (updated 2026-09-30)

Reference: DrewR_Adam_91526_DDA_OT_4 (Jurkat, Lumos, ion-trap HCD), PEAKS Studio, DL boost on:
7,524 PSM rows / 3,472 peptides / 990 protein groups. Concordance = exact modified sequence,
scan and charge (PSM); exact modified sequence (peptide); shared accession (protein group).

## Exactly reproduced from PEAKS's own tables
| item | state |
|---|---|
| CSV formats (all 8 tables, Java number formatting) | byte-exact round trip |
| protein inference, grouping, scores, coverage, areas | exact; group numbering among ties differs |
| PTM profile | byte-identical |
| modified-peptide table | all cells except 8 RT-shift values (0.0001 off) |
| precursor list (given PEAKS's features) | 98.8% of the precursors behind PEAKS PSMs |
| I/L variant rows and protein lists per PSM | reproduced |

## End-to-end from raw data (Comet + rescoring), Jurkat run
| configuration | PSMs reported | PSM recall | PSM precision | peptide recall | protein-group recall |
|---|---|---|---|---|---|
| Comet e-value only | 1,198 | 0.13 | 0.84 | | |
| + rescoring (mokapot), corrected precursors | 11,126 | 0.74 | 0.50 | 0.82 | 0.95 |
| + PEAKS-imitation model (cross-validated) | 10,722 | 0.75 | 0.52 | 0.83 | 0.95 |
| + two-pass search with protein support, imitation | 14,609 | 0.81 | 0.41 | 0.88 | 0.95 |
| same, cut at PEAKS's PSM count | 7,595 | 0.67 | 0.66 | 0.78 | 0.88 |

On scans identified by both, the peptide is identical in 97% of cases.

Entrapment validation (human + shuffled entrapment database of equal size): measured false
discovery proportion 0.5% at the PEAKS-style filter and 0.9% at nominal 1% PSM FDR (single
pass). The pipeline therefore reports more correct identifications than PEAKS at the same
nominal error rate; "precision" against PEAKS is low because PEAKS reports fewer PSMs, not
because the extra PSMs are wrong. Two-pass entrapment validation: pending.

## Findings that shaped the design
- Feature detection is not the bottleneck: our detector (70% recall of PEAKS features) and
  PEAKS's own feature list give the same candidate coverage (78%).
- The bottleneck is candidate ranking on weak/chimeric spectra. PEAKS rescues these with its
  protein feature; the two-pass search imitates that (coverage 83%).
- Restricting to zero isotope error loses recall (0.67 vs 0.75): keep Comet's -1..3 window.
- Sage needs ~300-500 GB for semi-tryptic + 3 variable mods; Comet needs 3 GB. Comet is the engine.
- MS2PIP 4.2 `correlate` misattaches vectors across PSMs sharing a spectrum; features are computed
  from `predict_batch` plus our own peak matching.

## Open items
- PEAKS unfiltered export with decoys (running on the Windows machine): gives PEAKS's score for
  every PSM, needed to model its score and threshold rather than only its accepted set.
- De novo for ion-trap spectra: Casanovo's released model reaches 12% full-sequence accuracy on
  this data; needs fine-tuning on ion-trap PSMs. The denovoOnly export is not yet concordant.
- features.csv export and Area columns from our own feature detector.
- HeLa benchmark (PXD004940) runs through the pipeline and through PEAKS.
- macOS packaging and the browser results viewer (in progress).

## HeLa benchmark (PXD004940, 90-min gradients, 1% FDR, pipeline defaults; 2026-09-30)
| run | MS2 | PSMs | peptides | protein groups | wall min |
|---|---|---|---|---|---|
| 03_01 500 ng | ion trap | 47,030 | 38,138 | 5,172 | 35 (single pass) |
| 03_02 500 ng | ion trap | 46,728 | 37,764 | 5,239 | 35 (single pass) |
| 03_03 500 ng | ion trap | 49,123 | 40,092 | 4,979 | 58 |
| 05_01 1000 ng | ion trap | 55,767 | 44,565 | 5,149 | 63 |
| 05_02 1000 ng | ion trap | 54,572 | 42,966 | 5,166 | 57 |
| 05_03 1000 ng | ion trap | 55,563 | 44,108 | 5,184 | 58 |
| 03_01 1 ug | Orbitrap 0.02 Da | 30,958 | 27,543 | 3,992 | 32 |
| 03_02 1 ug | Orbitrap 0.02 Da | 32,111 | 28,694 | 4,029 | 26 |
| 03_03 1 ug | Orbitrap 0.02 Da | 32,145 | 28,633 | 4,024 | 38 |
Exports under /drive_a/drew/peaks-parity-work/hela/<run>/export. PEAKS results for these runs are
pending on the Windows machine; they will be the training set for the score-imitation model.

## De novo (2026-09-30)
Casanovo fine-tuned on the HeLa ion-trap PSMs (docs/DENOVO.md): held-out HeLa accuracy 42%
(released model 17%); Jurkat 20% vs 18% (different scan mode: turbo vs rapid). PEAKS's own de novo
agrees with database identifications on only 4% of its de novo-only spectra on this data, so
sequence-level concordance of the denovoOnly export is not a meaningful target on ion-trap data;
the export is produced in the PEAKS layout with our model's calls and confidences.
De novo agreement features (docs/DENOVO_FEATURES.md): CAA / tag / beam features from Casanovo calls on every
(scan, precursor) query raise the PEAKS-score imitation from CV Pearson 0.620 (base + ion features) to 0.631
(score_prop_v2). The same features from PEAKS's own de novo reach 0.706. PEAKS's exported `CAA (%)` is the
share of residues with Positional Confidence >= 2 %, and `Tag Length` = `Length`; neither is the de novo CAA.

## Score imitation (2026-09-30)
Trained on the PEAKS propionamide export (decoys included, scores down to 15), cross-validated by
scan parity, evaluated by reproducing PEAKS's own accepted set (its 1% sequence-FDR threshold):
| model | CV Pearson with PEAKS -10lgP | Jurkat at PEAKS threshold: recall / precision |
|---|---|---|
| base features | 0.61 | 0.55 / 0.79 |
| + neutral-loss and 2+ ion evidence | 0.62 | 0.56 / 0.81 |
| + de novo agreement (fine-tuned Casanovo) | 0.63 | 0.60 / 0.84 (all 4 files) |
| ceiling: features from PEAKS's own de novo | 0.71 | 0.61 / 0.86 |
Reference point: two PEAKS runs of the same file agree with each other at Pearson 0.53-0.78 and
share 67-90% of accepted PSMs, depending on settings. Setting PEAKS's report filters to zero
changes its de novo output (CAA agreement 1.00 -> 0.67) and its scores; default filters plus
export-decoy is the right reference export.

## Ground-truth accuracy: ProteomeTools synthetic pools (2026-10-02)
24 pools of PRIDE PXD004732 (about 1,086 synthetic peptides each; truth = Supplementary Table 1 of
Zolg et al. 2017 plus the 87 spiked QC peptides), ion-trap runs that interleave HCD-IT, CID-IT and
HCD-OT scans of every precursor. Both methods searched the full human SwissProt database with the
same parameters; PEAKS Studio 13.5 with DL boost on. Means per pool over 23 pools (pool 123, a sparse
acquisition, excluded); "at 1%" = each method's own accepted list truncated at a measured 1% error.
| scan type | method | PSMs | measured FDP | true peptides | true peptides at 1% FDP |
|---|---|---|---|---|---|
| HCD ion trap | ours | 8,910 | 1.14% | 815 | 707 |
| HCD ion trap | PEAKS | 8,796 | 1.31% | 831 | 706 |
| CID ion trap | ours | 8,535 | 1.13% | 805 | 728 |
| CID ion trap | PEAKS | 9,056 | 1.34% | 830 | 720 |
| HCD Orbitrap | ours | 18,863 | 0.90% | 848 | 748 |
| HCD Orbitrap | PEAKS | 18,049 | 0.94% | 841 | 781 |
Conclusion: parity in sensitivity and accuracy with PEAKS + DL boost on synthetic ground truth.
Full tables: /drive_a/drew/peaks-parity-work/proteometools/SUMMARY_truth.{md,csv}; script
scripts/aggregate_proteometools.py.

## Ion-trap fragment-intensity model (2026-10-03, docs/IONTRAP_MODEL.md)
peptdeep MS2 models fine-tuned on 52 further ProteomeTools pools (HCD-IT and CID-IT; validation on 8
more, benchmark pools never used). Benchmark (23 pools, 3 mokapot seeds): HCD-IT true peptides 803 -> 838
(19 of 23 pools up), CID-IT 805 -> 818 (20 up), measured FDP unchanged within noise. Jurkat: no gain (the
model predicts turbo-scan Jurkat spectra no better than MS2PIP). Opt-in via `--fragment-model peptdeep:...`;
MS2PIP HCD2021 remains the default.

## PEAKS Deep Learning Boost on vs off, ProteomeTools (2026-10-03)
Means over 23 pools, each method at its own claimed 1% FDR:
| scan type | PEAKS boost off: PSMs / FDP / true peptides | PEAKS boost on | Proteoscan |
|---|---|---|---|
| HCD ion trap | 7,659 / 0.98% / 800 | 8,796 / 1.31% / 831 | 8,910 / 1.14% / 815 (fine-tuned model: 9,396 / 1.19% / 838) |
| CID ion trap | 7,899 / 0.96% / 798 | 9,056 / 1.34% / 830 | 8,535 / 1.13% / 805 (fine-tuned: 9,054 / 1.24% / 818) |
| HCD Orbitrap | 14,715 / 0.55% / 798 | 18,049 / 0.94% / 841 | 18,863 / 0.90% / 848 |
The boost adds about 4% true peptides on these clean synthetic spectra (82% of PSMs on the Jurkat
lysate); Proteoscan's rescoring matches or exceeds the boosted result and exceeds the unboosted one
by 2-6%. Full table: /drive_a/drew/peaks-parity-work/proteometools/SUMMARY_truth.csv.

## Corrections (2026-10-05, branch fix/comet-fixed-mods-score-tail-cv)
Three problems found in the 261003_proteoscan_test replication (propionamide config) and fixed:
- **P1, hidden static +57 on Cys (bug).** `build_comet_params` started from `comet -p` defaults
  (`add_C_cysteine = 57.021464`, `variable_mod01 = 15.9949 M`) and never zeroed them. Configs without
  explicit fixed CAM (`configs/peaks_propionamide.yml`) searched every Cys as C+57 plus the variable
  CAM/propionamide stack (C+114 / C+128): propionamide and free Cys were unsearchable, exported
  `C(+57.02)` / `C(+71.04)` were really +114 / +128, and mokapot learned negative `has_cys` weights.
  Fixed: all `add_*` and variable slots are zeroed first; `fixed_mod_deltas`; a variable mod on a
  fixed residue is refused. Fixed-CAM configs (facility, HeLa, ProteomeTools) were not affected:
  their comet.params differ from before only in already-zero lines, and Comet output is identical.
- **P2, -10lgP ceiling on small runs (design limit).** `minus10lgp` tied every target above the best
  decoy at 10*log10(n_decoys+1), and the FDR threshold could rest on 0 decoys with no warning. Now
  the decoy tail is extrapolated above the best decoy, and `run_stats.json` reports decoy counts at
  the threshold, (D+1)/T, the decoy null size and ceiling, and `fdr_estimable` (>= 10 decoy
  sequences at the threshold). README "Small defined mixtures".
- **P3, in-sample rescoring scores (methodological).** `rescore_tda` scored every PSM with the mean
  of mokapot's fold models, i.e. with models trained on it. It now uses mokapot's cross-validated
  (held-out, per-fold calibrated) scores, checked against mokapot's own; adds `mokapot_q`.
  `has_cys` is no longer a rescoring feature.

Superseded by these corrections: the "Score imitation (2026-09-30)" table above and
docs/DENOVO_FEATURES.md numbers (trained on propionamide candidates with the hidden +57; models
`score_prop_v1`, `v1b_ions`, `v2` invalidated, models/README.md). The HeLa, Jurkat and
ProteomeTools tables above were produced with the in-sample scores of P3.

ProteomeTools regression with the corrections (same 23 pools, full rerun, mokapot seed 0; same
configs as the baseline: configs/peaks_facility_human.yml for the ion-trap scan types,
configs/peaks_hela_hcd_ft.yml (0.02 Da) for HCD Orbitrap; /drive_a/drew/peaks-parity-work/fixcys_regress/
SUMMARY_truth.{md,csv} and SUMMARY_truth_ot002.{md,csv}, before_after*.csv). Means over pools, before -> after;
per-pool up / down. Noise: SD over three mokapot seeds of the 23-pool mean is 4.95 (HCD-IT) and 2.89 (CID-IT)
true peptides at 1% FDP and about 0.02 percentage points of FDP (iontrap_model/results/bench_seeds.csv); for HCD
Orbitrap, from the rescoring ablation below, 1.1-1.5 true peptides, 0.7-1.2 true peptides at 1% FDP and
0.01 percentage points of FDP.
| scan type | PSMs | measured FDP | true peptides | true peptides at measured 1% FDP |
|---|---|---|---|---|
| HCD ion trap | 8,910 -> 8,851 | 1.14% -> 1.11% | 815 -> 818 (9 up / 10 down) | 707 -> 738 (median change 0; 11 up / 9 down) |
| CID ion trap | 8,535 -> 8,365 | 1.13% -> 1.13% | 805 -> 801 (8 up / 11 down) | 728 -> 744 (median 0; 10 up / 10 down) |
| HCD Orbitrap | 18,863 -> 18,684 | 0.90% -> 0.88% | 848 -> 847 (8 up / 11 down; median 0) | 748 -> 787 (median 0; 10 up / 10 down) |
No regression on any scan type: measured FDP is unchanged or lower, true peptides are within noise,
and fewer PSMs with unchanged true peptides is the expected effect of removing in-sample scores (P3).
The mean gains in "true peptides at measured 1% FDP" come from a few pools (HCD-OT +589 and +328)
whose baseline list was truncated by an early false PSM tied at the old -10lgP ceiling (see the P2
finding below); the medians are 0.

An earlier version of this regression searched HCD Orbitrap with the facility config (0.5 Da, low-res
binning) instead of the baseline's 0.02 Da config and showed 848 -> 808 true peptides. That was the
configuration change, not the corrections: a rescoring-only ablation on cached candidate tables
(no re-search; 23 pools x 3 mokapot seeds; /drive_a/drew/peaks-parity-work/ot_ablation/arm_table.csv)
gives, as true peptides (SD over seeds) / measured FDP:
| candidates | master scoring (in-sample, with has_cys) | in-sample, no has_cys | cross-validated, with has_cys | cross-validated, no has_cys (branch) |
|---|---|---|---|---|
| baseline HCD-OT, 0.02 Da | 847.0 (1.5) / 0.90% | - | - | 845.5 (1.1) / 0.85% |
| HCD-OT searched at 0.5 Da | 812.4 (1.2) / 0.74% | 810.2 (5.2) / 0.75% | 809.8 (7.6) / 0.74% | 811.0 (2.0) / 0.72% |
Neither the cross-validated scores nor dropping has_cys explains a loss; the branch scoring costs
about 1.5 true peptides on matched candidates (about one seed SD), consistent with removing
in-sample optimism, and raises the mean HCD-OT threshold only from 26.53 to 26.65. On 8 HCD-IT pools
the four scoring arms do not differ systematically (pool 1 varies by up to 220 true peptides between
seeds in every arm).

P2 finding: re-filtering identical scores with the old flat ceiling instead of the extrapolated tail
leaves PSMs, true peptides and FDP unchanged in all 23 HCD-OT pools, but changes "true peptides at
measured 1% FDP" in 2 pools (+3, -307): with the flat ceiling every PSM above the best decoy ties,
so the score-ordered list is ordered arbitrarily among them and one early false PSM can truncate it.
This tie-breaking artefact is the source of the large single-pool swings in that metric in earlier
tables (STATUS rows before 2026-10-05, docs/IONTRAP_MODEL.md).

Also observed: mokapot brew with max_workers > 1 is not reproducible for a fixed seed (the fold
models share one generator across threads); even with max_workers = 1, repeated runs differed on
one pool (9,435 vs 9,073 PSMs). Single-seed comparisons carry that noise.
