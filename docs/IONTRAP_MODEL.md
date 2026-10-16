# Ion-trap fragment-intensity model (2026-10-03)

Question: does an AlphaPeptDeep (peptdeep) MS2 model fine-tuned on ion-trap spectra identify more
than the MS2PIP HCD2021 (Orbitrap HCD) model the pipeline uses for its fragment features?

Short answer: on the ProteomeTools ground-truth benchmark, yes for ion-trap HCD (+35 true peptides
per pool, +4.4%, 19 of 23 pools up, measured error unchanged) and modestly for ion-trap CID (+13 true
peptides, +1.6%, 20 of 23 up); on the Jurkat run (a different instrument and scan rate) no: the
fine-tuned model predicts those spectra no better than MS2PIP and returns 3% fewer PSMs at the same
recall of the PEAKS export. MS2PIP HCD2021 stays the default; the model is opt-in
(`--fragment-model peptdeep:...`).

## Test-set hygiene
- Never used for training, checkpoint choice or feature choice: the 24 benchmark pools of
  /drive_a/drew/public-data/PXD004732/pools.txt and DrewR_Adam_91526_DDA_OT_4 (Jurkat).
  `scripts/iontrap_prepare.sh`, `iontrap_val_pipeline.sh`, `iontrap_dataset.py` (without
  `--allow-benchmark`) and `iontrap_train.py` refuse benchmark pools.
- Training spectra of any peptide that is in a benchmark pool (or is one of the 87 QC peptides
  spiked into every pool) are dropped before training (3,642 ProteomeTools and about 17,000 HeLa
  spectra), so benchmark peptides are unseen sequences for the model.
- Model and feature choices were made on 8 held-out validation pools (TUM_first_pool_125, 83, 122,
  66, 41, 47, 86, 4) and, for the complex-sample check, on HeLa peptides held out by sequence hash.
  The benchmark and Jurkat were run once with the pre-selected configuration (plus the released
  peptdeep model as a reference row).

## Data
| source | use | activation, NCE, scan rate (filter string) | spectra | peptides |
|---|---|---|---|---|
| ProteomeTools PXD004732, 52 TUM_first_pool pools (2xIT_2xHCD-1h) | train | HCD-IT 28, rapid / CID-IT 35, rapid | 61,231 / 60,947 | 44,563 / 44,487 |
| same, 8 other pools | validation (early stopping, selection) | as above | 8,562 / 8,535 | 6,301 / 6,300 |
| HeLa PXD004940, 6 ion-trap runs, our PSMs with -10LgP >= 30, best per (peptide, charge) | train (mixed model only); 10% of peptides held out by hash | HCD-IT 28, rapid | 236,456 | 55,559 |
| 24 benchmark pools (labelled the same way) | test only | HCD-IT 28 / CID-IT 35, rapid | 25,856 / 25,767 | 19,492 |
| Jurkat, PEAKS facility export PSMs (7,478 of 7,524 parsed) | test only | HCD-IT 30, **turbo** | 7,478 | 3,283 |
The 2xIT_2xHCD method acquires every precursor as a block of four consecutive MS2 scans: CID-IT
(NCE 35), HCD-IT (NCE 28), HCD-OT (NCE 20) and HCD-OT (NCE 23). PXD045495 (Eclipse) was not used.

Downloads: 60 raw files (54 GB) from the PRIDE FTP over HTTPS into
/drive_a/drew/public-data/PXD004732/raw_train, every file SHA1-checked against PRIDE (8 files corrupted
by an early downloader race were deleted with everything derived from them and fetched again).
Conversion `ThermoRawFileParser -f 2` (vendor centroiding), split with `proteoscan.split_mzml`.

### Labelling the ProteomeTools spectra (scripts/iontrap_label.py)
Targeted matching against the pool's own peptide list (Supplementary Table 1 + QC peptides; C
carbamidomethyl, which the data confirm: Cys peptides are labelled at the same rate as others; each
single M oxidised as an extra form) at charge 1-4, 15 ppm, isotope offset 0/+1. Score = binomial
peak-match score of b/y (1+, plus 2+ for 3+ precursors) against the top 10 peaks per 100 m/z;
decoy = reversed peptide (same precursor mass) scored on the same scans; 1% threshold on the best
(peptide form, charge) per pool; one peptide per scan; >= 6 matched b/y ions in the labelled scan.

Deviation from the brief: scoring only the ion-trap HCD scan is too weak in some runs (dense
ion-trap spectra of short 2+ peptides; e.g. pool 4: 269 HCD-IT labels from IT-only scoring vs 735
CID-IT). The labeller therefore scores the whole four-scan block (sum of per-scan scores, Orbitrap
scans at 0.02 Da, ion-trap scans at 0.5 Da), applies the decoy threshold to the block, and labels the
block's HCD-IT and CID-IT scans that themselves match >= 6 b/y ions. Per pool: 1,165 accepted blocks
(913 peptides), 1,163 HCD-IT and 1,158 CID-IT spectra on average. The IT-only mode is kept
(`--split-dir`) and was used for comparison.

Validation against the pipeline's own PSMs on the 8 validation pools (pipeline run exactly as for
the benchmark; PSMs whose peptide is a pool member with -10LgP >= 25), `scripts/iontrap_label_check.py`:
| labels | scan type | labels per run | scans with both | same peptide (I=L) and charge | same modified peptide | pipeline (peptide, charge) pairs with a label | labelled scans contradicted by any pipeline PSM |
|---|---|---|---|---|---|---|---|
| block (used) | HCD-IT | 1,070 | 877 | 99.76% | 86.1% | 94.1% | 0.47% |
| block (used) | CID-IT | 1,067 | 814 | 99.76% | 89.8% | 94.5% | 0.63% |
| IT-only | HCD-IT | 890 | 743 | 99.90% | 87.1% | 78.7% | 0.36% |
| IT-only | CID-IT | 1,020 | 776 | 99.70% | 90.6% | 95.2% | 0.78% |
Of the modified-sequence disagreements, 92% (1,518 / 1,642) are the pipeline reporting an extra
deamidation (N/Q +0.98) with an isotope-error shift (precursor error about +11 ppm) where the label is
the unmodified peptide at about 0 ppm. On synthetic unmodified peptides the label is the plausible
one; this looks like a pipeline weakness (deamidation + isotope error competing with the correct
monoisotopic assignment) worth its own look.

## Model and training (scripts/iontrap_train.py)
- Base: peptdeep 1.5.1 released generic MS2 model (pretrained_models_v3, ModelMS2Bert). Every layer
  is loaded; the output head is widened from b/y 1+/2+ to b/y 1+/2+ plus b/y -H2O and -NH3 (1+), the
  new rows initialised from the parent ion's row (peptdeep supports these fragment types; its
  modification-loss head is dropped).
- Inputs: instrument "Lumos", the spectrum's own NCE from the filter string, charge, sequence and mods.
  peptdeep has no activation input, so HCD-IT and CID-IT are separate models (no joint model).
- Targets: max peak within 0.4 Da of each ion, normalised to the spectrum's highest annotated ion;
  ions that cannot be observed (fragment charge > precursor charge, m/z outside the scan window,
  e.g. the CID low-mass cut-off) are masked.
- Loss: masked spectral angle. AdamW, lr 1e-4, weight decay 1e-4, batches of 256 same-length peptides,
  first epoch output head only with linear warm-up, cosine decay over 30 epochs, early stopping
  (patience 4) on the median spectral angle of the validation pools (averaged with the held-out HeLa
  peptides for the mixed model). Seed 0. About 5 s per epoch on one A6000.
- Checkpoints (not in git): /drive_a/drew/proteoscan-work/models/iontrap_model/models/
  | model | training data | spectra | best epoch | val median SA (b/y 1+) |
  |---|---|---|---|---|
  | hcd_it_pt.pth | ProteomeTools HCD-IT | 57,589 | 30 of 30 | 0.934 |
  | cid_it_pt.pth | ProteomeTools CID-IT | 57,305 | 29 of 30 | 0.911 |
  | hcd_it_mix.pth | ProteomeTools HCD-IT + HeLa (60,000 HeLa per epoch) | 255,802 | 15 of 18 | 0.872 (HeLa held-out 0.867) |
  Each .pth has a .json (fragment types, NCE, data, settings) and .history.csv.

## Predictor accuracy (scripts/iontrap_eval_predictor.py)
Median per-spectrum spectral angle (1 - 2*angle/pi) and Pearson r on max-normalised linear
intensities of the observable b/y 1+ ions (the ions every predictor predicts); each spectrum predicted
at its own NCE.
| dataset | spectra | MS2PIP HCD2021 | MS2PIP CID | peptdeep released | fine-tuned PT | fine-tuned mixed |
|---|---|---|---|---|---|---|
| PT validation pools, HCD-IT | 8,562 | 0.745 / 0.885 | | 0.871 / 0.971 | **0.937 / 0.993** | 0.887 / 0.976 |
| PT benchmark pools (23), HCD-IT | 25,710 | 0.754 / 0.891 | | 0.902 / 0.983 | **0.947 / 0.995** | 0.917 / 0.988 |
| PT validation pools, CID-IT | 8,535 | 0.553 / 0.641 | 0.767 / 0.896 | 0.526 / 0.592 | **0.915 / 0.987** | |
| PT benchmark pools (23), CID-IT | 25,621 | 0.587 / 0.686 | 0.788 / 0.915 | 0.561 / 0.648 | **0.931 / 0.991** | |
| HeLa 05_03, held-out peptides | 4,128 | 0.623 / 0.745 | | 0.788 / 0.918 | 0.711 / 0.847 | **0.863 / 0.965** |
| Jurkat PEAKS PSMs (turbo, NCE 30) | 7,478 | **0.717 / 0.852** | | 0.687 / 0.819 | 0.687 / 0.814 | 0.619 / 0.723 |
Including all predicted fragment types (2+ and neutral losses) the fine-tuned PT models reach median
SA 0.880 (HCD) and 0.885 (CID) on the benchmark pools. On Jurkat the ranking does not depend on the
NCE given to the model (0.681 at NCE 28) or on PSM score quartile (MS2PIP best in every quartile).
The three sample types disagree on which model is best: each model is best on the data it resembles
(synthetic peptides on the 2016 Lumos; HeLa on the PXD004940 instrument; neither for Jurkat).

## Integration (proteoscan/fragpred.py, proteoscan/rescore_features.py)
`compute_features(..., ms2pip_model="peptdeep:<spec>")` computes the same columns as the MS2PIP path
(ms2_corr, ms2_cos, ms2_frac_pred_found, longest_b_obs, longest_y_obs, n_ions_obs, residue_coverage):
observed log2(intensity/TIC + 0.001) against the prediction renormalised to sum 1 over the same
ions, log2(x + 0.001). Spec: `peptdeep:<ckpt.pth>`, `peptdeep:released`, or per activation
`peptdeep:hcd=<ckpt>,cid=<ckpt>`; suffix `#all` uses every predicted, observable fragment type instead
of b/y 1+. Each PSM is predicted at its scan's NCE. Batched in chunks of 50,000 PSMs (bounded memory),
eval mode; repeated calls give identical values, changing the chunking changes values by < 3e-6.
The ion-series columns agree with the MS2PIP path to rounding (correlation 0.9997-0.9998).
Pipeline flag: `python -m proteoscan.pipeline ... --fragment-model peptdeep:hcd=<ckpt>,cid=<ckpt>`.
Default unchanged (MS2PIP HCD2021). Tests: tests/test_fragpred.py.

## Choice made on the validation pools (scripts/iontrap_benchmark.py on val_pipeline runs)
True peptides at measured 1% FDP, mean over the 8 validation pools:
| scan type | MS2PIP | peptdeep released | fine-tuned PT | fine-tuned mixed (HCD) |
|---|---|---|---|---|
| HCD-IT | 712.5 | 740.9 | **741.3** | 737.6 |
| CID-IT | 702.3 | 701.0 | **726.1** | (= PT) |
On 5 validation pools with earlier checkpoints, b/y 1+ features beat `#all` (HCD +42 vs +3 true
peptides over MS2PIP; CID +42 vs -4). Selected: fine-tuned PT models, b/y 1+ features.

## Benchmark: 24 ProteomeTools pools (pool 123 excluded from means, as in STATUS.md)
Cached candidate tables (`<run>/<type>/candidates_features_2pass.parquet`), no re-search; the seven
fragment columns replaced by the predictor's values; then `proteoscan.pipeline.run` on the cached
tables (rescore_tda, peaks_filter, export) and scoring of the export exactly as
scripts/aggregate_proteometools.py does (truth = Supplementary Table 1 + QC peptides). mokapot's
result depends on its random seed (seed-to-seed SD of true peptides up to 58 in pool 1), so every
configuration was rescored with seeds 0, 1, 2; numbers are means over seeds, then over 23 pools.
| scan type | predictor | PSMs | measured FDP | true peptides | PSMs at measured 1% FDP | true peptides at measured 1% FDP |
|---|---|---|---|---|---|---|
| HCD-IT | MS2PIP HCD2021 | 8,811 | 1.16% | 803 | 7,665 | 695 |
| HCD-IT | fine-tuned (hcd_it_pt) | 9,396 | 1.19% | **838** | 8,222 | **732** |
| CID-IT | MS2PIP HCD2021 | 8,547 | 1.17% | 805 | 7,560 | 725 |
| CID-IT | fine-tuned (cid_it_pt) | 9,054 | 1.24% | **818** | 7,938 | 719 |
Paired over pools (fine-tuned minus MS2PIP, seed means; Wilcoxon signed-rank):
| scan type | true peptides | true peptides at measured 1% FDP | PSMs |
|---|---|---|---|
| HCD-IT | +35.6 mean, +4 median, 19 up / 3 down, p = 0.001 | +37.4 mean, +4 median, 18 up / 4 down, p = 0.003 | +585, 21 up / 2 down |
| CID-IT | +13.3 mean, +8.3 median, 20 up / 3 down, p = 0.0009 | -5.7 mean, +8 median, 19 up / 3 down, p = 0.007 | +507, 23 up / 0 down |
- The HCD-IT mean is carried partly by three pools where the MS2PIP rescoring is weak (pool 1: 332
  vs 659 true peptides, pools 6 and 10: about +100); without pool 1 the gain is +20 true peptides.
- The CID-IT "at measured 1% FDP" mean is negative only because of pool 114, where one early false
  PSM truncates the list (133 vs 601); its median and sign test favour the fine-tuned model. This
  metric collapses to near 0 for pools 12 and 81 under both predictors (an early high-scoring false
  PSM), so the plain "true peptides" column is the steadier comparison.
- Measured FDP is unchanged within noise (HCD 1.16% vs 1.19%, CID 1.17% vs 1.24%; per-pool sign tests
  not significant). Both predictors overshoot 1% slightly, as in STATUS.md.
- The released peptdeep model (seed 0 only) matches the fine-tuned model on HCD-IT (true peptides at
  1% FDP 731 vs 732) but not on CID-IT (723, = MS2PIP): fine-tuning matters for CID, while for HCD-IT
  most of the gain comes from replacing MS2PIP by peptdeep.
- Rerunning the MS2PIP baseline through the same code reproduces the STATUS.md exports exactly for 16
  of 48 runs and within a few PSMs for most others; pool 1 HCD-IT differs strongly (525 true peptides
  originally, 332 / 307 / 221 in the three reruns). The comparisons above are all from reruns.

## Jurkat (DrewR_Adam_91526_DDA_OT_4, cached candidates, concordance with the facility PEAKS export)
Pipeline export at the PEAKS-style 1% filter vs the facility PEAKS export (7,478 PSMs, 3,472 peptides,
990 protein groups), `proteoscan.concordance`; means over mokapot seeds 0-2 (SD in brackets).
| predictor | PSMs at 1% | peptides | protein groups | PSM recall | PSM precision | peptide recall | protein-group recall |
|---|---|---|---|---|---|---|---|
| MS2PIP HCD2021 (default) | 15,425 (34) | 5,972 | 1,271 | **0.790** (0.002) | 0.383 | **0.865** | **0.927** |
| peptdeep released | 15,211 (189) | 5,874 | 1,251 | 0.787 | 0.387 | 0.862 | 0.921 |
| fine-tuned PT (hcd_it_pt, pre-selected) | 14,971 (28) | 5,757 | 1,238 | 0.787 (0.002) | 0.393 | 0.860 | 0.920 |
| fine-tuned mixed (hcd_it_mix) | 14,913 (61) | 5,693 | 1,214 | 0.785 | 0.394 | 0.857 | 0.911 |
The fine-tuned model reports 3% fewer PSMs with PSM recall of the PEAKS set unchanged within seed noise
(0.787 vs 0.790) and slightly higher precision; it does not help on this run, consistent with its
spectral angle on Jurkat (no better than MS2PIP). The rerun MS2PIP baseline (15,425) differs from the
cached STATUS.md run (15,398 PSMs, recall 0.790) by mokapot run-to-run variation.

## Files
- Code: proteoscan/fragpred.py (peptide parsing, fragment m/z, annotation, PeptdeepPredictor),
  proteoscan/rescore_features.py (peptdeep_correlations, SpectrumStore scan info), proteoscan/pipeline.py
  (`--fragment-model`), scripts/iontrap_{label,label_check,dataset,train,eval_predictor,benchmark,jurkat}.py,
  scripts/iontrap_{prepare,val_pipeline,build_pt}.sh.
- Work directory /drive_a/drew/proteoscan-work/models/iontrap_model: runs/ (converted, split, labelled
  training pools), data/ and test_data/ (fragment tables), models/, val_pipeline/ (pipeline runs of the
  validation pools), val_bench/, bench/, jurkat_bench/ (rescored exports per variant and seed),
  results/ (label_check_*.csv, predictor/*.csv with per-spectrum parquet files, val_bench.csv,
  bench.csv, bench_seeds.csv).

## What remains uncertain
- Transfer to other instruments: the model that wins on the synthetic pools does not win on Jurkat
  (turbo scan, NCE 30, 2026 Lumos) or HeLa (PXD004940); intensity patterns of ion-trap spectra appear
  instrument- and method-specific. A model for Drew's own data would need training PSMs from the same
  instrument and scan rate; PEAKS PSMs of other facility runs (not run 4) would be the natural source.
- The benchmark gain is concentrated in pools where the MS2PIP-based rescoring is weak; median
  per-pool gains are small (+4 HCD, +8 CID true peptides). mokapot seed variance is of the same order
  as the gain in single pools; three seeds per configuration were used.
- Labels are the best-matching scan per precursor block, i.e. biased towards clean spectra; models
  were never trained on weak or chimeric ProteomeTools spectra.
- peptdeep has no activation input: one checkpoint per activation; a joint model was not tried.
- Neutral-loss and 2+ outputs are trained but did not help as rescoring features (`#all`).
- PXD045495 (Eclipse) was not processed.
