# Status and known limitations (updated 2026-09-25 22:00 UTC)

Reference: DrewR_Adam_91526_DDA_OT_4 (Jurkat, Lumos, ion-trap HCD), PEAKS Studio DL boost on:
7,524 PSMs / 3,472 peptides / 990 protein groups.

| stage | state | evidence |
|---|---|---|
| export format round-trip | exact | 99 tests, all 4 reference dirs |
| protein inference / grouping / scores | exact except group numbering among ties | docs/PROTEIN_INFERENCE.md |
| featureless precursor correction | 92% agreement with PEAKS's corrected m/z | precursors_test |
| MS1 feature detection + chimera precursors | in progress | |
| candidate generation (Comet, instrument precursor) | 75% of PEAKS PSMs among top-5 candidates | comet_OT4 |
| rescoring (mokapot TDA) | 4,453 PSMs at 1% (recall of PEAKS 0.45) with misaligned MS2 feature; rerun pending | rescore_exp1 |
| rescoring (PEAKS-imitation, CV) | 5,471 PSMs at 1% (recall 0.55, precision 0.75), rerun pending | rescore_exp1 |
| de novo (Casanovo, high-res model) | 12% full-sequence accuracy on ion-trap data; unusable as is | denovo_casanovo |

De novo: Casanovo's released model is trained on Orbitrap MS2 and fails on ion-trap spectra
(2% agreement with PEAKS's de novo-only calls). Options: fine-tune Casanovo on the PEAKS-labelled
ion-trap PSMs once the HeLa PEAKS runs exist (tens of thousands of training spectra), or implement a
PEAKS-style dynamic-programming de novo for low-resolution data. Until then the denovoOnly export
is produced but not concordant.

## Session cut-off note (2026-09-25 ~22:15 UTC)
Usage limit reached mid-run. State at cut-off:
- Rescoring experiment 1 (instrument precursors, partly corrupted MS2 feature): mokapot 4,908 PSMs at 1%
  (recall 0.47); PEAKS-imitation CV 7,071 PSMs at 1% (recall 0.63, precision 0.67). Candidate ceiling 72%.
- Fixed: MS2PIP `correlate` misattaches vectors across PSMs sharing a spectrum; features are now
  computed from `predict_batch` plus our own matching (peaksparity/rescore_features.py).
- Experiment 2 launched on the corrected-precursor Comet candidates (ceiling 76% top-3):
  /drive_a/drew/peaks-parity-work/rescore_exp2/run.log, results.json when finished.
- Still running as background workers (write their own REPORT.md / repo files):
  features (/drive_a/drew/peaks-parity-work/features), monoiso tuning (monoiso_eval),
  PTM/modified-peptide tables (peaksparity/ptm_tables.py). Baseline Sage report is at
  /drive_a/drew/peaks-parity-work/baseline_sage/REPORT.md (unread).
Next steps: read the three worker reports and commit their files; wire features.py into
precursors.build_precursor_list (feature-derived entries) and rerun Comet on the expanded mzML
(expect candidate coverage >90%); rerun rescoring; run pipeline.py end-to-end with --ref for the
concordance report; ingest the Windows PEAKS exports (decoy-on) when Drew copies them.
