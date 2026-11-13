# proteoscan architecture

Goal: an open-source (MIT) DDA identification pipeline that runs on macOS and Linux and
reproduces PEAKS Studio's exports (PSMs, peptides, protein groups, features, de novo-only
spectra, PTM profile) as closely as measurable. PEAKS is closed source; we imitate its
behavior from its published algorithms and from its outputs on reference data. We never
inspect the PEAKS binary.

## Reference data

| set | instrument | runs | PEAKS search | role |
|---|---|---|---|---|
| DrewR_Adam_91526 | Lumos, OT MS1 / IT HCD MS2 | 4 (3 panel + 1 Jurkat) | human SP, CAM fixed, DL boost | primary concordance target |
| DrewR_Abate_081926 | same | 8 | BSA only | secondary |
| PXD004940 (HeLa, Lumos) | OT MS1 / IT HCD and OT HCD | 9 | to be run in PEAKS on Drew's desktop | large training/benchmark set |

## Stages

1. `convert`   Thermo .raw -> mzML (ThermoRawFileParser; see docs/INSTALL_MACOS.md for the macOS route).
2. `features`  MS1 feature detection: isotope envelopes, charge, monoisotopic m/z, RT
               range, area. Used to (a) correct precursor m/z and charge of each MS2,
               (b) split chimeric MS2 scans into one query per co-isolated feature
               (PEAKS: 262k MGF entries from 195k MS2 scans on Jurkat), (c) report
               "Feature Id" and "Area".
               Correction to (a)/(b), 2026-11-09: Comet takes the precursor m/z of each query
               from the Thermo trailer "Monoisotopic M/Z" userParam, which every expanded copy
               inherited unchanged, not from the selected-ion m/z Proteoscan writes. Before
               branch feat/comet-terminal-mods-precursor-offset every entry of a scan was
               therefore searched at the instrument's trailer m/z; the per-entry (feature-based
               or monoisotope-corrected) m/z changed only the charge tried. From that branch,
               `comet_precursor_mz_source: entry` (default) sets each copy's trailer to the
               entry m/z so the corrected m/z is searched; `trailer` reproduces the old
               behaviour. See docs/STATUS.md "Corrections (2026-11-09)".
               Every feature in the isolation window becomes an entry (`feature_linking: all`); each
               entry records whether it is the instrument's own precursor (`entry_primary`), its
               share of the co-isolated feature area and the number of entries of the scan, and
               the rescorer uses these. `feature_linking: primary` keeps only the instrument's
               precursor.
3. `denovo`    De novo sequencing per query (Casanovo or PEAKS-style tag graph) giving
               a full sequence, ALC, and per-residue local confidence. Provides the
               `denovoOnly` export and the "tag" features used by DB search rescoring.
4. `search`    Database search with Comet in the end-to-end pipeline (`proteoscan.pipeline`,
               `comet_search.py`: semi-tryptic, mods, wide isotope errors, several candidates
               per query, then a protein-restricted second pass). Sage is available only as
               the standalone `proteoscan.cli search` command (`sage_search.py`).
               Precursor calibration (`precursor_recalibration: auto`, `pipeline.search_run`): a
               quick estimation search (fully specific digest, one candidate per query, tolerance
               at least 20 ppm) gives the median and robust SD of the precursor error of confident
               PSMs; the real search removes an offset of 1 ppm or more and uses a tolerance of
               max(5 ppm, 10 robust SD), capped at the configured tolerance, in both passes. Values
               per run: `precursor_recalibration.json`.
5. `rescore`   Feature engineering per PSM in the pipeline: Comet scores, predicted fragment
               intensity correlation and observed b/y ion-series coverage from the predictor's
               annotation (MS2PIP by default; AlphaPeptDeep ion-trap models opt-in), predicted
               RT delta (DeepLC), sequence features and the second-pass protein-support feature;
               mokapot (semi-supervised target/decoy, cross-validated scores) or a supplied
               PEAKS-imitation model (`--mode imitate`); outputs a -10lgP-scaled score.
               Rescoring is deterministic: rows are put in a canonical order and both mokapot and
               its Percolator model are seeded (`mokapot.brew` otherwise builds an unseeded model
               whose inner cross-validation differs from call to call).
               De novo agreement (tag length, CAA%) and the extended ion features of
               `ion_features.py` are not used by the pipeline; they feed only the
               score-imitation training in `scripts/train_score_model.py`.
6. `fdr`       Decoy-fusion-style target/decoy FDR at PSM, peptide sequence and protein
               group level, thresholds chosen exactly the way PEAKS reports them
               (peptide -10lgP threshold at 1% peptide FDR; protein -10lgP >= 15 and
               >= 1 unique peptide). `run_stats.json` says whether the FDR is estimable
               (`fdr_estimable`, decoy counts at the threshold); for small defined mixtures
               see README "Small defined mixtures".
7. `proteins`  Protein score, grouping and "Top" flag following PEAKS rules.
   Cached stages in an output directory carry a fingerprint of the config, inputs, database and
   options (`cache_fingerprint.json`); a mismatch is an error unless `--overwrite` is given.
8. `export`    Write the PEAKS CSV set with identical columns (`proteoscan.peaks_export`).
9. `concord`   Compare any export directory to a PEAKS reference; the optimization loop
               maximizes PSM recall/precision at exact modified-sequence level.
10. `ui`       Browser-based results viewer (later).

## Optimization loop

`scripts/optimize.py` runs stages 4-8 with a config, scores concordance, and logs to
`runs/<id>/`. Subagents propose config/model changes; each run is accepted only if the
concordance harness improves on held-out runs.
