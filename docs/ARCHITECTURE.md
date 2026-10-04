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
3. `denovo`    De novo sequencing per query (Casanovo or PEAKS-style tag graph) giving
               a full sequence, ALC, and per-residue local confidence. Provides the
               `denovoOnly` export and the "tag" features used by DB search rescoring.
4. `search`    Database search with Sage (semi-tryptic, mods, wide isotope errors,
               multiple candidates per query). Optionally tag-filtered candidate lists.
5. `rescore`   Feature engineering per PSM: Sage scores, predicted fragment intensity
               correlation (AlphaPeptDeep / MS2PIP, ion-trap models), predicted RT delta
               (DeepLC), de novo agreement (tag length, CAA%), positional confidence.
               A learned model (LightGBM / mokapot) trained to imitate PEAKS acceptance
               and rank order; outputs a -10lgP-scaled score.
6. `fdr`       Decoy-fusion-style target/decoy FDR at PSM, peptide sequence and protein
               group level, thresholds chosen exactly the way PEAKS reports them
               (peptide -10lgP threshold at 1% peptide FDR; protein -10lgP >= 15 and
               >= 1 unique peptide). `run_stats.json` says whether the FDR is estimable
               (`fdr_estimable`, decoy counts at the threshold); for small defined mixtures
               see README "Small defined mixtures".
7. `proteins`  Protein score, grouping and "Top" flag following PEAKS rules.
8. `export`    Write the PEAKS CSV set with identical columns (`proteoscan.peaks_export`).
9. `concord`   Compare any export directory to a PEAKS reference; the optimization loop
               maximizes PSM recall/precision at exact modified-sequence level.
10. `ui`       Browser-based results viewer (later).

## Optimization loop

`scripts/optimize.py` runs stages 4-8 with a config, scores concordance, and logs to
`runs/<id>/`. Subagents propose config/model changes; each run is accepted only if the
concordance harness improves on held-out runs.
