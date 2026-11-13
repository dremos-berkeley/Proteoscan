# Proteoscan

Open-source (MIT) DDA proteomics identification pipeline for Thermo Orbitrap data, on macOS and
Linux. It writes PEAKS Studio-compatible exports (PSMs, peptides, protein groups, features,
de novo-only spectra, PTM profile). On synthetic-peptide ground truth (ProteomeTools, 23 pools) it finds
as many true peptides as PEAKS Studio 13.5 with its Deep Learning Boost on HCD scans (825 vs 831 per pool
on ion-trap MS2, 852 vs 841 on Orbitrap MS2) and 3% fewer on ion-trap CID (805 vs 830); on a Jurkat lysate
it reports 91% of the PSMs and peptides and 95% of the protein groups PEAKS reports. Numbers, error rates
and limitations: docs/STATUS.md "Current results".

Pipeline: ThermoRawFileParser -> MS1 feature detection and PEAKS-style precursor list -> Comet ->
protein-restricted second pass -> mokapot rescoring with MS2PIP fragment-intensity (including observed
b/y ion-series coverage), DeepLC retention-time, sequence and protein-support features -> PEAKS-style
FDR filter -> protein inference -> exports and a self-contained HTML results viewer. (De novo agreement
from fine-tuned Casanovo and the extended ion features are used only for training score-imitation
models, `scripts/train_score_model.py`, not by the pipeline.)

See `docs/ARCHITECTURE.md` for the stage design and
`docs/PEAKS_EXPORT_FORMAT.md` for the export format we reproduce.

Small defined mixtures (a few standard proteins, bead digests): the target-decoy FDR needs a decoy
null. Search a restricted FASTA of the expected proteins plus contaminants *and* a background set
(e.g. the human proteome) so that the decoy null keeps on the order of 10^4 spectrum winners, then
read `fdr_estimable` and the sequence counts in `run_stats.json`: when fewer than 10 decoy
sequences pass the threshold (or none, `threshold_above_all_decoys`), the reported 1% is not an
estimate; quote the conservative `sequence_fdr_conservative` = (D+1)/T instead. For PSM-level
claims use `mokapot_q` (mokapot's cross-validated PSM q-value, in `psms_scored.parquet`).
See docs/STATUS.md "Corrections" (P2).

Linux below; macOS (Apple Silicon or Intel): see `docs/INSTALL_MACOS.md`, then check any
install with `python scripts/check_install.py`.

```bash
conda env create -f env/environment.yml
conda activate proteoscan
proteoscan run --config configs/peaks_facility_human.yml --input data/sample.raw --out work/sample
proteoscan run --config configs/peaks_facility_human.yml --input data/sample.raw --out work/sample \
    --ref <PEAKS export dir>      # also writes work/sample/concordance/summary.md
```

`--input` takes Thermo raw or mzML files of one sample. The output directory holds the PEAKS-format
tables and HTML pages in `export/`, the scored PSM table `psms_scored.parquet` and the run's quality
figures in `run_stats.json`.

## Settings that differ from a plain database search

All are fields of the search config (`proteoscan/config.py`); the defaults are what the benchmarks in
docs/STATUS.md used.

| setting | default | meaning |
|---|---|---|
| `precursor_recalibration` | `auto` | A quick first search measures the run's precursor mass error. The real search removes the offset (if at least 1 ppm) and narrows the precursor tolerance to max(5 ppm, 10 robust SD), never wider than `parent_mass_error_tolerance_ppm`. `offset` only removes the offset; `off` searches the measured m/z at the configured tolerance. The values used are in `precursor_recalibration.json` and `run_stats.json`. |
| `feature_linking` | `all` | Every MS1 feature in the isolation window becomes a precursor entry of the MS2 scan, and the rescorer is told which entry is the instrument's own precursor. `primary` keeps only the instrument's precursor: fewer second identifications per scan, lists closer in size to PEAKS's. |
| `comet_precursor_mz_source` | `entry` | Comet searches each entry at its own m/z. `trailer` reproduces searches made before 2026-11-09. |
| `clip_nterm_methionine` | `true` | Proteins are also searched without the initiator Met, as PEAKS does. |

Why the tolerance is narrowed: at 15 ppm a 13C isotope peak of an unmodified peptide matches the same
peptide with a deamidation (the two differ by 19 mDa), and target-decoy rescoring cannot remove such
matches because the sequence is right. Narrowing after calibration removes them (docs/STATUS.md).

Results are reproducible: the same input, config and output directory give identical tables. A
directory written with another config or input is refused; use a new directory or `--overwrite`.

PSM and spectral counts are not comparable with PEAKS. Proteoscan reports every PSM of an accepted
peptide down to -10lgP 20 for every co-isolated precursor entry, about twice as many PSMs per peptide
as PEAKS on the same file. Compare peptides and protein groups, or use `feature_linking: primary`.

## Other commands

```bash
proteoscan convert data/*.raw --out work/mzml                     # Thermo raw -> mzML only
proteoscan concordance --ref <PEAKS export dir> --test work/sample/export --out work/concordance
proteoscan ui -h                                                  # self-contained HTML results viewer
```

`proteoscan search` is the earlier standalone Sage search; the pipeline (`proteoscan run`) uses Comet.
