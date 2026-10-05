# Proteoscan

Open-source (MIT) DDA proteomics identification pipeline for Thermo Orbitrap data, on macOS and
Linux. It writes PEAKS Studio-compatible exports (PSMs, peptides, protein groups, features,
de novo-only spectra, PTM profile) and, on synthetic-peptide ground truth, matches PEAKS Studio 13.5
with its Deep Learning Boost in sensitivity and accuracy (docs/STATUS.md).

Pipeline: ThermoRawFileParser -> MS1 feature detection and PEAKS-style precursor list -> Comet ->
protein-restricted second pass -> mokapot rescoring with MS2PIP fragment-intensity (including observed
b/y ion-series coverage), DeepLC retention-time, sequence and protein-support features -> PEAKS-style
FDR filter -> protein inference -> exports and a self-contained HTML results viewer. (De novo agreement
from fine-tuned Casanovo and the extended ion features are used only for training score-imitation
models, `scripts/train_score_model.py`, not by the pipeline.)

Status: under construction. See `docs/ARCHITECTURE.md` for the stage design and
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
python -m proteoscan.cli convert data/*.raw --out work/mzml
python -m proteoscan.cli search work/mzml/*.mzML --config configs/peaks_facility_human.yml --out work/search
python -m proteoscan.cli concordance --ref <peaks export dir> --test <our export dir> --out work/concordance
```
