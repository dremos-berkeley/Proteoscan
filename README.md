# peaks-parity

Open-source (MIT) DDA proteomics identification pipeline built to reproduce PEAKS Studio's
exports (PSMs, peptides, protein groups, features, de novo-only spectra, PTM profile) from
Thermo raw files, on macOS and Linux.

Status: under construction. See `docs/ARCHITECTURE.md` for the stage design and
`docs/PEAKS_EXPORT_FORMAT.md` for the export format we reproduce.

```bash
conda env create -f env/environment.yml
conda activate peaksparity
python -m peaksparity.cli convert data/*.raw --out work/mzml
python -m peaksparity.cli search work/mzml/*.mzML --config configs/peaks_facility_human.yml --out work/search
python -m peaksparity.cli concordance --ref <peaks export dir> --test <our export dir> --out work/concordance
```
