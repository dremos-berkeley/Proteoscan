# peaks-parity

Open-source (MIT) DDA proteomics identification pipeline that imitates PEAKS Studio's
exports. See docs/ARCHITECTURE.md. Python package `peaksparity/`, CLI entry points via
`python -m peaksparity.<module>`.

Rules
- Never inspect, decompile or copy the PEAKS binary. Only its published papers and its
  exported results are used.
- Reference exports live under /drive_a/drew/vendor-exports (read-only) and
  /drive_a/drew/public-data (PXD downloads). Heavy intermediates go to
  /drive_a/drew/peaks-parity-work/<stage>/, never into the repo.
- Conda env: `conda activate peaksparity` (see env/environment.yml). Until it exists,
  `conda activate openms` has pyopenms/pyteomics/sage/comet/percolator.
- Concordance is measured with `python -m peaksparity.concordance`; PSM key is
  (source file, raw scan, exact modified sequence, charge). Report recall and precision,
  never "concordance" alone.
- Shared server (128 threads): use <= 32 threads per job, and never more than 2 GPUs.
- Tests: `pytest tests/` (reference-dependent tests skip when data is absent).
