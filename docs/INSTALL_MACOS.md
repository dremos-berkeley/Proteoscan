# Installing peaks-parity on macOS

Tested state (2026-09-30): the recipe below was assembled on Linux and has **not yet been run
on a Mac**. What was checked, and how, is listed in [What was verified](#what-was-verified).
Run `python scripts/check_install.py` after installing; it executes every tool and a small
end-to-end search, so a Mac-specific problem shows up there first.

## Requirements

| | Apple Silicon (M1-M4) | Intel |
|---|---|---|
| macOS | 14 Sonoma or later (Comet binary); 12-13 only with Comet built from source | 12 or later |
| environment file | `env/environment-macos.yml` | `env/environment-macos-intel.yml` |
| Comet | bioconda `comet-ms=2026011` (native arm64) | build from source (bioconda's osx-64 package contains an arm64 binary) |
| ThermoRawFileParser | 2.0.0-dev self-contained, via `scripts/install_trfp.sh` | same, x86_64 build |
| Rosetta 2 | not needed | n/a |

- RAM: with the current code a 3 GB raw file peaks at about 34 GB (MS2PIP feature stage, see
  [Run time and memory](#run-time-and-memory-for-a-3-gb-raw-file)). 64 GB runs it without
  swapping; 32 GB works with some swapping; 16 GB swaps heavily in that stage until it is
  batched (proposal below). Comet itself needs about 3 GB.
- Disk: about 5 GB for the environment, 1 GB for the MS2PIP model (`~/.ms2pip`), and about
  4-5x the raw-file size per run (a 3.4 GB raw file produced a 3.0 GB mzML and a 7.8 GB
  expanded mzML).
- Internet on the first run: conda packages, the ThermoRawFileParser zip and the MS2PIP
  HCD2021 model (955 MB, downloaded on first use into `~/.ms2pip`).

## 1. Install Miniforge

Open Terminal (not a Rosetta Terminal: `uname -m` must print `arm64` on Apple Silicon).

```bash
curl -fL -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash "Miniforge3-$(uname)-$(uname -m).sh" -b -p "$HOME/miniforge3"
"$HOME/miniforge3/bin/conda" init zsh      # zsh is the macOS default shell; use bash if you changed it
exec zsh                                   # or open a new Terminal window
conda info | grep platform                 # osx-arm64 on Apple Silicon, osx-64 on Intel
```

If you already have Anaconda or Miniconda, that works too, provided `conda info` reports
`osx-arm64` on Apple Silicon and the `conda-forge` channel is used as below. Do not install
into a conda that was itself installed as x86_64 on an Apple Silicon Mac.

## 2. Create the environment

```bash
git clone <peaks-parity repository URL> peaks-parity
cd peaks-parity
conda env create -f env/environment-macos.yml          # Intel: env/environment-macos-intel.yml
conda activate peaksparity
```

This takes 3-10 minutes (about 220 conda packages plus five pure-Python pip packages). All
compiled packages come from conda-forge/bioconda on purpose: the PyPI wheels of torch and
scikit-learn each bundle a private `libomp.dylib`, and the PyPI wheels of lightgbm and
xgboost need Homebrew's libomp. Two OpenMP runtimes in one process abort with
`OMP: Error #15`. Do not `pip install torch`, `lightgbm` or `xgboost` into this env.

### Intel Macs: Comet

bioconda's `comet-ms` for osx-64 (and upstream's `comet.macos.exe` since v2026.01.0) is an
arm64 binary. Build the same Comet version from source:

```bash
xcode-select --install                    # once; installs clang and make
git clone --depth 1 --branch v2026.01.1 https://github.com/UWPR/Comet.git /tmp/Comet
make -C /tmp/Comet
cp /tmp/Comet/comet.exe "$CONDA_PREFIX/bin/comet"
comet -p && rm comet.params.new           # prints "Comet version 2026.01 rev. 1"
```

Fallback without a compiler: Comet v2025.03.0 was the last release with an x86_64 macOS
binary and writes the same `comet.params` keys the pipeline sets.

```bash
curl -fL -o "$CONDA_PREFIX/bin/comet" https://github.com/UWPR/Comet/releases/download/v2025.03.0/comet.macos.exe
chmod +x "$CONDA_PREFIX/bin/comet"
```

The same source build also serves Apple Silicon Macs on macOS 12 or 13.

## 3. ThermoRawFileParser (only needed for .raw input)

```bash
bash scripts/install_trfp.sh
ThermoRawFileParser --version             # 2.0.0.0
```

The script downloads the self-contained .NET 8 build of ThermoRawFileParser 2.0.0-dev from
the CompOmics GitHub release (`osx-arm64` or `osx` x86_64 zip), checks its sha256, unpacks it
into `$CONDA_PREFIX/share/thermorawfileparser-2.0.0-dev`, and puts a wrapper in
`$CONDA_PREFIX/bin/ThermoRawFileParser`. It needs neither mono nor a separate .NET install.
Two macOS details:

- the arm64 executable in the release zip is not code-signed, and Apple Silicon kills
  unsigned binaries (`Killed: 9`). The script ad-hoc signs it with `codesign --force --sign -`
  (`/usr/bin/codesign` ships with macOS);
- it removes the `com.apple.quarantine` attribute, so Gatekeeper does not block it.

This is the same ThermoRawFileParser version bioconda installs on Linux (2.0.0.dev). Why not
conda on the Mac: bioconda's macOS package is 1.4.5, a `noarch` package that depends on
mono, and conda-forge builds mono only for osx-64. Alternatives if the script fails:

| route | commands |
|---|---|
| framework-dependent build + conda .NET | `conda install -c conda-forge dotnet-runtime=8.0` then `bash scripts/install_trfp.sh --route dotnet --force` |
| Homebrew mono + 1.4.5 | `brew install mono`; unzip `ThermoRawFileParser1.4.5.zip` from the GitHub release; wrapper `exec mono /path/ThermoRawFileParser.exe "$@"` |
| Rosetta env (Apple Silicon) | `CONDA_SUBDIR=osx-64 conda create -n trfp -c conda-forge -c bioconda thermorawfileparser=1.4.5`; put that env's `bin` on PATH after `peaksparity` |

Or convert on Linux/Windows and copy the mzML files over: the pipeline accepts mzML directly.

## 4. Install peaks-parity itself

```bash
pip install --no-deps -e .     # editable: the pipeline reads configs/ from the checkout
peaksparity -h                 # convert, search, concordance, run, ui
```

`--no-deps` keeps pip from touching the conda packages (plain `pip install -e .` also works in
this env: every dependency is already satisfied and nothing else is installed). Keep the
install editable: `pipeline.py` loads `configs/features_lumos_120k.yml` relative to the
source tree.

## 5. Check the installation

```bash
python scripts/check_install.py
```

It prints the platform (and whether Python runs under Rosetta), versions of every Python
package and whether more than one OpenMP runtime got loaded, runs `comet` and
`ThermoRawFileParser` (catching wrong-architecture and unsigned binaries), then generates a
tiny synthetic data set (80-protein FASTA, 2,035-spectrum mzML with 450 planted peptides)
and runs the real pipeline on it: MS1 features, precursor list, Comet, `read_comet`, MS2PIP
and DeepLC features, mokapot, PEAKS-style filter, protein inference, CSV export. It passes if
at least 80% of planted peptides and none of the decoy-only "foreign" peptides are reported.
Expected: about 30 s after the first run (the first run downloads the 955 MB MS2PIP model).
Output on the Linux server with the same package versions:

```
== Smoke test
  [ok  ] two-pass stage                     skipped on the synthetic set (use --two-pass to include it)
  [ok  ] stage: MS1 features + precursor list 1.7s
  [ok  ] stage: Comet search                1.0s
  [ok  ] stage: read_comet                  0.1s
  [ok  ] stage: MS2PIP features             4.5s
  [ok  ] stage: DeepLC features             1.1s
  [ok  ] stage: mokapot rescoring           3.8s
  [ok  ] stage: PEAKS CSV export            0.0s
  [ok  ] planted peptides identified        450/450 (100%); PSM rows 467; foreign peptides reported 0 (must be 0)
== PASS: 0 failure(s), 8 warning(s), 17s total
```

The warnings are optional tools that are absent (sage, percolator, dotnet, mono, casanovo,
peptdeep, polars). The pipeline's second, protein-restricted search pass is skipped on the
synthetic set: it re-scores only candidates that are new in pass 2, and on 700 spectra there
are fewer than the 50 confident PSMs DeepLC needs to calibrate them (`--two-pass` shows the
resulting `ReferenceSelectionError`; it is a pipeline limitation on small inputs, identical on
Linux, not an installation problem).

## 6. First run

Edit the FASTA path in the config (the committed one points at the Linux server), set
`threads:` to the number of performance cores, then:

```bash
cp configs/peaks_facility_human.yml my_search.yml
# edit my_search.yml:  database: /Users/<you>/fasta/human_sp.fasta   threads: 8
peaksparity run --config my_search.yml --input /path/run.raw --out work/run1 --processes 8
# lower-memory variant: no MS1 feature detection, one monoisotope-corrected precursor per MS2
peaksparity run --config my_search.yml --input /path/run.raw --out work/run1 --features none
```

`--input` also takes mzML. `--processes` sets the MS2PIP threads; `threads:` in the config sets
Comet's. The run includes the two-pass search (second Comet pass on the proteins accepted in
pass 1; `--single-pass` skips it). Every stage caches its output in `--out`, so an interrupted run resumes. Results: `work/run1/export/db.*.csv` (PEAKS layout) and `work/run1/run_stats.json`.
HTML viewer: `peaksparity ui --export work/run1/export --mzml work/run1/mzml/run.mzML --out report.html`.

## Run time and memory for a 3 GB raw file

Measured on the Linux server (Xeon Gold 6338, shared and heavily loaded at the time) with
the same package versions as the macOS env, for DrewR_Adam_91526_DDA_OT_4: 3.4 GB raw, Lumos,
211 min gradient, 22,785 MS1 / 194,508 MS2 scans, human SwissProt, semi-tryptic, 3 variable
mods. Mac figures are estimates, not measurements (an M2/M3 Pro/Max performance core is
roughly 1.5-2x faster single-threaded than this Xeon; Comet scales with cores).

| stage | server: wall time | server: peak RSS | Mac estimate (8-12 cores) |
|---|---|---|---|
| raw -> mzML (ThermoRawFileParser 2.0.0-dev release zip, 1 thread) | 12.6 min | 4.1 GB | 6-12 min |
| MS1 features (`detect_features`, committed config) | 7.2 min (numba grouping 5.7 min, single-threaded) | 2.7 GB alone | 4-7 min |
| precursor list + expanded mzML (`stage_precursors`, features in-process) | 15.3 min total | **10.7 GB** | 10-15 min |
| same with `--features none` | 4.9 min | **6.7 GB** | 4-6 min |
| Comet, 482,826 queries (32 threads) | 16.7 min | about 3 GB | 25-60 min |
| MS2PIP + DeepLC features, 1.03 M candidates (8 threads, CPU only) | 22.5 min (MS2PIP 15.8, DeepLC 6.6) | **34.3 GB** (MS2PIP) | 15-30 min; same memory, so it swaps on a 16-32 GB Mac |
| mokapot (1.03 M candidates) | 2.3 min | 19.3 GB process (+7.5 GB over the 11.8 GB the process still holds after MS2PIP) | 2-4 min |
| second (protein-restricted) Comet pass + features for new candidates | not measured | | +20-40% of Comet + features |
| filter, proteins, PTM tables, export | not measured separately (seconds on the smoke test) | not measured | |

The memory figures quoted before this measurement ("about 3 GB for the search, up to 27 GB for
MS1 feature detection") need two corrections for the current code:

- Comet itself does stay near 3 GB, but `peaksparity run` is one Python process whose peak is
  set by `stage_precursors` (10.7 GB with MS1 features, 6.7 GB with `--features none`, both
  dominated by holding the whole mzML plus the whole expanded mzML in memory) and by the
  MS2PIP feature stage (see table).
- `detect_features` with the committed `configs/features_lumos_120k.yml` peaks at **2.7 GB**
  on this file, not 27 GB; with `boundary: walk` it peaks at 9.3 GB (the re-extraction over
  +/-1 min windows). The 27 GB in `peaks-parity-work/features/REPORT.md` was measured through
  the evaluation harness, which holds PEAKS's features and extra lookups as well.

Whole run on the server: about 80-90 min of compute for this file (conversion 13, precursors
15, Comet 17, features 23, mokapot 2, plus the second pass). Peak memory of `peaksparity run`
is the MS2PIP stage, 34 GB, not feature detection.

Practical guidance until the MS2PIP stage is batched: 64 GB Macs run the default
configuration; 32 GB Macs run it with some swapping (macOS compresses and swaps rather than
killing the process, so expect the MS2PIP stage to take several times longer); on 16 GB use
`--features none` (saves 4 GB in the precursor stage and about 60% of the candidates, since
feature-derived precursor entries are 324 k of the 483 k queries) and close other
applications, or run the search on a larger machine and only view results on the Mac.
MS2PIP memory scales with the number of candidates, i.e. with queries x `--candidates`.

## Memory: bounding each stage (proposal)

This is a proposal; none of it is implemented. Target: every stage of `peaksparity run`
under about 8 GB on a 3-4 GB raw file, so a 16 GB Mac runs the default configuration.

### Where the memory goes now (S4, measured above)

| component | peak | cause |
|---|---|---|
| `detect_features` (hull config) | 2.7 GB | flat MS1 arrays (32.2 M centroids x 24 B), trace/segment tables, profiles |
| `detect_features` (`boundary: walk`) | 9.3 GB | `_reextract` materialises every (feature, MS1 scan) pair of a +/-1 min window at once: about 200 k features x 250 scans x ~15 float64 columns |
| `stage_precursors` | 10.7 GB | `load_experiment` keeps the whole mzML (3.3 GB), features run in the same process (+2.6 GB), `write_expanded_mzml` builds all 482,826 expanded spectra in an `MSExperiment` before storing (+5.8 GB) |
| MS2PIP features | **34.3 GB** | `predict_batch` returns one result object (theoretical m/z and predicted intensities per ion series) for all 1.03 M PSMs at once; the process keeps 11.8 GB afterwards |
| mokapot | +7.5 GB | full feature matrix and cross-validation copies |

### m/z-slab chunking of MS1 feature detection

The traces detector and the re-quantification are both local in m/z: a trace links centroids
within `trace_ppm` (6 ppm, < 0.015 Th at 2,000 m/z); an isotope chain reaches at most
`(max_isotopes - 1) x 1.00286 / charge_low` = 5.01 Th above its monoisotopic peak; the area
uses isotopes up to `(area_isotopes - 1) x 1.00336 / z` <= 1.01 Th above it. So:

1. **Pass 0, streaming.** Read MS1 spectra one at a time (`pyopenms.MzMLFile().transform`
   with a consumer, or `OnDiscMSExperiment`) and write the flat centroid arrays once to
   disk-backed `numpy.memmap` files: `mz` float64, `int` float32, `scan` int32 (16 B per
   centroid, 0.5 GB for S4), plus the per-scan offsets and RTs. Build a 0.01-Th histogram
   of centroid m/z on the fly.
2. **Choose K slabs** of equal centroid count from the histogram (as `_run_centroided`
   already does with quantiles). Slab k has core `[a_k, b_k)` and is loaded as
   `[a_k - 0.05, b_k + 5.1]` Th: per scan, `searchsorted` on the memmapped m/z (spectra are
   sorted) gives the index range, so only the slab's rows are copied into RAM.
3. **Per slab**: `run_traces` (unchanged numba code) on the slab arrays, then `_reextract`
   on the slab arrays; keep only features whose monoisotopic m/z is in the core `[a_k, b_k)`.
   Every isotope and every XIC such a feature needs lies inside the loaded extension.
4. **Seams.** The one non-local step is the greedy acceptance in `_group`: an elution peak in
   the 5 Th above `b_k` can be claimed as an isotope by a slab-k feature and as a
   monoisotopic peak by a slab-(k+1) feature, and the unsplit run resolves that by global
   score order. Re-run the detector on a seam window `[b_k - 10.2, b_k + 10.2]` around every
   edge and let its features replace both slabs' features with monoisotopic m/z in
   `[b_k - 5.1, b_k + 5.1]`. Inside that zone every competing hypothesis is fully contained
   in the seam window, so the result equals the unsplit run except for interaction chains
   longer than 5 Th (rare; measure it). Seams add about `K x 20 Th / 1,500 Th` of extra work
   (about 20% for K = 16).
5. **Budget.** Resident memory per slab is about `(N / K) x B` with N centroids and
   B ~ 80 bytes per centroid for detection (arrays, `tid`, `order`, profiles, segment
   tables: 2.7 GB / 32 M on S4) and ~150 B per (feature, scan) row for `_reextract` in walk
   mode (7.2 GB for about 48 M rows on S4). Choose `K = ceil(N x B / budget)` with
   `budget = max_memory_gb` (new `detect_features` argument, default 8 GB, 4 GB on 16 GB
   machines). S4 needs K = 1 in hull mode and K = 2
   in walk mode; a 4x larger file needs K = 2-8. Slabs are independent, so they can also run
   in a `spawn` process pool (`workers x slab memory <= budget`), which cuts the 5.7 min
   single-threaded grouping on S4 roughly by the worker count.
6. **Validation**: on S1 and S4, the chunked feature table must equal the unchunked one
   (same `feature_id` set after sorting by m/z and RT; areas identical) for K = 1, 2, 8, 32;
   report any seam differences.

A smaller, exact first step: run `_reextract` over batches of features sorted by m/z, each
batch limited to ~5 M (feature, scan) rows (about 1.3 GB), and concatenate. This alone takes
the walk configuration from 9.3 GB to about 3 GB and needs no seam logic.

### The two larger levers for `peaksparity run`

- **Stream the expanded mzML.** `write_expanded_mzml` should write each expanded spectrum
  through `pyopenms.PlainMSDataWritingConsumer` as it is created instead of collecting an
  `MSExperiment` (-5.8 GB), and `stage_precursors` should read MS2 spectra through
  `OnDiscMSExperiment` and keep only MS1 peaks in memory for the monoisotopic correction
  (about -2 GB). Expected `stage_precursors` peak: about 3-4 GB with features, 2 GB without.
- **Batch the MS2PIP features.** `ms2pip_correlations` can process the candidate table in
  chunks of 100 k PSMs (sorted by scan so the spectrum cache stays effective) and concatenate
  the feature columns; results are identical because predictions are per PSM. Estimated
  peak about 3-4 GB per 100 k-PSM chunk (34 GB / 1.03 M PSMs, measured unbatched), and the
  11.8 GB the process keeps afterwards would go too, taking mokapot's stage to about 8 GB.
  DeepLC can be chunked the same way after calibration. This is the change that matters most
  for laptops.

## Troubleshooting

| symptom | cause and fix |
|---|---|
| `zsh: killed ThermoRawFileParser` or `Killed: 9` | unsigned arm64 binary: `codesign --force --sign - "$CONDA_PREFIX"/share/thermorawfileparser-*/ThermoRawFileParser` (the install script does this) |
| `"ThermoRawFileParser" cannot be opened because the developer cannot be verified` | quarantine flag from a browser download: `xattr -dr com.apple.quarantine "$CONDA_PREFIX"/share/thermorawfileparser-*` |
| `Bad CPU type in executable: comet` | Intel Mac with the bioconda arm64 Comet: build from source (step 2) |
| `dyld: ... built for macOS 14.0 which is newer than running OS` (comet) | Apple Silicon on macOS 12/13: build Comet from source (step 2) |
| `OMP: Error #15: Initializing libomp.dylib, but found libomp.dylib already initialized` | a pip wheel of torch/lightgbm/xgboost/scikit-learn is in the env: `pip list | grep -Ei 'torch|lightgbm|xgboost|scikit'`, `pip uninstall` it, `conda install -c conda-forge <pkg>`; `check_install.py` shows the loaded OpenMP libraries. Last resort: `export KMP_DUPLICATE_LIB_OK=TRUE` |
| `Library not loaded: @rpath/libomp.dylib` (lightgbm/xgboost) | pip wheel of lightgbm/xgboost; reinstall from conda-forge as above (or `brew install libomp` if you insist on pip) |
| `conda env create` fails to solve | use `mamba`/`conda` >= 23.10 (libmamba solver), keep `channels: conda-forge, bioconda` in that order, `conda config --set channel_priority strict`; do not add `defaults` |
| `conda info` says `osx-64` on Apple Silicon | Terminal or conda runs under Rosetta; reinstall Miniforge from a native Terminal (`uname -m` = `arm64`) |
| `UserWarning: OPENMS_DATA_PATH environment variable not found` | Python run without `conda activate peaksparity`; activate the env |
| MS2PIP download fails (firewall, offline Mac) | copy `~/.ms2pip/model_20210416_HCD2021_{B,Y}.xgboost` from a machine that has them |
| `RuntimeError: Failed to calibrate scores ... no target PSMs could be found below 'test_fdr'` (mokapot) | too few confident PSMs (tiny or failed run); check `comet` output and the FASTA path |
| `comet.params keys not found for this Comet version` | a different Comet release is first on PATH: `which -a comet`; use 2026.01.1 (or 2025.03.0 on Intel) |
| process killed / system swaps during a run | lower `--processes` and `threads:`; use `--features none`; see the memory section |
| `ModuleNotFoundError: peaksparity` | `pip install --no-deps -e .` from the repository root, inside the activated env |

## Optional extras

- **Sage**: `conda install -c bioconda sage-proteomics=0.14.7` (osx-arm64 and osx-64 builds
  exist). Not used by `peaksparity run`; note that Sage needs 300-500 GB RAM for the
  semi-tryptic human search, so it is not practical on a Mac anyway.
- **Casanovo** (de novo, not used by `peaksparity run` yet): it requires numpy < 2, which
  conflicts with pyopenms 3.5.0 from bioconda on macOS. Use a separate env:
  `conda create -n casanovo -c conda-forge python=3.11 pytorch=2.10 "numpy<2" && conda activate casanovo && pip install casanovo==5.2.1`.
  (Not solved or tested on macOS.)
- **Percolator**: no macOS conda build; mokapot replaces it.

## Package matrix

Verified 2026-09-30 against the conda-forge/bioconda repodata for osx-arm64 and osx-64 and
the PyPI JSON API. "solve" = part of the `conda create --dry-run` solution for that platform.

| package | version | osx-arm64 (Apple Silicon) | osx-64 (Intel) | notes / fallback |
|---|---|---|---|---|
| python | 3.11.16 | conda-forge, solve | conda-forge, solve | |
| numpy | 2.3.5 | conda-forge, solve | conda-forge, solve | numpy < 2 is impossible with bioconda pyopenms 3.5.0 on macOS (needs numpy >= 2); the codebase passes its 149 tests on numpy 2.3.5 |
| pandas | 2.3.3 | conda-forge, solve | conda-forge, solve | ms2pip/deeplc need < 3 |
| scipy | 1.17.1 | conda-forge, solve | conda-forge, solve | |
| scikit-learn | 1.9.1 | conda-forge, solve | conda-forge, solve | PyPI wheel exists but bundles its own libomp |
| lightgbm | 4.7.0 | conda-forge, solve | conda-forge, solve | PyPI `macosx_12_0_arm64` wheel exists but needs `brew install libomp` |
| xgboost | 2.1.4 | conda-forge (cpu), solve | conda-forge (cpu), solve | required by ms2pip (< 3); PyPI wheel needs Homebrew libomp |
| pyarrow | 21.0.0 | conda-forge, solve | conda-forge, solve | pinned by pyopenms (libarrow 21.0) |
| numba / llvmlite | 0.67.0 / 0.49.0 | conda-forge, solve | conda-forge, solve | PyPI has no macOS x86_64 wheel for these versions |
| pyopenms | 3.5.0 | bioconda py311 `h8a561bf_1`, solve | bioconda py311 `h049568f_1`, solve | PyPI 3.5.0 wheels are `macosx_15_0` only (macOS 15+) |
| pyteomics | 4.7.5 | conda-forge noarch | conda-forge noarch | ms2pip/psm-utils need < 5 |
| pyahocorasick | 2.3.1 | conda-forge, solve | conda-forge, solve | PyPI arm64 and universal2 wheels also exist |
| pyyaml, click, pytest, tqdm, lxml, matplotlib-base | 6.0.3, 8.5.0, 9.1.1, 4.70.1, 6.0.2, 3.10.9 | conda-forge | conda-forge | |
| pytorch | 2.10.0 cpu | conda-forge `cpu_generic`, solve | conda-forge `cpu_mkl` (MKL on llvm-openmp), solve | PyPI: torch 2.6.0 has a `macosx_11_0_arm64` wheel (bundles libomp.dylib); no x86_64 wheel after 2.2.2, and DeepLC 4.5 needs >= 2.6. conda-forge pytorch 2.6 cannot be co-installed with pyopenms 3.5.0 (abseil/icu pins) |
| llvm-openmp | 23.1.2 | conda-forge | conda-forge | the single OpenMP runtime of the env |
| comet-ms | 2026011 (= 2026.01.1) | bioconda, arm64 Mach-O, min macOS 14.0 | bioconda package contains an **arm64** binary: unusable | Intel: build v2026.01.1 from source, or upstream v2025.03.0 `comet.macos.exe` (x86_64, same params keys, same results on the smoke test) |
| ThermoRawFileParser | 2.0.0-dev | GitHub `osx-arm64` zip: self-contained .NET 8, arm64, min macOS 11, apphost unsigned (ad-hoc signed by `install_trfp.sh`) | GitHub `osx` zip: x86_64, signed, min macOS 10.15 | bioconda macOS build is 1.4.5 noarch needing mono; conda-forge mono 6.12.0.199 is osx-64 only; Homebrew mono 6.14.1 has arm64 bottles; conda-forge `dotnet-runtime` 8.0.31 exists for osx-arm64 (for the `net8` zip) |
| ms2pip | 4.2.0 | PyPI pure wheel (+ `ms2rescore-rs` 0.5.0 `macosx_11_0_arm64` wheel); bioconda py311 osx-arm64 build also exists | PyPI pure wheel (+ `ms2rescore-rs` `macosx_10_12_x86_64`); bioconda osx-64 only 4.1.0 | |
| deeplc | 4.5.0 | PyPI pure wheel; bioconda noarch | same | needs torch >= 2.6 (satisfied by conda pytorch; pip does not reinstall it) |
| mokapot | 0.10.0 | PyPI pure wheel; bioconda noarch | same | pip resolves triqler 0.9.2 (Linux reference env has 0.7.3) |
| psm-utils | 1.5.5 | PyPI pure wheel; bioconda noarch | same | |
| other pip-resolved wheels | pydantic-core 2.46.5, rapidfuzz 3.14.6, sqlalchemy 2.1.1, ms2rescore-rs 0.5.0 | arm64 wheels | x86_64 wheels (sqlalchemy: pure fallback) | none use OpenMP |
| sage-proteomics (optional) | 0.14.7 | bioconda, solve | bioconda | |
| casanovo (optional) | 5.2.1 | PyPI pure wheel; needs numpy < 2, pylance 0.15.0 (arm64/x86_64 abi3 wheels) | same | separate env |
| percolator (optional) | - | no macOS conda build | no macOS conda build | mokapot replaces it |

## What was verified

Checked on 2026-09-30 from the Linux server; nothing was run on macOS.

Verified:

- **Conda solves**: `env/environment-macos.yml` solves for osx-arm64 and (as a Rosetta env)
  osx-64, and `env/environment-macos-intel.yml` for osx-64, with
  `CONDA_SUBDIR=<platform> CONDA_OVERRIDE_OSX=14.0 conda create --dry-run` (also with
  `CONDA_OVERRIDE_OSX=12.0` and `11.0` for the conda part). Exact builds are in the package
  matrix above.
- **pip side**: `uv pip compile --python-platform aarch64-apple-darwin` and
  `x86_64-apple-darwin` (`--only-binary :all:`, with the conda-provided packages as
  constraints) resolve ms2pip 4.2.0, deeplc 4.5.0, mokapot 0.10.0, psm-utils 1.5.5 and their
  dependencies to wheels for both architectures; the two lock files are identical except
  numba/llvmlite (no x86_64 wheels; conda provides them).
- **Same pins on Linux**: an env built from the same pins on linux-64 (`conda env create`,
  1.5 min) installed nothing extra with pip (conda's pytorch satisfies DeepLC's `torch>=2.6`),
  passes `scripts/check_install.py` (450/450 planted peptides, 17-30 s) with multiprocessing
  forced to `spawn` (the macOS default), passes the full test suite (149 passed, numpy 2.3.5),
  and `pip install -e .` in it installs only peaksparity.
- **Binaries, by inspecting the Mach-O headers** (architecture, code signature, minimum OS):
  bioconda comet-ms 2026011 osx-arm64 = arm64, signed, macOS >= 14.0; the osx-64 package and
  upstream `comet.macos.exe` (v2026.01.0 onward) = arm64 too; v2025.03.0 `comet.macos.exe` =
  x86_64. ThermoRawFileParser 2.0.0-dev `osx-arm64` zip = self-contained .NET 8.0.28,
  arm64, macOS >= 11.0, apphost unsigned (all dylibs signed); `osx` zip = x86_64, signed,
  macOS >= 10.15. Both link only `/usr/lib` and system frameworks.
- **Comet 2025.03.0** (the Intel fallback, Linux build of the same release) writes the same
  `comet.params` keys and gives identical output on the synthetic set, and the smoke test
  passes with it.
- **ThermoRawFileParser 2.0.0-dev release zip** (Linux build, installed with
  `scripts/install_trfp.sh`) converts the 3.4 GB raw file in 12.6 min / 4.1 GB RSS to an mzML
  within 31 bytes of the one the pipeline has been using.
- **Wheel contents**: PyPI torch 2.6.0 and scikit-learn 1.9.1 macOS arm64 wheels bundle
  `libomp.dylib`; lightgbm 4.7.0 and xgboost 2.1.4 macOS wheels do not bundle it (need
  Homebrew libomp); PyPI pyopenms 3.5.0 wheels are `macosx_15_0` only.

Not verified (needs a Mac):

- that the conda env actually installs and imports on macOS (solver success does not prove
  the dylibs load), and that only one OpenMP runtime ends up loaded (`check_install.py`
  reports it);
- the ad-hoc `codesign` fix for the unsigned ThermoRawFileParser apphost, and that
  ThermoRawFileParser converts a raw file on macOS (the .NET Thermo reader is the same
  managed code as on Linux);
- the Comet source build on Intel and on macOS 12/13;
- run time and memory on Apple Silicon (the table gives server measurements and estimates);
- Casanovo and Sage on macOS;
- `bash scripts/install_trfp.sh` under macOS's bash 3.2 (the script avoids bash-4 features
  and was run on Linux only).
