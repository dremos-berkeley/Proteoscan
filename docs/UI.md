# Results viewer (`proteoscan.ui`)

`proteoscan.ui` turns a PEAKS-format DB-search export (a PEAKS Studio export or one written
by our pipeline) into **one self-contained HTML file**. The CSS, JavaScript, tables and peak
lists are all inline. The page makes no network requests and uses no CDN or external
library, so it opens straight from disk (`file://`) in Safari, Chrome or Firefox. It shows
the same result views as a DB-search result in PEAKS Studio, but with its own layout and
styling.

```bash
conda activate proteoscan
python -m proteoscan.ui --export <export dir> \
    [--mzml run.mzML ...] [--fasta db.fasta] [--ref <PEAKS export dir>] \
    --out report.html
# same thing through the main CLI
python -m proteoscan.cli ui --export <export dir> ... --out report.html
```

| option | meaning |
|---|---|
| `--export` | Export directory with `db.psms.csv`, `db.peptides.csv`, `db.proteins.csv` and `db.protein-peptides.csv`. `db.ptmProfile.csv` and `db.summary.html` are used when present. Any missing table is skipped. |
| `--mzml` | Indexed mzML file(s) for the spectrum view. Each file is matched to the `Source File` column by base name without extension (`X.raw` pairs with `X.mzML`). A single mzML is used for a single source file whatever its name. |
| `--fasta` | Protein FASTA for the coverage maps. Defaults to `<export>/db.proteins.fasta` when that file exists (PEAKS writes one). |
| `--ref` | Reference PEAKS export. Adds the Concordance panel and an `In ref` column to the PSM, peptide and protein tables. |
| `--top-peaks` | Peaks kept per embedded spectrum, most intense first (default 150). |
| `--max-spectra-mb` | Cap on embedded peak data (default 40 MB). |
| `--frag-tol`, `--frag-unit` | Fragment tolerance for annotation. The default comes from the export's search parameters, or is 0.5 Da if the export has none. It can also be changed in the page. |
| `--title` | Page title (default: the export directory name). |

Example reports (Jurkat sample 4):

```bash
MZ=/drive_a/drew/openms-data/260921-DrewR_Adam_91526/mzml/DrewR_Adam_91526_DDA_OT_4.mzML
REF="/drive_a/drew/vendor-exports/DrewR_Adam_91526_DDA_OT/Analysis 1 - Sample 4"
FA=/drive_a/drew/proteoscan-work/fasta/human_sp_peaks_merged.fasta
python -m proteoscan.ui --export "$REF" --mzml $MZ --fasta $FA \
    --out /drive_a/drew/proteoscan-work/ui/report_peaks_ref.html          # 12.3 MB, 3398 spectra, ~25 s
python -m proteoscan.ui --export /drive_a/drew/proteoscan-work/exp_peaksfeat/export_tda \
    --mzml $MZ --fasta $FA --ref "$REF" \
    --out /drive_a/drew/proteoscan-work/ui/report_pipeline.html           # 19.7 MB, 4882 spectra, ~30 s
```

Most of the build time goes into opening the 3 GB mzML index (about 20 s). Peak memory is
about 0.7 GB.

## What the page shows

**Summary**
* Identification counts computed from the CSVs:
  * PSMs, as unique (scan, feature) pairs, the way the PEAKS summary counts them.
  * Peptides, I/L-collapsed.
  * Sequences.
  * Protein groups, all proteins and top proteins.
  * Identified scans.
  * Embedded spectra.
  When `db.summary.html` is present, its own values are shown next to the computed ones. For
  the sample 4 reference, every count matches the summary.
* Filter thresholds and FDRs, all search parameters, and the fragment tolerance the viewer
  is using.
* Plots:
  * PSM and peptide -10lgP histograms, with the peptide threshold marked.
  * Precursor error (ppm) histogram.
  * ppm vs m/z scatter.
* Missed-cleavage table, computed with the export's enzyme rule and set beside the summary's
  Table 5 when that exists.
* PTM profile: `db.ptmProfile.csv`, or computed from the peptide tables when that file is
  missing.
* **Concordance** (with `--ref`), built from `proteoscan.concordance.compare`:
  * Recall and precision at PSM, peptide and protein-group level, plus the diagnostic levels.
  * The breakdown of why PSMs differ.
  * Browsable lists of reference-only and test-only PSMs, peptides, protein groups and top
    accessions.
  * Clicking a listed PSM opens its scan annotated with that list's peptide. This works for
    reference-only PSMs too, provided the scan's peaks are embedded.

**Proteins.** A table of all protein rows showing group, top flag, accession, gene, -10lgP,
coverage, #peptides, #unique, #spec, area and description. A "top proteins only" switch hides
the rest; non-top proteins are dimmed. Clicking a protein opens three things:
* The sequence coverage map. Covered residues are shaded, modified residues are underlined in
  orange, and supporting peptides are drawn as bars under each line (green if unique, grey if
  shared) with dots at modified positions. Clicking a bar opens that peptide's spectrum.
* The other proteins in its group.
* The supporting peptides from `db.protein-peptides.csv`.

**Peptides, PSMs.** Every export column is shown; the column chooser can hide some. Features:
* Sort on any column header.
* Text filter on sequence, accession, gene or PTM. Several words are ANDed. Peptide filters
  also match the unmodified sequence.
* Numeric filter on -10lgP.
* Pagination (25 to 1000 rows per page).
* Export of the filtered rows to CSV.

Clicking a peptide shows the spectrum of its best PSM, its PSM list and its proteins.
Clicking a PSM shows its spectrum. The PSM table can be limited to PSMs whose peaks are
embedded.

**Spectrum view.** Everything here is drawn as inline SVG by `report.js`:
* The annotated MS2 spectrum, with b ions in blue, y ions in red and neutral losses in
  purple. Ion charges are 1+, plus 2+ (on by default when the precursor charge is 3 or more).
  Neutral losses (H2O on STED, NH3 on RKNQ) can be switched on.
* Drag to zoom, double-click to reset. Hovering over a peak shows its annotation and mass
  error.
* A fragment-error plot.
* A sequence ladder marking the matched b and y ions.
* A fragment table with every theoretical ion; matched cells are coloured and show the
  observed m/z and error on hover.
* Tolerance and unit can be edited.

Theoretical masses use exact modification deltas: the rounded PEAKS notation is mapped
through `seqnorm.MOD_MASSES`, so `M(+15.99)` becomes +15.994915. If a scan's peaks are not
embedded, the view still lists the theoretical fragments.

URL hash deep links:
* `#tab=proteins&protein=<accession>`
* `#tab=peptides&peptide=<sequence>`
* `#tab=psms&psm=<row index>`, where `psm=best` picks the best PSM with an embedded spectrum
* `#spectrum=<psm row>`, which opens the spectrum in a dialog
* `#bench`, which prints table and spectrum timings for the current browser

## Spectra

Peaks are read with pyopenms `OnDiscMSExperiment`, with random access by raw scan number. The
reader first tries spectrum index `scan - 1` and checks that the native id ends in
`scan=<N>`. If that check fails, it indexes every native id once.

Only the best PSM of each peptide is embedded: the peptide table's representative
`(Source File, Scan)`, taken in descending peptide score. Other PSMs on the same scan
(chimeric feature duplicates, I/L variants) reuse it. Each spectrum keeps its `--top-peaks`
most intense peaks, with m/z to 4 decimals and intensity as % of the base peak. Embedding
stops once the peak data reach `--max-spectra-mb`, which drops the lowest-scoring peptides
first. The Summary notes say how many spectra were skipped.

## Scale

The page renders one page of rows at a time, and table filters are debounced, so large
tables stay responsive. The spectrum JSON is parsed only when a spectrum is first needed.

Measured with `#bench` in headless Firefox on the 15.5k-PSM pipeline report:

| operation | time |
|---|---|
| parse the 10 MB spectrum block | 132 ms |
| sort 15.5k PSMs | 17–30 ms |
| text filter | 5–28 ms |
| render a 1000-row page | 70 ms |
| draw an annotated spectrum | 5 ms |

## Files

* `proteoscan/ui.py`: data collection (tables, counts, PTM and missed-cleavage summaries,
  sequences, spectra, concordance) and HTML assembly. The JSON sits in
  `<script type="application/json">` blocks with `<` escaped.
* `proteoscan/ui_assets/report.js`, `report.css`: the page, inlined at build time. The
  script is ES2017, with no modules and no dependencies.
* `tests/test_ui.py`: builds a synthetic export, a tiny mzML and a FASTA, then checks:
  * the embedded JSON, tables, spectra, top-N, the size cap and the concordance flags;
  * the CLI subcommand;
  * that the script parses in node;
  * that the fragment calculator in `report.js` matches b/y masses computed in Python.

# PEAKS-style HTML export pages (`proteoscan.summary`)

`proteoscan.summary` writes the static HTML pages that PEAKS Studio puts next to its CSVs, with
the same structure, so a PEAKS user reads them the same way:

* `db.summary.html` has the PEAKS sections "1. Notes", "2. Result Statistics",
  "3. Experiment Control" and "4. Other Information". It carries Figures 1 to 6 and Tables 1 to 6
  with PEAKS's own numbering and captions.
* `db.proteins.html` opens with the protein list, ordered by group with top proteins first and
  non-top proteins dimmed. Each protein then gets one block, as in PEAKS, holding:
  * the description;
  * the group, -10lgP, coverage, #peptides and #unique;
  * a coverage map image;
  * the supporting-peptide table;
  * the annotated spectrum of the best unique PSM.
* `db.peptides.html` holds the peptide list with the PEAKS columns. Each peptide that has an
  image links to a block below the table with the annotated MS2 spectrum of its best PSM.
* `export.css` and `img/*.png` hold our own stylesheet and figures (matplotlib). Nothing comes
  from PEAKS, and the pages load nothing from the network.

The pages and images follow PEAKS's layout, but every file is produced by our own code.
`db.summary.html` can be read back with `peaks_export.parse_summary_html`, the same parser
that reads a PEAKS summary (Tables 1, 2, 3, 5 and 6).

```bash
python -m proteoscan.summary --export <export dir> --out <dir> \
    [--mzml run.mzML ...] [--fasta db.fasta] [--psms-scored psms_scored.parquet] \
    [--features ms1_features.parquet ...] [--precursors precursors.parquet ...] \
    [--config search.yml] [--mgf run.mgf] [--max-spectra 2000] [--processes 8] [--no-images]
python -m proteoscan.cli summary ...        # same thing
```

`pipeline.run` calls the module at the end of every run (inside a try/except) and writes the
pages into the export directory. It also saves `search_config.yml` in the run directory. Use
`--no-summary` to skip both.

## Inputs

All inputs except `--export` are optional. When an input is missing, the figure or number that
needs it is left out, and the page says why in "1. Notes" or in place of the figure.

* **Run directory.** If the export sits in a pipeline run directory (`<run>/export`), these
  files are picked up automatically:
  * `<run>/psms_scored.parquet`
  * `<run>/<stem>/precursors.parquet`
  * `<run>/<stem>/ms1_features.parquet`
  * `<run>/search_config.yml`
* **PEAKS export.** A `*.mgf` inside the export directory (PEAKS writes one) supplies the
  `#Chimera` count when there is no precursor list.
* **Search parameters** (Tables 2 and 6) come from `--config` or the run's `search_config.yml`.
  Failing that, they come from the export's own `db.summary.html`, whose lines are copied
  verbatim.

| item | source | definition |
|---|---|---|
| Fig. 1 FDR curve | `psms_scored.parquet` | Two panels, FDR = D/T against the number of targets kept. Left: unique sequences (I = L, best PSM each), the peptide-FDR unit. Right: PSMs, one per (file, scan, peptide, z). The dashed line marks the current threshold. Omitted without scored decoys. |
| Fig. 2 (a) score histogram, (b) score vs ppm | `psms_scored.parquet`, otherwise `db.psms.csv` (targets only) | Targets and decoys stacked in 1-unit bins. The scatter is downsampled to 60k points per class. |
| Table 1 #MS1, #MS/MS | mzML (pyopenms metadata) | Number of MS1 and MS2 spectra. |
| Table 1 #Chimera | `precursors.parquet`, or the PEAKS MGF | Scans with more than one precursor entry (more than one MGF entry with the same `RAWSCANS`). |
| Table 1 #Features | `ms1_features.parquet` | Rows (every detected feature). |
| Table 1 Identified, #Peptides, #Sequences, #Proteins | export CSVs | Definitions in docs/PEAKS_EXPORT_FORMAT.md. #PSMs = unique (file, scan, Feature Id). #Features = distinct non-zero Feature Id. #Peptides = distinct modified sequences with I = L. #Sequences = distinct stripped sequences with I = L. Groups / All / Top from `db.proteins.csv`. |
| Fig. 3 Venn | per-sample columns | "Not applicable to only one sample" (PEAKS wording) for single-sample exports. Drawn for 2 or 3 samples. |
| Fig. 4, 5 feature m/z, RT, log10 area | `ms1_features.parquet` ("all"), `<Sample>.db.features.csv` ("de novo+db", "db") | Overlaid histograms of the nested sets. |
| Table 2 | `psms_scored` threshold, otherwise the export summary, otherwise the lowest peptide score | Peptide -10lgP is the lowest score with sequence FDR ≤ `peptide_fdr` (same rule as `pipeline.peaks_filter`). Protein -10lgP, unique peptides and de novo ALC come from the config. |
| Table 3 FDRs | `psms_scored.parquet` | PSM: decoy over target PSMs of the accepted sequences, down to -10lgP 20 (`peaks_filter` definition). Sequences: D/T at the threshold. Protein group: target and decoy groups are built separately from the sequences above the threshold, with the PEAKS rules (Σ s_k/k score, greedy grouping, ≥ 1 group-unique sequence, protein -10lgP ≥ 15). Without decoys: `N/A`, plus the value the export's own summary reports. |
| Table 3 De Novo Only Spectra | `<Sample>.db.denovoOnly.csv` | Number of rows. |
| Table 4 | `db.ptmProfile.csv` | Sites printed as "confident (all)", as PEAKS does. |
| Fig. 6 (a), (b) | `db.psms.csv` | `ppm` column. When it differs from the error computed from `m/z` and `Mass` by more than 1 ppm (PEAKS recalibrates), both are drawn as "Before / After s/w calibration". |
| Table 5 | `db.peptides.csv` rows | Missed cleavages under the export's enzyme rule (trypsin, not before P). |

Spectrum images use b/y ions (2+ as well when z ≥ 3) matched to the most intense peak within
the fragment tolerance, together with a residue ladder and a fragment-error panel. They are
drawn for the `--max-spectra` (default 2000) highest-scoring peptides, plus the best unique
peptide of every protein for the protein blocks. Both pages say how many peptides were left
without an image. Coverage maps show 100 residues per line:
* covered residues are shaded and bold;
* peptides are drawn as bars (blue unique, grey shared);
* variable modifications are coloured squares, with a legend;
* in proteins longer than 1500 residues, three or more consecutive uncovered lines collapse
  into a single "no identified peptides" line.

Images are drawn in up to `--processes` (≤ 8) spawned worker processes and palette-quantised
PNGs. Old `spec_*`/`cov_*` images in `img/` are deleted before redrawing.

## Jurkat sample 4 results

Inputs:
* Pipeline: `/drive_a/drew/peaks-parity-work/pipeline_jurkat_features/export`, written to
  `/drive_a/drew/peaks-parity-work/summary_html/jurkat/`.
* PEAKS reference: written to `.../summary_html/peaks_ref/`.

| | pipeline | PEAKS reference export |
|---|---|---|
| time (8 processes) | 82 s | 85 s |
| images | 4578 (2955 spectra, 1615 coverage maps) | 3630 (2385 spectra, 1238 coverage maps) |
| img/ | 53 MB | 42 MB |
| db.proteins.html / db.peptides.html | 6.5 / 3.8 MB | 4.4 / 2.5 MB |

On the reference export, every number we compute matches the PEAKS `db.summary.html`:
* Table 1: MS1 22785, MS/MS 194508, #Chimera 53263 (from the MGF), PSMs 7424, scans 7322,
  features 3386, peptides 3418, sequences 3230, groups/all/top 990/1238/1065.
* Table 3: De Novo Only Spectra 101551.
* Table 5: 3280/183/9/0/0.
* Table 2: Peptide -10lgP 23.02.

Two items cannot be computed from the reference:
* `#Features` (81484): PEAKS's full feature list is not exported.
* The Table 3 FDRs: the reference has no scored decoys. They are shown as N/A, together with
  the values PEAKS reports.

Known differences from PEAKS:
* No "_compare" image of the PSM against the predicted spectrum.
* No help-icon links.
* Coverage maps and spectra are drawn in our own style.
* The protein blocks follow the table order, not PEAKS's apparently random order.
* The FDR curve has a second (PSM) panel.
