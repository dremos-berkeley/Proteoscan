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
