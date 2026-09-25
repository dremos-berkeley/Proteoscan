# PEAKS Studio DB-search export format

Reverse-documented from the four reference exports in
`/drive_a/drew/vendor-exports/DrewR_Adam_91526_DDA_OT/Analysis 1 - Sample {1..4}/`
(PEAKS Studio, DB search with Deep Learning Boost, human Swiss-Prot 20,597 entries,
semi-specific trypsin, CAM fixed, Deamidation (NQ) and Oxidation (M) variable).
Everything below was checked against those files. `peaksparity.peaks_export` reads and writes
this format, and a test confirms that read then write reproduces all 32 CSVs byte for byte.

## General encoding

* UTF-8 (all ASCII in the references), no BOM, `\n` line endings, trailing newline after the
  last row, no embedded newlines or quote characters in any value.
* PEAKS uses its own CSV writer. **Whether a value is quoted depends on the column, not on the
  value.** The header row quotes every column name, with one exception: `tag(>=0%)` in
  `denovoOnly.csv` is written bare.
  * `always`: every value is quoted, including empty ones (`""`). This applies to `Source File`
    in psms/features/denovo, `Accession` in psms/proteins/protein-peptides, `PTM`
    everywhere, and `Gene`, `Description`, `Name`, `Position`.
  * `never`: values are bare, even text with spaces.
  * `nonempty`: quoted when there is a value and bare when empty. This applies only to
    `AScore` in psms, peptides and protein-peptides. In modifiedPeptide `AScore` is always
    non-empty, so we can't tell whether it is `always` or `nonempty`.
  * Values that contain commas (`Description`, e.g. `"Insulin, isoform 2 ..."`) sit in
    `always` columns. The writer also quotes a `never` value if it contains `,` or `"`.
* The same logical column is quoted differently in different files. `Source File` is quoted in
  psms but bare in peptides and protein-peptides. `Accession` is quoted in psms and proteins but
  bare in peptides, features and denovoOnly.
* Missing values:
  * A missing number is an empty field. The one exception is `MS2 Correlation`, which writes
    `NaN` instead (seen in samples 2 to 4).
  * A missing text value is empty: `""` in an `always` column, nothing in a `never` column.
  * `-` is never used as a missing marker.
  * `0` means "none" for some columns: `Feature Id`, `-10LgP` and `ALC (%)` in features,
    and `raw m/z` = `0.00000`.
* Number formats (the `kind` in `peaks_export`):

| kind | Java equivalent | examples |
|---|---|---|
| `int` | `Integer.toString` | `12571` |
| `fixedN` | `String.format("%.Nf")`: HALF_UP on the shortest decimal | `36.28`, `1281.6313`, `0.00000`, `-0.6` |
| `sci` | `DecimalFormat("0.##E0")`: at most 3 significant digits, trailing zeros dropped, no `+` | `6.58E9`, `1.1E5`, `3E7` |
| `jfloat` | `Float.toString` from JDK 18 or earlier: plain for 1e-3 <= abs(x) < 1e7, otherwise `d.dddE<n>`; integer-valued floats print every integer digit | `129882.984`, `4531940.0`, `3.0000496E7`, `2.95233664E8` |
| `dec2` | `DecimalFormat("#.##")`: no leading zero | `.98`, `57.02` |
| `bool` | | `true`, `false` |

The `jfloat` columns show that PEAKS runs on an older JDK. `2.95233664E8` is not the shortest
round-trip form (`2.9523366E8`), which JDK 19+ would print. Our emulation reproduces all 129,858
jfloat cells in the references.

* **Protein accessions** are `ACCESSION|ENTRY_NAME` (the `sp|`/`tr|` prefix is dropped).
  The separator for several accessions depends on the file:
  * `;` in psms and denovoOnly
  * `:` in peptides, features and modifiedPeptide
  * protein-peptides and proteins have one accession per row.
* **PTM name lists** also use different separators:
  * psms: `;` with no space, one entry per modified site, repeats kept, e.g.
    `Oxidation (M);Deamidation (NQ)`, `Deamidation (NQ);Deamidation (NQ)`.
  * proteins and denovoOnly: `; ` with a space, e.g. `Deamidation (NQ); Oxidation (M)`.
  * `Carbamidomethylation` is listed as a PTM even though it is a fixed modification.
* **Peptide strings.** Mass deltas follow the residue, with a sign and 2 decimals:
  `M(+15.99)`, `N(+0.98)`, `Q(+0.98)`, `C(+57.02)`. Every C carries `(+57.02)` explicitly.
  The references contain no N-terminal modifications. The expected form, `(+42.01)PEPTIDE`,
  is supported by `seqnorm` but not observed. `Length` counts residues only.
* **I/L.** When the database has both I and L variants of a peptide with identical scores,
  PEAKS reports both as separate rows: the same scan and Feature Id, e.g.
  `LIFAGK`/`LLFAGK` or `VEIIAN(+0.98)DQGNR`/`VEILAN(+0.98)DQGNR`. The summary's
  `#Peptides` and `#Sequences` count these I/L-collapsed, so they are smaller than the number
  of rows in peptides.csv (sample 4: 3418 vs 3472 rows).
* **Scan numbers.** `Scan` is the Thermo raw scan number (MGF `RAWSCANS=`). MGF `SCANS=` is a
  PEAKS-internal query index, not the scan.

## db.psms.csv

This file has one row per (scan, feature, peptide). It holds **every PSM of every peptide that
passed the peptide-level filter**, including PSMs that score below the peptide -10lgP threshold.
In sample 1, 149 of 834 rows are below 28.36; the lowest score is about 20.

The same scan and peptide appear once per co-isolated feature (chimeric queries), so rows
repeat with a different `Feature Id`/`m/z`/`Area`. The summary's `#PSMs` equals the number of
unique (Scan, Feature Id) pairs: 794 for sample 1 against 834 rows. The concordance key
(file, scan, sequence, z) collapses these repeats: 825 keys for sample 1.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Peptide | str | str / never | modified sequence, PEAKS mod notation |
| -10LgP | float | fixed2 / never | PSM score |
| Mass | float | fixed4 / never | theoretical monoisotopic peptide mass (neutral) |
| Tag Length | int | int / never | de novo tag length supporting the PSM |
| CAA (%) | float | fixed1 / never | % confident amino acids (de novo agreement), 0 to 100 |
| Length | int | int / never | residues |
| Delta RT | float | fixed4 / never | observed minus predicted RT (min, DL boost); range about -137 to +49 |
| MS2 Correlation | float | fixed2 / never, NA=`NaN` | correlation with the predicted spectrum; `NaN` when not computed |
| ppm | float | fixed1 / never | precursor error, apparently **after PEAKS recalibration**. Computing it from `m/z` gives a near-constant ~+8 ppm offset in sample 4 (r = 0.9997) |
| m/z | float | fixed5 / never | observed precursor m/z (feature-based, not recalibrated) |
| raw m/z | float | fixed5 / never | instrument precursor m/z from the scan header. `0.00000` in about 54% of rows. When non-zero it usually equals `m/z` but is sometimes one isotope off (differs by about ±1.003/z) |
| z | int | int / never | charge |
| RT | float | fixed4 / never | feature apex RT (min). Equals `raw RT` when `Feature Id`=0 |
| raw RT | float | fixed4 / never | RT of the MS2 scan (min) |
| Area | float | sci / never | feature area. Empty when `Feature Id`=0 |
| Scan | int | int / never | raw scan number |
| Source File | str | str / **always** | e.g. `DrewR_Adam_91526_DDA_OT_1.raw` |
| Feature Id | int | int / never | id in features.csv. `0` = no feature |
| Accession | str | str / **always** | `;`-joined |
| PTM | str | str / **always** | `;`-joined PTM names per site, `""` if none |
| AScore | str | str / **nonempty** | `<AA><pos>:<PTM name>:<score>` joined with `;`, e.g. `M3:Oxidation (M):1001.00;N15:Deamidation (NQ):0.00`. `1001.00` is the unambiguous maximum here (modifiedPeptide uses `1000.00`) |
| Positional Confidence | str | str / never | space-separated per-residue confidence, one token per residue |
| tag(>=0.0%) | str | str / never | the de novo tag at the 0% cutoff; equals `Peptide` in all references. The header embeds the tag-confidence setting |

## db.peptides.csv

This file has one row per accepted modified peptide. I/L variants get separate rows. The values
come from the peptide's best PSM: `-10LgP` matches the best PSM score in 100% of rows and `Scan`
in 97.6%, with the remainder being ties.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Peptide, -10LgP, Mass, Tag Length, CAA (%), Length, Delta RT, MS2 Correlation, ppm, m/z, z, RT | | as psms | from the representative PSM. There is no `raw m/z`/`raw RT` |
| Area Sample 1 | float | sci / never | peptide area in that sample. Empty if none. The column name contains the sample name |
| Scan | int | int / never | scan of the representative PSM |
| Source File | str | str / **never** | bare here, unlike psms |
| #Feature | int | int / never | number of features |
| #Feature Sample 1 | int | int / never | per-sample count; equals `#Feature` for single-sample runs |
| Accession | str | str / never | **`:`-joined** |
| PTM | str | str / always | as psms |
| AScore | str | str / nonempty | as psms |
| Positional Confidence, tag(>=0.0%) | str | str / never | as psms |
| Found By | str | str / never | `DB Search` in all references |

## db.proteins.csv

This file has one row per protein. Rows are **not sorted by group**, and `Protein Group` ids are
not contiguous in file order. A group can have several `Top=true` proteins (indistinguishable
top proteins). In sample 4, 51 of 990 groups have 2 to 10 tops. Summary Table 1: Groups = the
number of distinct `Protein Group`, All = rows, Top = rows with `Top=true`.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Protein Group | int | int / never | group id |
| Top | bool | bool / never | top (representative) protein of its group |
| Accession | str | str / always | single accession |
| Gene | str | str / always | gene name from the FASTA `GN=` |
| -10LgP | float | **jfloat** / never | protein score, e.g. `133.34937`, `28.5408` |
| Coverage(%) | float | fixed2 / never | sequence coverage |
| Coverage(%) Sample 1 | float | fixed2 / never | per-sample coverage |
| #Peptides | int | int / never | |
| Area Sample 1 | float | sci / never | empty if none |
| #Unique | int | int / never | always <= #Peptides |
| #Spec Sample 1 | float | **jfloat** / never | a spectrum count, but a float: `491.0` |
| PTM | str | str / always | **`; `-joined** unique PTM names |
| Average Mass | float | **jfloat** / never | e.g. `14728.311`, `529021.0`, `15388.0205` |
| Description | str | str / always | FASTA description (contains commas) |

## db.protein-peptides.csv

This file has one row per (protein, peptide).

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Protein Group | int | int / never | |
| Accession | str | str / always | |
| Gene | str | str / always | |
| Top | bool | bool / never | |
| Peptide | str | str / never | |
| Extended Peptide | str | str / never | `R.GFFYTPK.T`. Flanking residues from the protein. `-` is expected at a protein terminus but was not seen |
| Unique | bool | bool / never | peptide unique to this group |
| -10LgP | float | fixed2 / never | |
| Mass | float | fixed4 / never | |
| Length | int | int / never | |
| ppm | float | fixed1 / never | |
| m/z | float | fixed5 / never | |
| z | int | int / never | |
| RT | float | fixed4 / never | |
| Scan | int | int / never | |
| Source File | str | str / never | |
| Area Sample 1 | float | sci / never | |
| #Feature | int | int / never | |
| Start, End | int | int / never | 1-based inclusive positions in the protein |
| PTM | str | str / always | |
| AScore | str | str / nonempty | |
| Found By | str | str / never | `DB Search` |

## <Sample>.db.features.csv

Only features that carry a DB or de novo assignment are exported. Every row has exactly one of
`DB Peptide` or `Denovo Peptide`. Sample 4 has 23,540 rows, while the summary says 81,484
features were detected. `Feature Id`s are unique, and every non-zero psms `Feature Id` appears
here.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Source File | str | str / always | |
| Feature Id | int | int / never | |
| m/z | float | fixed5 / never | monoisotopic m/z |
| z | int | int / never | 1 to 9 observed |
| RT | float | fixed4 / never | apex (min) |
| RT Begin, RT End | float | fixed4 / never | elution window |
| Area | float | **jfloat** / never | e.g. `129882.984`, `3.0000496E7` |
| Quality | float | fixed2 / never | always `0.00` in the references |
| DB Peptide | str | str / never | **`;`-joined alternatives** (I/L variants, e.g. `LHLVIR;LHLVLR;IHIVLR`); empty if de novo |
| -10LgP | float | fixed2 / never | `0.00` when no DB peptide |
| Denovo Peptide | str | str / never | empty if DB |
| ALC (%) | float | fixed1 / never | `0.0` when DB |
| ppm | float | fixed1 / never | |
| Accession | str | str / never | **`:`-joined**; empty for de novo |

## <Sample>.db.denovoOnly.csv

This file has one row per de novo-only spectrum with ALC >= the de novo threshold (50%). Its row
count equals the summary's "De Novo Only Spectra". About 1% of its scans also appear in psms,
presumably through other features of chimeric scans.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Source File | str | str / always | |
| Scan | int | int / never | raw scan |
| Peptide | str | str / never | de novo sequence, PEAKS mod notation |
| Tag length | int | int / never | lower-case `l`. Equals Length in the references |
| ALC (%) | float | fixed1 / never | average local confidence |
| Length | int | int / never | |
| m/z | float | fixed5 / never | |
| z | int | int / never | |
| RT | float | fixed4 / never | |
| Area | float | **jfloat** / never | empty when `Feature Id`=0 |
| Mass | float | fixed4 / never | |
| ppm | float | fixed1 / never | |
| PTM | str | str / always | **`; `-joined**, repeats kept |
| local confidence (%) | str | str / never | space-separated per residue |
| mode | str | str / never | fragmentation, `HCD` |
| tag(>=0%) | str | str / never | **header not quoted**; the cutoff is written `0%` (psms/peptides use `0.0%`) |
| Accession | str | str / never | `;`-joined, usually empty (de novo sequence found in the database) |
| Feature Id | int | int / never | `0` = none |

## db.ptmProfile.csv

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Name | str | str / always | PEAKS PTM name, e.g. `Deamidation (NQ)`, `Carbamidomethylation` |
| Mass (Da) | float | **dec2** / never | `.98` (no leading zero), `15.99`, `57.02` |
| Position | str | str / always | residues, e.g. `NQ` |
| # peptides, # proteins | int | int / never | |
| # confident modification sites | int | int / never | the number in parentheses in summary Table 4 |
| # modification sites | int | int / never | |

## db.modifiedPeptide.csv

This file has one row per modified peptide. The left half (`... Mod` columns) is the modified
form; the right half is the unmodified counterpart. The right half is **entirely empty** when no
unmodified form was identified, which turns the integer columns into nullable ints.

| column | dtype | kind / quote | meaning, quirks |
|---|---|---|---|
| Peptide Mod | str | str / never | |
| -10LgP Mod | float | fixed2 / never | |
| Mass Mod | float | fixed4 / never | |
| Length Mod | int | int / never | |
| PPM Mod | float | fixed1 / never | note the upper-case `PPM` |
| m/z Mod | float | fixed5 / never | |
| z Mod | int | int / never | |
| RT Mod | float | fixed4 / never | |
| Scan Mod | int | int / never | |
| Source File Mod | str | str / never | |
| Area Sample 1 Mod | float | sci / never | |
| #Feature Sample 1 Mod | int | int / never | |
| Found By Mod | str | str / never | |
| PTM | str | str / always | **different format**: `<AA><pos>:<name>:<true/false>`, e.g. `Q4:Deamidation (NQ):false` |
| AScore | str | str / always | `Q4:Deamidation (NQ):1000.00` (max 1000 here and 1001 in psms) |
| Accession | str | str / never | `:`-joined |
| RT shift | float | fixed4 / never | `RT Mod - RT`; empty with no counterpart |
| Peptide, -10LgP, Mass, Length, PPM, m/z, z, RT, Scan, Source File, Area Sample 1, #Feature Sample 1, Found By | | as above / never | unmodified counterpart; empty when none (`Int64` dtype for the int columns) |

## db.summary.html

`parse_summary_html` / `PeaksExport.params` reads the following tables:

* **Table 1** (statistics per sample and total) goes to `params["statistics"][sample]` with the
  keys `ms1_scans`, `ms2_scans`, `chimera_scans`, `features`, `psms` (= unique scan x feature),
  `identified_scans`, `identified_features`, `peptides` and `sequences` (both I/L-collapsed),
  `protein_groups`, `proteins_all`, `proteins_top`.
* **Table 2** (filters) gives `peptide_min_10lgp` (e.g. 28.36, set by the 1% peptide FDR),
  `protein_min_10lgp` (15), `protein_min_unique_peptides` (1) and `denovo_min_alc_pct` (50).
* **Table 3** gives `fdr_psm_pct`, `fdr_peptide_pct`, `fdr_protein_group_pct` and
  `denovo_only_spectra`. The FDRs are percentages with one decimal.
* **Table 5** (missed cleavages per sample) goes to `params["missed_cleavages"]`.
* **Table 6** (search parameters) gives the following keys:
  * `search_engine`
  * `precursor_tolerance` and `precursor_tolerance_unit` (15.0 ppm)
  * `fragment_tolerance` and `fragment_tolerance_unit` (0.5 Da)
  * `precursor_mass_type`, `enzyme`, `max_missed_cleavages`, `digest_mode`
  * `peptide_length_min` and `peptide_length_max`
  * `fixed_modifications` and `variable_modifications` (lists of `{name, delta}`)
  * `max_variable_ptm_per_peptide`, `database`, `taxon`, `searched_entries`
  * `deep_learning_boost`, `fdr_estimation`, `peaks_ptm`, `spider` (booleans)
  * `confident_aa_threshold_pct`

  The raw label/value strings are kept under `params["raw"]`. In the HTML, modification list
  items are `<br>&nbsp Name (+delta)`, with `&nbsp` lacking its semicolon.

## Files not parsed

These files are left as they are:

* `db.mzidentml_1_2_0.xml` (mzIdentML 1.2)
* `db.proteins.fasta`
* `db.peptides.html` and `db.proteins.html`
* `img/` (figures)
* `export.css`
* the MGF: `SCANS=` holds the PEAKS query index and `RAWSCANS=` the raw scan. There are more
  MGF entries than MS2 scans because chimeric scans are split per feature.
