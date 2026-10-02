# PEAKS PTM tables: db.ptmProfile.csv and db.modifiedPeptide.csv

`proteoscan/ptm_tables.py` builds both tables from PEAKS-format psms, peptides and
protein-peptides tables. Every rule below was inferred from the four DrewR_Adam_91526
reference exports (PEAKS Studio 12/13; Carbamidomethylation fixed, Deamidation (NQ) and
Oxidation (M) variable, at most 3 variable PTMs per peptide). The inputs were PEAKS's own tables,
read with `peaks_export.read_export`. The PEAKS binary was never used. Tests are in
`tests/test_ptm_tables.py`.

## Final agreement

Each output was written with `peaks_export.write_table` and compared with the reference file as
text.

| sample | ptmProfile | modifiedPeptide rows | modifiedPeptide lines byte-identical | differing cells |
|---|---|---|---|---|
| 1 | byte-identical | 3/3 | 3/3 | none |
| 2 | byte-identical | 72/72 | 69/72 | 3 `RT shift` cells, each off by 0.0001 |
| 3 | byte-identical | 27/27 | 27/27 | none |
| 4 (Jurkat) | byte-identical | 741/741 | 736/741 | 5 `RT shift` cells, each off by 0.0001 |

modifiedPeptide needs two inputs that the CSV exports do not contain:

- **Full-precision PSM scores**, passed as `scores=`. For the references they come from
  `db.mzidentml_1_2_0.xml` through `peaks_reference_scores`.
- **Isolation-window target m/z per MS2 scan**, passed as `isolation_mz=`. For the references
  they come from the ThermoRawFileParser mzML through `read_isolation_targets`.

Without `scores`, the representative PSM and the row order fall back to the 2-decimal
`-10LgP`. With the isolation targets supplied, the table then has:

- every row present, with the right unmodified counterpart;
- rows identical by content: 3/3, 68/72, 27/27 and 733/741;
- rows in the right position: 3/3, 50/72, 25/27 and 538/741.

Without the isolation targets, the `m/z` columns are wrong in most rows.

The only unresolved discrepancy is `RT shift`: 8 of 187 cells differ in the last digit (see
"RT shift"). ptmProfile needs only psms, peptides and protein-peptides.

## Terms

- **Site position** is 1-based within the peptide. PTM strings name a site as `<residue><position>`,
  e.g. `M5`. Sites are listed in position order.
- **Site names** come from the peptide's `PTM` column: one `;`-separated entry per modified site,
  in the same order as the mass deltas in the sequence.
- **Variable vs fixed.** A PTM is variable if its name ever appears in an `AScore` string.
  Fixed Carbamidomethylation never does. The `variable_ptms=` argument overrides this.
- **Confident site.** A PSM's modification at peptide position p is confident when both hold:
  - its `Positional Confidence` token p is ≥ 2;
  - its AScore for p is ≥ 15.

  Fixed modifications are always confident.
  - Evidence: 1006 variable sites in `modifiedPeptide.PTM`, over all samples, with 0 errors.
  - Among sites with positional confidence ≥ 2, the largest AScore flagged `false` is 14.63
    and the smallest flagged `true` is 15.62.
  - Among sites with AScore ≥ 15, the largest positional confidence flagged `false` is 1 and
    the smallest flagged `true` is 2. Of the false ones, 189 have 0 and 43 have 1.
  - Neither criterion works alone:
    - positional confidence ≥ 2 alone gives 58 false positives;
    - AScore ≥ 15 alone gives 232 false positives.
  - The thresholds are therefore AScore ∈ (14.63, 15.62] and positional confidence ∈ (1, 2].
    We use 15 and 2.
  - The PEAKS manual describes the confident-PTM filter as an AScore threshold plus an
    ion-intensity threshold (docs/PEAKS_ALGORITHMS.md §3, "≥ 2 %"). The per-residue positional
    confidence behaves exactly like the second criterion here. That the two are the same
    quantity is a guess.
  - The one Carbamidomethylation site in modifiedPeptide (sample 2) is `true`, even though no
    AScore is reported for it.

## db.modifiedPeptide.csv

### Rows

There is one row for each `db.peptides.csv` peptide that carries at least one variable PTM. The
counts are 3, 72, 27 and 741 rows, all of which are reproduced. Peptides whose only PTM is the
fixed C(+57.02) are left out.

### The "Mod" half

- `-10LgP Mod`, `Mass Mod`, `Length Mod`, `Area <S> Mod`, `#Feature <S> Mod`, `Found By Mod`
  and `Accession` are copied from the `db.peptides.csv` row. This matches 843/843 rows.
- `PPM Mod`, `z Mod`, `RT Mod`, `Scan Mod`, `Source File Mod` and `AScore` come from a
  **different representative PSM** than db.peptides.csv uses:
  - db.peptides.csv takes the *first* PSM (in db.psms.csv row order) among those with the
    best full-precision score.
  - modifiedPeptide takes the *last* one.
  - Evidence: 54 peptides have exact full-score ties (3 in sample 2, 51 in sample 4):
    - the last-row rule matches the modifiedPeptide choice in 54/54;
    - the first-row rule matches db.peptides.csv in 54/54;
    - "last by scan number" gets only 42/54 and "last by MGF query index" 49/54.
  - Exact ties are common because PEAKS scores are quantised. Every PSM at a sample's ceiling
    has the same float: 44.80417 × 1534 in sample 4. Values below the ceiling repeat as well,
    for example 43.59689 and 39.2851.
  - The CSV's 2-decimal scores therefore cannot tell exact ties from near ties (7 cases each
    in samples 2 and 4). Hence the `scores=` argument.
  - Without it, a 2-decimal tie is trusted only at the sample's score ceiling (the maximum
    `-10LgP` in psms). Otherwise the db.peptides.csv PSM is used. This leaves 0/1/0/3
    `Scan Mod` errors.
- `AScore` is the representative PSM's AScore string with values capped at 1000.00; psms.csv
  uses 1001.00 for an unambiguous site. It matches 843/843 rows. In 5 rows it differs from
  db.peptides.csv, because the representative differs.
- `PTM` is `<AA><pos>:<name>:<true|false>` for every modified site, fixed ones included, joined
  with `;`. The flag is the confident-site rule applied to the representative PSM alone, not to
  any PSM.
- `m/z Mod` and `m/z` are **the MS2 scan's isolation-window target m/z**. This is the m/z in the
  Thermo filter string (`ITMS + c NSI t d Full ms2 702.0225@hcd30.00`), a float32, printed to 5
  decimals.
  - It is neither the PSM `m/z` (the feature or monoisotope-corrected PEPMASS) nor `raw m/z`
    (the Thermo "Monoisotopic M/Z").
  - Evidence: 843/843 Mod values and 187/187 unmodified values equal the mzML
    `isolation window target m/z` of `Scan Mod` or `Scan`.
  - By comparison, only 219/843 equal the PSM `m/z`. Also, `raw m/z` is 0 in psms for most
    rows, and when non-zero it is often one isotope away (e.g. 701.68842 vs 702.02252).
  - Rounding is **HALF_EVEN** on the float value. The one exact tie, float 848.765625, is
    printed `848.76562`. `peaks_export`'s `fixed5` would give `848.76563`, so the builder
    pre-rounds this column.
  - Without `isolation_mz=`, the representative PSM's `raw m/z` is used, or `m/z` when that is
    0. This value is not what PEAKS writes.

### The unmodified half

- `Peptide` is the modified peptide with its **variable** modifications removed and fixed ones
  kept, e.g. `NM(+15.99)QDM(+15.99)VEDYR` → `NMQDMVEDYR`. The half is filled only when that
  exact string, with I and L kept as they are, is a db.peptides.csv row. Otherwise it is empty.
  - Evidence: 741/741 decisions agree (3, 26, 8 and 150 rows have a counterpart). No I/L-only
    counterpart occurs.
  - A partially modified form such as `NMQDM(+15.99)VEDYR` never counts as the counterpart.
- `-10LgP`, `Mass`, `Length`, `PPM`, `z`, `RT`, `Scan`, `Source File`, `Area <S>` and
  `Found By` are copied from the counterpart's db.peptides.csv row. This matches 187/187.
  - Unlike the Mod half, this half uses the db.peptides.csv (first) representative.
- `#Feature <S>` = the number of distinct `Feature Id`s over all psms rows of the counterpart,
  **including 0**. PSMs without a feature count as one extra "feature".
  - Evidence: 187/187 rows. The db.peptides.csv `#Feature`, which excludes 0, matches only
    119/187.
  - The Mod half does *not* do this: its value equals db.peptides.csv in 843/843 rows. The
    rule with 0 included would give 509/843.
- `RT shift` = `RT Mod − RT` (see below). It is empty when there is no counterpart.
- The integer columns become nullable (`Int64`) because of the empty rows.

### RT shift

This is the difference of the two full-precision feature RTs, rounded to 4 decimals.

- From the 4-decimal CSV RTs, 179/187 cells are exact. Samples 2 and 4 have 3 and 5 cells that
  are 0.0001 off.
- Other sources of RT do not do better:
  - The MGF `RTINSECONDS`, whose 3 decimals of seconds are more precise, gives 177/187. Every
    residual there is an exact ...5 tie at the fifth decimal.
  - Intersecting both precision intervals leaves 29 of 187 cases undetermined.
  - mzML scan start times are not the source: at most 90/187 match (nearest scan of any MS
    level). Feature apex RTs are interpolated.
- The true RTs are therefore not recoverable from the exports. This is the only unresolved
  discrepancy. Pass `rts=` (full-precision minutes per psms row) when the pipeline has them.

### Row order

Rows are sorted by the following keys:

1. The peptide's best full-precision score, descending. There are 0 violations in 4 samples.
2. The stripped sequence as a Java/ASCII string, **descending**. This holds for 361/361
   adjacent full-score-tied pairs with different sequences.
3. For the same stripped sequence, residue by residue, descending, with this order at a
   residue: deamidated N or Q < unmodified < any other variable modification (Oxidation).
   - This is the same as comparing strings in which N(+.98)→D, Q(+.98)→E and M(+15.99)→M~.
   - It holds for 7/7 cases. For example, `M(+15.99)SVQPTVSLGGFEITPPVVLR` comes before
     `MSVQ(+0.98)PTVSLGGFEITPPVVLR`, and `GNLGAGNGN(+0.98)LQGPR` before
     `GNLGAGN(+0.98)GNLQGPR`.
   - Plain descending order of the modified string explains 6/7. With only 7 cases, this key
     has limited evidence.

db.peptides.csv itself is in ascending modified-string order. Without `scores`, keys 1-2 use
the 2-decimal `-10LgP`. That misplaces rows among 2-decimal ties: 538/741 rows are in place in
sample 4.

## db.ptmProfile.csv

This file is the 7-column PTM summary, identical to Table 4 of db.summary.html. It is **not**
the per-site "PTM profile" view (Protein Position, AScore, Ion Intensity, Modified/Unmodified
area, "Modified Peptide Sequence Window") described in docs/PEAKS_ALGORITHMS.md §7.7. None of
those columns appear in this PEAKS version's CSV, so no sequence window, area ratio or per-site
AScore/ion-intensity column has to be reproduced here.

| column | rule | evidence |
|---|---|---|
| rows | one per PTM name found in any db.peptides.csv `PTM`; a PTM with no peptide is left out (sample 1 has no Oxidation row) | 11/11 rows |
| order | `# peptides` descending (ties: variable before fixed, then name; tie order untested) | 4/4 files |
| `Name` | PEAKS PTM name | |
| `Mass (Da)` | the 2-decimal delta from the peptide notation, printed `dec2` (`.98`, `15.99`, `57.02`) | |
| `Position` | the residues in the name's trailing parentheses (`NQ`, `M`); otherwise the sorted residues that carry it (`C`) | |
| `# peptides` | distinct **stripped** sequences (I ≠ L) in db.peptides.csv carrying the PTM | 11/11. Counting rows gives 52 vs 40 (sample 2 Oxidation) and 570 vs 533; I/L-collapsed gives 529 vs 533; protein-peptides Top+Unique gives 36 vs 40 |
| `# proteins` | distinct accessions among db.protein-peptides.csv rows with **Top = true and Unique = true** whose peptide carries the PTM | 11/11. All rows give 17 vs 11; Top only gives 13 vs 11; Unique only gives 318 vs 293 (sample 4 Oxidation) |
| `# modification sites` | distinct (accession, `Start` + position − 1) over the same Top+Unique rows | 11/11 |
| `# confident modification sites` | those sites where **any** db.psms.csv PSM (below-threshold PSMs included) of a supporting peptidoform is confident at that position; fixed PTMs are always confident | 11/11 |

Notes on the site columns:

- **Positions use `Start`, which is the first valid occurrence** (docs/PROTEIN_INFERENCE.md).
  Counting every occurrence of the peptide in the merged FASTA
  (`/drive_a/drew/proteoscan-work/fasta/human_sp_peaks_merged.fasta`) overcounts:
  - sample 1 Deamidation: 15 vs 5 (the UBC/UBB repeats);
  - sample 2 Deamidation: 80 vs 40;
  - sample 4 Oxidation: 583 vs 582.

  `Start` equals the first verbatim FASTA occurrence for every Top+Unique row. The `fasta`
  argument is read only when protein-peptides has no `Start` column.
- **Confident sites use any PSM, not the modifiedPeptide flag.** Using the modifiedPeptide
  representative-PSM flags instead gives:
  - sample 1: 0 vs 5 (all three flags there are `false`);
  - sample 2 Deamidation: 18 vs 20;
  - sample 3 Deamidation: 1 vs 2;
  - sample 4: 451 vs 471 and 83 vs 87.

  Treating Carbamidomethylation by the positional-confidence rule would give 0/1 and 3/5.
- `proteins` (db.proteins.csv) is not needed. Protein-peptides already holds only reported
  proteins, with `Top` and `Unique`.
- `# confident modification sites` is the first number in summary Table 4 and
  `# modification sites` is the one in parentheses, e.g. `471 (582)`. The Table 4 header prints
  only "# modification sites". docs/PEAKS_EXPORT_FORMAT.md describes the parentheses the other
  way round.

## Inputs for the pipeline

```python
from proteoscan.ptm_tables import build_ptm_profile, build_modified_peptides, read_isolation_targets

ptm_profile = build_ptm_profile(psms_tbl, peptides_tbl, proteins_tbl, protein_peptides_tbl,
                                cfg.database, sample_name)
iso = {}
for mz in mzmls:                      # key by the Source File name written in psms
    raw = Path(mz).name.replace(".mzML", ".raw")
    iso.update({(raw, s): v for s, v in read_isolation_targets(mz).items()})
modified = build_modified_peptides(psms_tbl, peptides_tbl, proteins_tbl, protein_peptides_tbl,
                                   sample_name, scores=full_scores, isolation_mz=iso)
```

- `full_scores` must be aligned with the rows of `psms_tbl`. Note that `build_psms` sorts its
  output, so the unrounded `minus10lgp` must be carried through the same order.
- "Last tied PSM in psms.csv order" is defined in terms of PEAKS's own row order. Our psms table
  is ordered differently, so for our own searches the choice among exact ties is only
  PEAKS-like.

## Reference helpers

- `read_isolation_targets(mzML)` returns {scan: isolation target} by streaming the mzML text.
  It takes about 5 s for a 1.3 GB file.
- `peaks_reference_scores(export_dir, psms, cache_dir)` gives full-precision scores for the
  rows of a PEAKS psms table:
  - `spectrumID="index=N"` in the mzIdentML is MGF entry `SCANS=N+1`;
  - that entry's `RAWSCANS`, `PEPMASS` and `CHARGE`, together with the peptide, give the psms
    key (Scan, m/z, z, Peptide).
  - This matches 100 % of mzIdentML items. PSMs below the peptide threshold have no item and
    get NaN, but they never tie with the best PSM.
  - The test caches the join in `/drive_a/drew/proteoscan-work/ptm_tables/cache/`.
