# Feature columns and `<Sample>.db.features.csv`

Code: `peaksparity/features_export.py` (`build_features_table`, `feature_concordance`),
wired into `pipeline.run` as `stage_features_table`. Tests: `tests/test_features_export.py`.

## Data flow

1. `precursors.build_precursor_list` gives each feature-derived precursor entry its MS1
   `feature_id` and `area`. Featureless entries get `feature_id` 0 and no area.
2. `pipeline.attach_entries` carries these onto every search candidate as `feature_id` /
   `feature_area`, and `entry_rt` becomes the PSM RT (the feature apex).
3. `export_builder.build_psms` writes psms `Feature Id` / `Area`. `build_peptides` writes
   `Area <Sample>` / `#Feature`. `proteins.infer_proteins` writes protein `Area` from the
   distinct features of the unique peptides.
4. `features_export.build_features_table` writes `<Sample>.db.features.csv` from the run's
   feature table (`<run>/ms1_features.parquet`, or the `--features` parquet). It needs
   the psms table after `apply_protein_mapping`.

Steps 1 to 3 already worked before this change. In the cached Jurkat run, 71.2% of psms rows
had a Feature Id and every one of them had an Area. Two things were wrong or missing:

* **Peptide `Area` and `#Feature` summed PSM rows**, so a feature hit by several scans was
  counted several times. PEAKS sums over *distinct* features: we now get 3470/3472 areas and
  3472/3472 `#Feature` when we rebuild sample 4's peptides from its own psms. The row sum gave
  76%. On Jurkat the median peptide area ratio to PEAKS went from 1.34x to 1.02x.
* **`features.csv` was not written at all.**

## PEAKS rules for features.csv

These were checked by rebuilding each reference file from that sample's own psms.csv,
features.csv (as the feature map) and denovoOnly.csv. The test
`test_reproduces_peaks_features_from_its_psms` covers this.

* **Which rows.** Only features that have a DB PSM (any psms row with that Feature Id) or a de
  novo result with ALC >= 50 are listed. Rows are sorted by Feature Id. A feature with a DB PSM
  never also shows a de novo peptide. Every non-zero psms Feature Id appears, for PEAKS and for
  us.
* **Feature columns.** m/z, z, RT, RT Begin, RT End and Area come from the feature. Area is a
  Java float (jfloat). Quality is always `0.00`.
* **`-10LgP`** is the maximum over the feature's PSMs (100%).
* **`DB Peptide`** lists the peptides of the PSM rows that reach that maximum, `;`-joined.
  These are the I/L variants that `apply_protein_mapping` splits into separate rows. The set
  matches in 99.94 to 100% of rows. The remaining 15 of 3386 DB rows in S4 also list a
  different, non-I/L peptide, and we could not find the rule for that. PEAKS's order of the
  variants looks like hash order (in S4, ascending 17/48 and descending 27/48), so we sort them.
* **`Accession`** holds only the proteins of the lexicographically first variant, in psms
  order, `:`-joined. This matches 66/66 multi-variant features across S1 to S4, and 99.99% of
  all rows.
* **`ppm`** is the uncalibrated error of the feature m/z against the peptide mass,
  `((m/z - 1.007276) * z - M) / M * 1e6`. psms `ppm` is recalibrated, so in S4 it is
  8.0 to 8.1 ppm lower. With the exact peptide mass, 99.8 to 100% of rows land within 0.1 ppm,
  but only 74 to 81% print the same digit. PEAKS's value runs about 0.03 ppm high and we don't
  know why.
* **De novo rows** take the best-ALC de novo result of the feature: 99.9% peptide and 100% ALC
  agreement with denovoOnly.csv. They carry its own ppm, `-10LgP` is `0.00` and `DB Peptide`
  is empty. On DB rows, `ALC (%)` is `0.0`.

Rebuilt this way, 96.2% (S4) to 99.7% (S1) of the reference lines are byte-identical. The rest
differ only in the ppm digit or the variant order.

We have no de novo stage yet, so our features.csv has DB rows only. `build_features_table`
accepts `denovo_df` (a denovoOnly-style table with Feature Id, Peptide, ALC (%), Mass, ppm and
Accession) for when we do.

## Jurkat (sample 4) concordance

Run: the cached `pipeline_jurkat` intermediates, re-run through rescoring and export into
`/drive_a/drew/peaks-parity-work/pipeline_jurkat_features/`. Mokapot rescoring is not
bit-reproducible, so this run has 15,427 PSM rows against 15,398 in the earlier run. Numbers
are in `feature_concordance.json` and `concordance_features/`.

Matching follows `peaksparity.concordance --features`: same z, m/z within 10 ppm, RT within
0.5 min, one-to-one greedy.

| comparison | n PEAKS | n ours | matched | recall | precision |
|---|---|---|---|---|---|
| features.csv, all rows (CLI `--features`) | 23,540 | 6,646 | 3,965 | 16.8% | 59.7% |
| **features with a DB peptide, both sides** | 3,386 | 6,646 | 2,176 | **64.3%** | 32.7% |
| PEAKS DB features vs our whole feature map (203,862) | 3,386 | - | 2,770 | 81.8% | - |
| all PEAKS features.csv rows vs our whole map | 23,540 | - | 16,486 | 70.0% | - |

The all-rows recall is low because 20,154 of PEAKS's 23,540 rows are de novo only.

Area agreement on the 2,176 matched DB features:

* Spearman of log area: **0.973** (Pearson 0.974).
* Median log10(ours / PEAKS): **-0.003**, a ratio of 0.99. The IQR is 0.13 log10.
* 91.2% of areas are within 2x.
* Median m/z difference is +0.012 ppm and median RT difference is 0.000 min.
* The two DB peptides share an I/L-collapsed peptide in 94.5% of matched features.

Against our whole map (2,770 pairs), Spearman is 0.970 and the median log ratio is 0.000.

**PSMs with Feature Id > 0:** ours 71.9% (11,086 of 15,427 rows) vs PEAKS 71.9% (5,413 of
7,524).

Why the other 1,210 PEAKS DB features are missing from our features.csv:

| reason | features |
|---|---|
| not detected by our feature finder (no feature within 10 ppm / 0.5 min / same z) | 616 |
| detected, but we have no PSM at PEAKS's scans | 263 |
| detected, and we have the same peptide, but our PSM went to another (nearby or split) feature | 203 |
| detected, and we have the same peptide, but our PSM came from the featureless entry | 103 |
| detected, but our PSM at those scans is a different peptide | 25 |

Downstream areas, compared with PEAKS on shared keys:

| table | Spearman | median ratio before | median ratio after | within 2x |
|---|---|---|---|---|
| peptides (2,331 shared) | 0.853 | 1.34x | 1.02x | 77% |
| proteins (835 PEAKS top accessions) | 0.796 | 1.59x | 1.58x | 53% |

Protein Area was already computed over distinct features, so this change did not affect it.
It stays about 1.6x high, most likely because we report twice as many PSMs and therefore link
more features per protein.

## Open issues

* **psms `m/z`.** For feature-linked PSMs, PEAKS's psms m/z equals the feature m/z in 100% of
  cases. Ours is computed from Comet's `exp_mass` minus the isotope error. It is within about
  0.5 ppm of the feature m/z, but about 860 rows in the earlier run are shifted by ±1 isotope while still pointing at the
  feature.
* **Multi-file samples.** `features.csv` is keyed by (Source File, Feature Id). However,
  `proteins.infer_proteins` keys feature areas by Feature Id alone, and feature ids restart in
  every run.
* **Unknown rules.** The `ppm` last digit, the order of I/L variants in `DB Peptide`, and the 15
  non-I/L extra peptides (see above).
