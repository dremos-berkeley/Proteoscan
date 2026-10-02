# PEAKS protein inference and FDR: rules inferred from the reference exports

Stages 6-7 (`proteoscan/fdr.py`, `proteoscan/proteins.py`). The published rules are in
docs/PEAKS_ALGORITHMS.md, sections 2.2 and 4. This file records what the four
DrewR_Adam_91526 exports (PEAKS Studio 12/13, human SwissProt, semi-specific trypsin,
2 missed cleavages) show, and how each rule was checked. The PEAKS binary was never used.

Harness: `python scripts/eval_proteins.py`. For each sample it feeds PEAKS's own
`db.peptides.csv` into `infer_proteins`, together with `db.psms.csv`, `<sample>.db.features.csv`
and the full-precision PSM scores from `db.mzidentml_1_2_0.xml`. It then diffs the output
against `db.proteins.csv` and `db.protein-peptides.csv`. Every field is compared as the exact
string PEAKS would print. Results are written to
`/drive_a/drew/proteoscan-work/proteins/eval/`. Reference tests are in `tests/test_proteins.py`
and `tests/test_fdr.py`.

## Final agreement

Mapping mode `accession` takes each peptide's protein list from the peptide `Accession` column
and reads positions from the FASTA. Mapping mode `fasta` maps peptides by sequence, using the
digestion rules below.

| sample | groups identical | proteins (ref/test) | Top, -10lgP (exact float), #Peptides, #Unique, Coverage, #Spec, Area, Average Mass, PTM, Gene, Description | protein-peptides rows | Unique, Start/End, Extended Peptide, Mass, Area, #Feature | Protein Group number exact | pair order |
|---|---|---|---|---|---|---|---|
| 1 | 7/7 | 20/20 | 100 % each | 114/114 | 100 % | 85 % | 0.857 |
| 2 | 26/26 | 44/44 | 100 % | 397/397 | 100 % | 75 % | 0.963 |
| 3 | 16/16 (fasta: 15/16) | 34/34 (fasta: 34/35) | 100 % | 250/250 | 100 % | 68 % | 0.950 |
| 4 (Jurkat) | 990/990 (fasta: 988/990) | 1238/1238 (fasta: 1238/1240) | 100 % | 4463/4463 | 100 % | 3.7 % | 0.881 |

The same numbers hold in `fasta` mode except where marked. The extra proteins found in
`fasta` mode are P0DOX7 IGK_HUMAN (sample 3), and Q13765 NACA_HUMAN and P42167 LAP2B_HUMAN
(sample 4). Each is a database-release difference (see "FASTA" below).

If the CSV's 2-decimal peptide -10lgP is used instead of the mzIdentML scores (`--no-mzid`),
every field above is unchanged except the protein -10lgP. The rounded input limits it to
74-91 % agreement at 2 decimals and 0 % as an exact float.

Only one thing is unresolved: the order of Protein Group numbers among groups whose first
protein has the same number of peptide rows (see "Group numbering").

## FASTA

PEAKS's database, "Homo_sapiens" with 20597 entries, is an older UniProt release than the
local `human_sp.fasta`, which has 20431 entries. For evaluation,
`/drive_a/drew/proteoscan-work/proteins/human_sp_peaks_merged.fasta` is built by
`scripts/eval_proteins.py:build_merged_fasta`. It is the local FASTA with every entry
replaced or added from the union of the four `db.proteins.fasta` files, which hold PEAKS's
own headers and sequences for the reported proteins.

The differences found:

**Accessions in the exports that are missing locally**
- A0A0G2JMH3|A0A0G2JMH3_HUMAN: appears only in peptide Accession lists (TTILYK). Its sequence
  is unknown and it is never reported.
- A0A140TA62 (tr) and O60361 NDK8_HUMAN: restored from `db.proteins.fasta`.

**Entries renamed**
- P38646 GRP75_HUMAN is now HSPA9_HUMAN.
- Q9HCU5 PREB_HUMAN is now SEC12_HUMAN.
- Q8IVW1 ARL17_HUMAN is now AL17A_HUMAN.

**Sequences changed**
- C2TA_HUMAN, 4F2_HUMAN, ATD3A_HUMAN, PAP1L_HUMAN, CRBG1_HUMAN, ABCAD_HUMAN, WDR36_HUMAN.
- SELH_HUMAN differs only in U44 versus X44. PEAKS stores selenocysteine as X: no exported
  sequence contains U. `ProteinDB` therefore converts U to X.

**Headers changed**
- 81 headers changed. The changes include descriptions and genes; for example, H4_HUMAN has
  GN=H4C16 in PEAKS and H4C1 locally.

**Local entries that PEAKS did not have, or had with a different sequence**
These peptides are absent from PEAKS's Accession lists although the local sequence
contains them as valid digestion products:
- P0DOX7 IGK_HUMAN
- Q13765 NACA_HUMAN
- P42167 LAP2B_HUMAN
- P63092 GNAS2_HUMAN (the protein is filtered out in either case)

## Peptide to protein mapping

- **Semi-specific trypsin, K/R not before P, at most 2 missed cleavages.**
  - Evidence (all four samples): exact substring mapping reproduces PEAKS's Accession lists
    for 3458/3472 peptides in sample 4.
  - Requiring at least one enzymatic terminus with the "not before P" rule raises this to
    3464/3472 (47/47, 261/261 and 168/169 in samples 1-3).
  - All remaining differences are database-release differences.
  - Examples of proteins that are no longer assigned: AVFPSIVGR to ACTBM (followed by P),
    TTLLYK to AAKB1/2, VDLQAK to CA213, and VTIAQGGVLPNIQ to H2AX/H2A2B/H2A1A.
  - The missed-cleavage histogram of db.summary.html is reproduced exactly with the "not
    before P" rule: 36/11, 227/32/2, 146/21/2 and 3280/183/9.
- **A protein N-terminus, or position 2 after an initiator Met, is an enzymatic terminus.**
  Evidence: PQ(+0.98)LGG...AAG in TF7L1_HUMAN starts at residue 2 after M, and its C-terminus
  is not tryptic.
- **I and L are equal for domination and uniqueness; protein-peptides rows list verbatim
  matches only.**
  - PEAKS emits one PSM and peptide row per database I/L variant, with the same scan and
    score. For example, sample 4 has AILFLPR (CCL21) and ALLFIPR (HS90B), and VVVLTK
    (RHG35) and VVVITK (HS90B).
  - Only the verbatim variant is listed under each protein. Yet CCL21 and RHG35 are
    dominated members of HS90B's group 8, and AILFLPR is `Unique=true`.
  - So grouping works on I/L-collapsed "peptide families".
  - If a caller's peptide table is not I/L-closed (e.g. our own search), a protein that
    matches a family only through I/L gets the input row at the I/L position. This never
    triggers on PEAKS input.
- **Start/End use the first valid occurrence.** Coverage uses only that occurrence.
  - Evidence: UBC_HUMAN repeats; 121 multi-occurrence rows.
- **Extended Peptide omits the flank and its dot at a protein terminus.** Examples:
  `MAGSGAGVR.C` and `R.TLYGFGG`.

## Protein -10lgP

- **Score = Σ_k s_k / k over distinct residue sequences, not peptidoforms.**
  - The same sequence with different PTM sets counts once, at its best score. This
    contradicts the published "different PTMs count separately".
  - Evidence, with scores from mzIdentML: 1336/1336 protein rows (all samples) match to
    within 1e-4 relative under the sequence rule. Under the published PTM-set rule, only
    1143/1336 match. For example, K2C1_HUMAN in sample 2 scores 173.04 under the sequence
    rule and 181.57 under the PTM-set rule.
- **Float32 arithmetic: `t = (float)(t + (float)(s_k / k))`, with s_k sorted descending.**
  - This matches the Java float printed value bit for bit for 1336/1336 rows. Double
    precision matches 1098/1336.
- **Peptide scores must be full precision.** The CSVs round -10lgP to 2 decimals; the
  mzIdentML `PEAKS:peptideScore` has the float value. All peptides' mzIdentML maximum rounds
  to the CSV value (0 disagreements).
  - The peptide score is the maximum over its PSMs.

## Grouping, Top and numbering

- **Domination** (published): X dominates Y if every family of Y with score ≥ T (T = 15) is
  in X. In these data every peptide scores above 15, so T is not exercised.
- **Greedy order: most supporting peptide rows (peptidoforms) first.**
  - Evidence: in sample 2, K2C78_HUMAN is dominated by both K2C3 and K2C4.
  - PEAKS puts it in K2C3's group, which has 5 rows and -10lgP 77.5. K2C4 has 4 rows and a
    higher -10lgP of 80.4, so a score-ordered greedy would pick K2C4.
  - With row-count order, all pre-filter groups in the mzIdentML (including filtered
    proteins) are reproduced in `accession` mode:
    - 7/7, 31/31 and 19/19 in samples 1-3;
    - 997/998 in sample 4. The one miss is the ARF group, which contains A0A0G2JMH3 and
      ARL17, both missing from the FASTA.
  - In `fasta` mode, 994/998 match in sample 4. The misses are NACA, LAP2B, and a
    filtered G-protein group that GNAS2 changes.
- **Top = same family set as the group's first protein (mutual domination).**
  - Evidence: 1336/1336 rows, e.g. the ACTB/ACTG and RSMB/RSMN pairs.
  - "Most unique peptides in the group" gives only 1216/1238 in sample 4.
- **A dominated protein joins the first group, in processing order, whose leader dominates
  it.** An example is H90B3_HUMAN, which has 16 rows, in HS90B's group 8.
- **Group numbering.**
  - `Protein Group` = the mzIdentML `PROTEINAMBIGUITYGROUP_n` + 1.
  - Numbers are assigned before the protein filter, in greedy creation order: by the first
    protein's number of peptide rows, descending. The reference numbering violates this key
    0 times in 4 samples.
  - Groups whose proteins are all filtered (all `passThreshold="false"` in mzIdentML) leave
    gaps. Examples: group 19 in sample 4 is the tubulin TBA1B/TBA1C group, with no unique
    peptides; groups 11 and 15-17 in sample 2.
  - The order among ties looks hash-based. No tested secondary key does better than chance,
    over 114729 within-tie pairs:
    - -10lgP, #Spec, area, coverage, accession, sequence length or mass;
    - first scan or RT;
    - best or summed peptide score;
    - mzIdentML DBSequence order;
    - Java `String.hashCode` bucket order for 16-32768 buckets.
  - Group size is weakly informative, at 64 %.
  - We break ties by -10lgP, then accession. This is the only unresolved discrepancy: the
    exact group number matches for 85 %, 75 %, 68 % and 3.7 % of proteins in samples 1-4.
    Sample 4 has 960 adjacent ties among 990 groups.

## Unique, #Peptides, #Unique, #Spec, coverage, PTM

All of the following reproduce 100 % of rows in all samples.

- **Unique:** a family is unique if all proteins it maps to, before filtering, are in one
  group.
- **#Peptides / #Unique:** distinct (Start, End) among the protein's rows, or among its
  unique rows. This is the "same start/end counted once regardless of PTM" rule.
- **#Spec:** the number of `db.psms.csv` rows for all supporting peptidoforms, unique or not.
  This includes PSMs below the peptide threshold. Counting only unique peptides matches
  1078/1238 in sample 4.
- **Coverage(%) and Coverage(%) <sample>:** the union of [Start, End] divided by the sequence
  length. The two columns are identical for a single sample.
- **PTM:** the sorted unique PTM names of all supporting rows, joined with "; ". Fixed
  Carbamidomethylation is included when a supporting peptide carries it.

## Area

- **Protein Area:** the float32 sum over the *distinct* MS1 features of the unique
  supporting peptidoforms. The feature areas come from `<sample>.db.features.csv`, where
  Area is a Java float.
  - Summing per peptide instead double-counts shared features: 1232/1238 against 1237/1238.
  - It is printed with 3 significant digits rounding **HALF_UP**. With float32 and HALF_UP,
    1165/1165 non-empty areas match; any other combination gets 1164.
  - The single tie is ERF3A_HUMAN: 1065000.0f is printed "1.07E6". `peaks_export`'s
    formatter rounds HALF_EVEN, so `infer_proteins` returns the next double above an exact
    tie.
  - RALY_HUMAN needs the float32 sum: 671499.98 in double becomes 671500.0f, printed 6.72E5.
  - Area is empty when no unique peptide has a feature.
- **Protein-peptides `Area <sample>` / `#Feature`:** counts only features that carry a PSM
  with -10lgP ≥ 20.
  - Evidence: DVDAAYM(+15.99)NK in sample 3 has three features. Its `db.peptides.csv` row
    counts all 3 (6.92E5), while its protein-peptides rows count 2 (6.31E5).
  - The excluded feature's only PSM scores 19.73; the others score 20.75-27.18. Over all
    samples, the consistent threshold lies in (19.73, 20.02], so it is 20 (p = 0.01).
  - The same filter is applied to the protein Area; no reference case distinguishes it.

## Masses, Description, Gene

- **Average Mass:** the sum of classic average residue masses (A 71.0788, R 156.1875, ...;
  fitted to 1336 proteins), plus **18.01056** (monoisotopic water with 5 decimals), in
  double, cast to float32.
  - This matches 1260/1260 distinct proteins. With average water (18.01528), 247 match.
  - X is the mean of the 20 residues (SELH_HUMAN).
- **Protein-peptides Mass:** the exact decimal sum of 5-decimal monoisotopic residue masses
  (A 71.03711 ...), plus 18.01056, plus exact modification deltas (15.994915, 0.984016,
  57.021464), rounded **HALF_DOWN** to 4 decimals.
  - This matches 3932/3932 distinct peptides.
  - `db.peptides.csv` prints a float32 search mass instead, so the two files disagree in the
    4th decimal for about 30 % of peptides.
- **Description** is the FASTA header after the first token. **Gene** is the `GN=` value, or
  "" when absent. Both match 1336/1336, using PEAKS's headers.

## Filter

- -10lgP ≥ 15 and #Unique ≥ 1, applied after grouping and uniqueness.
- Filtered proteins stay visible in mzIdentML with `passThreshold="false"`.
  - Sample 4 has 1326 hypotheses, of which 1238 pass.
  - Our filter agrees with mzIdentML's `passThreshold` for every hypothesis in all samples
    in `accession` mode.

## FDR (stage 6)

The published rule (docs/PEAKS_ALGORITHMS.md §2.2) is FDR = D/T at score ≥ t. The exports
contain no decoys, so the threshold itself cannot be checked. The counting units were checked
against db.summary.html (`tests/test_fdr.py::test_reference_counting_units`):

| summary count | equals |
|---|---|
| #PSMs (e.g. 7424) | distinct (Scan, Feature Id) in db.psms.csv: one per query (spectrum × precursor feature); I/L variants and the per-feature duplication collapse |
| #Peptides (3418) | distinct modified sequences with I = L |
| #Sequences (3230) | distinct stripped sequences with I = L; this is the unit of "FDR (Peptide Sequences)" |

`threshold_at_fdr(scores, is_decoy, alpha)` returns the lowest observed score t with
D(≥t)/T(≥t) ≤ alpha. Ties are accepted together, and there is no monotonisation, following
docs/PEAKS_ALGORITHMS.md §2.2 [I]. `q_values` gives the conventional monotone q-value.

The three level-specific functions are:
- `psm_fdr` keeps the best hit per query.
- `peptide_fdr` scores each I/L-collapsed stripped sequence by its best PSM.
- `protein_group_fdr` uses each group's Top protein.

The peptide thresholds reported by PEAKS are 28.36, 26.79, 24.08 and 23.02 for samples 1-4,
each at 1.0 % peptide FDR. They are a target for the rescoring stage, not for this module.
