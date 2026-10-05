# De novo agreement features for rescoring (CAA, tag length)

> **Superseded (2026-10):** the validation and score-imitation numbers below were produced on
> candidates of the propionamide search, which carried a hidden static +57 on every Cys
> (docs/STATUS.md "Corrections", P1). The feature definitions stand; the numbers need re-running on
> corrected candidates.

PEAKS DB scores every candidate partly by how well it agrees with the de novo sequence of the
spectrum (docs/PEAKS_ALGORITHMS.md 2.1, feature 1 "CAA score"; section 6). This file documents
our imitation: Casanovo calls for every (scan, precursor) query of a run, per-candidate
agreement features, their validation against the PEAKS propionamide rerun export, and their
effect on the score-imitation model (scripts/train_score_model.py).

Code: `proteoscan/denovo_features.py` (features + sharded Casanovo runner),
`scripts/add_denovo_features.py` (join onto a candidate table),
`scripts/eval_denovo_features.py` (validation against a PEAKS export),
`tests/test_denovo_features.py`.
Working data: `/drive_a/drew/proteoscan-work/denovo_features/`.

## Features

For a candidate c and the de novo call d of the same query (precursor-consistent pick from 5
beams, docs/DENOVO.md), with confidence threshold 0.3 and mass tolerance 0.02 Da:

| column | definition |
|---|---|
| `caa` | confident residues of d (Casanovo per-residue score >= 0.3) aligned to the same residue of c. Low-confidence residues act as PEAKS's mass segments: they are never counted and break tags, but their mass keeps the alignment. A residue is aligned when its start and end cumulative masses coincide (0.02 Da) with a residue boundary of c and the residue identity is equal. Identity = letter with I = L plus the modification delta rounded to 0.01 Da; N-terminal mods are folded into residue 1. Masses are taken from the N terminus (prefix masses) before a split point and from the C terminus (suffix masses) after it. The split that maximises the count is used, and each residue of c is used at most once. This way a mass offset inside the peptide (a wrong precursor, an isotope error, a residue outside Casanovo's vocabulary such as propionamide-Cys) costs only the residues around it. |
| `caa_pct` | 100 * caa / len(c) (PEAKS `CAA (%)` scale) |
| `tag_len` | longest run of consecutive aligned confident residues of d |
| `max_caa_any_beam` | max `caa` over all (up to 5) beam candidates of the query |
| `denovo_conf_mean` | mean de novo confidence of the aligned residues (0 if none) |
| `denovo_exact` | 1 if the chosen call equals c residue by residue (I = L, mods as rounded mass) |

Queries without a Casanovo call (charge > 4, too few peaks) get NaN. The training script and
rescoring fill NaN with 0.

Performance: `denovo_features` takes 1.6 M candidates x 330 k queries in about 30 s: 27 s to
parse the unique peptides, under 1 s for the numba kernel (16 threads).

API:

```python
from proteoscan.denovo_features import run_denovo_queries, add_denovo_features
# queries: DataFrame(id=query_id, mz, charge) of the pipeline's precursors.parquet;
# the expanded mzML has one spectrum per query, native id scan=<query_id>
dn = run_denovo_queries(expanded_mzml, queries, ckpt, out_dir, stem, n_shards=2, gpus=("2", "3"),
                        precision="16-mixed")          # every beam candidate, `picked` flagged
cands = add_denovo_features(cands, dn.rename(columns={"id": "query_id"}))
```

## De novo runs (Casanovo ft5runs, beam 5, top 5, precursor-consistent pick)

Model `/drive_a/drew/proteoscan-work/denovo_finetune/models/casanovo_iontrap_hcd_ft5runs_step15000.ckpt`
(copy: `models/casanovo_iontrap_hcd_ft5runs_step15000.ckpt`), `n_beams 5`, `top_match 5`,
`max_peptide_len 50`, `max_charge 4` (Casanovo skips higher charges), fp16 autocast
(`precision: 16-mixed`), batch 512. Input: the pipeline's expanded mzML (one spectrum per
(scan, precursor) query, native id `scan=<query_id>`), with m/z and charge from
`precursors.parquet`, split into 2 MGF shards run in parallel on GPUs 2 and 3.

| run | queries sequenced | rows (all beams) | picks within 15 ppm | wall time |
|---|---|---|---|---|
| OT_4, queries with candidates | 330,223 / 330,480 | 1,650,574 | 12.3 % | 56 min GPU + 50 min Casanovo mzTab writing (before the fix below) |
| OT_4, other z <= 4 queries (no candidates) | 98,493 / 98,720 | 492,462 | 1.6 % | 24.5 min |
| OT_1, queries with candidates | 140,754 / 141,628 | 701,018 | 11.6 % | 26.5 min |
| OT_2, queries with candidates | 159,047 / 159,181 | 795,017 | 9.5 % | 29.9 min |
| OT_3, queries with candidates | 94,622 / 94,754 | 472,863 | 5.9 % | 25.5 min |

Throughput is about 45 spectra/s per A6000 at beam 5 (fp16); the GPU is the bottleneck, and
4 processes on one GPU were no faster than 1. Choices, measured on 2,000 OT_4 spectra:
- fp16 autocast: 51 s against 85 s in fp32. 95 % of picked calls are identical (98 % in the
  top-scoring half) and the precursor-consistent fraction is unchanged (14.4 vs 14.6 %).
- bf16: 41 s, but only 73 % identical. Not used.
- TF32: 64 s, 96 % identical.
- batch 2048: no faster.

Casanovo 5.2.1 parses every output PSM with spectrum_utils' ProForma parser after prediction,
at ~3.5 ms per PSM, which is ~50 min for 825k beam rows. `proteoscan/casanovo_fastio.py`
(`python -m proteoscan.casanovo_fastio sequence ...`, used by default by
`run_denovo_queries`) replaces that with a regex. Predictions and scores are byte-identical, and
the mzTab `modifications` column becomes `null`.

Outputs (`/drive_a/drew/proteoscan-work/denovo_features/runs/`):
- `<run>.cand.denovo_all.parquet`: every beam candidate for the queries that have database
  candidates. Columns: query_id, scan (raw), precursor_idx, rank, sequence (ProForma), peptide
  (PEAKS notation), score, per_residue_scores, charge, exp_mz, ppm, precursor_ok, picked, shard.
- `DrewR_Adam_91526_DDA_OT_4.rest.denovo_all.parquet`: the remaining z <= 4 queries of OT_4.
- `DrewR_Adam_91526_DDA_OT_4.denovo_all.parquet` / `.denovo.parquet`: all OT_4 queries, all
  beams / picked call only.
- per-run shard MGFs (~10 GB, re-creatable), mzTabs and Casanovo logs under `runs/<run>/`.

OT_4 in total: 428,716 of the 429,200 z <= 4 queries (192,090 raw scans) have a call. The
53,626 queries with z 5-8 are outside Casanovo's range. Wall times include MGF writing (1-3 min per
run).

Candidate tables with the features, next to the originals:
`/drive_a/drew/proteoscan-work/prop/<run>/candidates_features_2pass_denovo.parquet`
(`scripts/add_denovo_features.py`; OT_4 1.6 M candidates: 71 s wall, 11 GB RAM). The query
coverage is 99.4-99.9 %.

## Validation against PEAKS (propionamide rerun export, 4 runs)

**What PEAKS 13.5 exports is not the de novo CAA score.** In `db.psms.csv`:
- `Tag Length` equals `Length` in 100 % of rows. The tag threshold is 0 %, as in
  `tag(>=0.0%)`.
- `CAA (%)` equals the percentage of residues whose `Positional Confidence` is >= 2 in 100 % of
  the 13,594 rows. 2 % is the report's "Confident Amino Acid Threshold" in
  `parameters.txt`.

So both columns describe the database peptide's own positional confidence (fragment evidence).
Neither is the de novo-tag CAA of the PEAKS DB paper, which feeds the score but is not
exported. `CAA (%)` itself barely tracks PEAKS's score (Spearman with -10LgP 0.05).

On the 9,920 PEAKS target PSMs that are among our candidates (same file, scan, modified
sequence and charge, I = L):

| our feature | PEAKS column | Spearman |
|---|---|---|
| tag_len | Tag Length (= Length) | -0.23 |
| caa_pct | CAA (%) | 0.36 |
| caa | CAA (%) | 0.24 |
| denovo_conf_mean | mean Positional Confidence | 0.22 |
| caa_pct | -10LgP | 0.29 |
| tag_len | -10LgP | 0.38 |
| max_caa_any_beam | -10LgP | 0.44 |

The negative tag_len correlation is peptide length: de novo agreement falls with length.

The features separate PEAKS-listed candidates from the rest (3.57 M candidates, all ranks):

| group | n | caa mean (median) | caa > 0 | caa_pct mean | tag_len mean | max_caa_any_beam mean | conf mean | denovo_exact |
|---|---|---|---|---|---|---|---|---|
| PEAKS-listed targets | 14,104 | 5.66 (6) | 84.9 % | 49.6 | 5.03 | 6.44 | 0.77 | 23.7 % |
| other targets | 1,770,176 | 0.50 (0) | 33.8 % | 5.3 | 0.46 | 0.70 | 0.28 | 0.16 % |
| decoys | 1,688,825 | 0.45 (0) | 33.0 % | 4.7 | 0.42 | 0.64 | 0.27 | 0.01 % |

Rank-1 candidates only: listed 6.23 / 87.5 % / 27.0 % exact, other targets 0.57 / 34.5 % /
0.46 %, decoys 0.47 / 33.2 % / 0.02 %. Unlisted targets sit only slightly above decoys, which fits
mostly-wrong candidates. The ~33 % of wrong candidates with caa > 0 are short chance matches
(1-2 residues).

The same features computed from **PEAKS's own de novo** (`<file>.raw.denovoAllCandidates.csv`,
highest ALC per scan and z, LC / 100 as confidence): 84.5 % of the listed PSMs have a PEAKS de
novo row (the export keeps ALC >= 50 only). On those rows:

| | PEAKS de novo | ours |
|---|---|---|
| exact call = the listed peptide (I = L) | 12.7 % | 25.9 % |
| mean caa_pct | 47.5 | 53.0 |

Spearman(PEAKS-de novo caa_pct, ours) is 0.61. Our calls equal PEAKS's de novo call for only
1.3 % of 316k queries (9.9 % where PEAKS ALC >= 80), the same disagreement as in docs/DENOVO.md.

## Score model (scripts/train_score_model.py, CV by scan parity)

The training script already lists `caa, caa_pct, tag_len, max_caa_any_beam, denovo_conf_mean,
denovo_exact` among its optional features. It also loads `ion_features.parquet` from the
candidate file's directory when present (neutral-loss / 2+ fragment features, added in
parallel). The official v2 therefore contains base + ion + de novo features. Two ablations
separate the contributions; their pkl/json are in `/drive_a/drew/proteoscan-work/denovo_features/models/`.
Absent feature columns are now filled with 0 when a model is applied (`predict` in the script,
`rescore.apply_imitation` with a warning), so a v2 model still scores runs without de novo.

The table sums the 4 runs. Recall and precision are against PEAKS's accepted set (its own 1 %
sequence-FDR threshold).

| model | features | CV Pearson (matched) | MAE | own 1 % FDR: n / recall / precision | at PEAKS threshold: n / recall / precision |
|---|---|---|---|---|---|
| score_prop_v1 | base | 0.609 | 8.59 | 13,678 / 0.780 / 0.594 | 7,508 / 0.595 / 0.825 |
| score_prop_v1b_ions | base + ions | 0.620 | 8.30 | 13,039 / 0.779 / 0.622 | 7,499 / 0.603 / 0.838 |
| denovo_noions (ablation) | base + our de novo | 0.626 | 8.42 | 13,540 / 0.783 / 0.603 | 7,568 / 0.603 / 0.830 |
| **score_prop_v2** | base + ions + our de novo | **0.631** | **8.22** | 13,094 / 0.780 / 0.620 | 7,508 / 0.605 / 0.839 |
| peaksdn_ions (ablation) | base + ions + de novo features from PEAKS's own de novo | 0.706 | 7.39 | 12,398 / 0.781 / 0.656 | 7,416 / 0.613 / 0.860 |

Our de novo features add +0.017 Pearson over base and +0.011 over base + ions. The fold-level
PSM sets barely move: the gains are in ranking PEAKS's scores, not in which PSMs pass. In
importance, `denovo_conf_mean` ranks with ms2_cos and ΔRT, followed by `max_caa_any_beam` and
`caa_pct`.

The same features computed from PEAKS's own de novo add +0.086. Most of the missing PEAKS
signal is therefore specific to PEAKS's de novo. It is not de novo accuracy: ours matches
the listed peptide twice as often. PEAKS DB scores agreement with *its own* tag, and whether a
spectrum has an ALC >= 50 PEAKS de novo at all (30-52 % of candidates) is itself informative.
That ablation cannot be used on new runs (it needs PEAKS). It sets the ceiling a PEAKS-like de
novo would give.

## Limitations / next steps

- Casanovo has no propionamide-Cys token. Such residues cannot match (the C-anchored alignment
  limits the damage), and they bias `denovo_exact` to 0 for C-containing propionamide peptides.
- Only ~10 % of picked calls fit the precursor within 15 ppm (PEAKS: 100 % of its de novo
  rows). Mass-constrained decoding would make the calls, and the CAA, more PEAKS-like.
- A PEAKS-like local confidence (agreement across many candidates) and an ALC >= 50 presence
  flag could narrow the gap to the PEAKS-de novo ceiling.
- The pipeline (`proteoscan.pipeline`) does not yet run de novo. Integrating
  `run_denovo_queries` + `add_denovo_features` before rescoring costs ~45 spectra/s per GPU.
