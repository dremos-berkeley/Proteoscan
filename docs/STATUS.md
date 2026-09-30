# Status and known limitations (updated 2026-09-30)

Reference: DrewR_Adam_91526_DDA_OT_4 (Jurkat, Lumos, ion-trap HCD), PEAKS Studio, DL boost on:
7,524 PSM rows / 3,472 peptides / 990 protein groups. Concordance = exact modified sequence,
scan and charge (PSM); exact modified sequence (peptide); shared accession (protein group).

## Exactly reproduced from PEAKS's own tables
| item | state |
|---|---|
| CSV formats (all 8 tables, Java number formatting) | byte-exact round trip |
| protein inference, grouping, scores, coverage, areas | exact; group numbering among ties differs |
| PTM profile | byte-identical |
| modified-peptide table | all cells except 8 RT-shift values (0.0001 off) |
| precursor list (given PEAKS's features) | 98.8% of the precursors behind PEAKS PSMs |
| I/L variant rows and protein lists per PSM | reproduced |

## End-to-end from raw data (Comet + rescoring), Jurkat run
| configuration | PSMs reported | PSM recall | PSM precision | peptide recall | protein-group recall |
|---|---|---|---|---|---|
| Comet e-value only | 1,198 | 0.13 | 0.84 | | |
| + rescoring (mokapot), corrected precursors | 11,126 | 0.74 | 0.50 | 0.82 | 0.95 |
| + PEAKS-imitation model (cross-validated) | 10,722 | 0.75 | 0.52 | 0.83 | 0.95 |
| + two-pass search with protein support, imitation | 14,609 | 0.81 | 0.41 | 0.88 | 0.95 |
| same, cut at PEAKS's PSM count | 7,595 | 0.67 | 0.66 | 0.78 | 0.88 |

On scans identified by both, the peptide is identical in 97% of cases.

Entrapment validation (human + shuffled entrapment database of equal size): measured false
discovery proportion 0.5% at the PEAKS-style filter and 0.9% at nominal 1% PSM FDR (single
pass). The pipeline therefore reports more correct identifications than PEAKS at the same
nominal error rate; "precision" against PEAKS is low because PEAKS reports fewer PSMs, not
because the extra PSMs are wrong. Two-pass entrapment validation: pending.

## Findings that shaped the design
- Feature detection is not the bottleneck: our detector (70% recall of PEAKS features) and
  PEAKS's own feature list give the same candidate coverage (78%).
- The bottleneck is candidate ranking on weak/chimeric spectra. PEAKS rescues these with its
  protein feature; the two-pass search imitates that (coverage 83%).
- Restricting to zero isotope error loses recall (0.67 vs 0.75): keep Comet's -1..3 window.
- Sage needs ~300-500 GB for semi-tryptic + 3 variable mods; Comet needs 3 GB. Comet is the engine.
- MS2PIP 4.2 `correlate` misattaches vectors across PSMs sharing a spectrum; features are computed
  from `predict_batch` plus our own peak matching.

## Open items
- PEAKS unfiltered export with decoys (running on the Windows machine): gives PEAKS's score for
  every PSM, needed to model its score and threshold rather than only its accepted set.
- De novo for ion-trap spectra: Casanovo's released model reaches 12% full-sequence accuracy on
  this data; needs fine-tuning on ion-trap PSMs. The denovoOnly export is not yet concordant.
- features.csv export and Area columns from our own feature detector.
- HeLa benchmark (PXD004940) runs through the pipeline and through PEAKS.
- macOS packaging and the browser results viewer (in progress).
