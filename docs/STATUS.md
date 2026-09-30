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

## HeLa benchmark (PXD004940, 90-min gradients, 1% FDR, pipeline defaults; 2026-09-30)
| run | MS2 | PSMs | peptides | protein groups | wall min |
|---|---|---|---|---|---|
| 03_01 500 ng | ion trap | 47,030 | 38,138 | 5,172 | 35 (single pass) |
| 03_02 500 ng | ion trap | 46,728 | 37,764 | 5,239 | 35 (single pass) |
| 03_03 500 ng | ion trap | 49,123 | 40,092 | 4,979 | 58 |
| 05_01 1000 ng | ion trap | 55,767 | 44,565 | 5,149 | 63 |
| 05_02 1000 ng | ion trap | 54,572 | 42,966 | 5,166 | 57 |
| 05_03 1000 ng | ion trap | 55,563 | 44,108 | 5,184 | 58 |
| 03_01 1 ug | Orbitrap 0.02 Da | 30,958 | 27,543 | 3,992 | 32 |
| 03_02 1 ug | Orbitrap 0.02 Da | 32,111 | 28,694 | 4,029 | 26 |
| 03_03 1 ug | Orbitrap 0.02 Da | 32,145 | 28,633 | 4,024 | 38 |
Exports under /drive_a/drew/peaks-parity-work/hela/<run>/export. PEAKS results for these runs are
pending on the Windows machine; they will be the training set for the score-imitation model.

## De novo (2026-09-30)
Casanovo fine-tuned on the HeLa ion-trap PSMs (docs/DENOVO.md): held-out HeLa accuracy 42%
(released model 17%); Jurkat 20% vs 18% (different scan mode: turbo vs rapid). PEAKS's own de novo
agrees with database identifications on only 4% of its de novo-only spectra on this data, so
sequence-level concordance of the denovoOnly export is not a meaningful target on ion-trap data;
the export is produced in the PEAKS layout with our model's calls and confidences.
