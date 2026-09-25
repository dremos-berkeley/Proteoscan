# Status and known limitations (updated 2026-09-25 22:00 UTC)

Reference: DrewR_Adam_91526_DDA_OT_4 (Jurkat, Lumos, ion-trap HCD), PEAKS Studio DL boost on:
7,524 PSMs / 3,472 peptides / 990 protein groups.

| stage | state | evidence |
|---|---|---|
| export format round-trip | exact | 99 tests, all 4 reference dirs |
| protein inference / grouping / scores | exact except group numbering among ties | docs/PROTEIN_INFERENCE.md |
| featureless precursor correction | 92% agreement with PEAKS's corrected m/z | precursors_test |
| MS1 feature detection + chimera precursors | in progress | |
| candidate generation (Comet, instrument precursor) | 75% of PEAKS PSMs among top-5 candidates | comet_OT4 |
| rescoring (mokapot TDA) | 4,453 PSMs at 1% (recall of PEAKS 0.45) with misaligned MS2 feature; rerun pending | rescore_exp1 |
| rescoring (PEAKS-imitation, CV) | 5,471 PSMs at 1% (recall 0.55, precision 0.75), rerun pending | rescore_exp1 |
| de novo (Casanovo, high-res model) | 12% full-sequence accuracy on ion-trap data; unusable as is | denovo_casanovo |

De novo: Casanovo's released model is trained on Orbitrap MS2 and fails on ion-trap spectra
(2% agreement with PEAKS's de novo-only calls). Options: fine-tune Casanovo on the PEAKS-labelled
ion-trap PSMs once the HeLa PEAKS runs exist (tens of thousands of training spectra), or implement a
PEAKS-style dynamic-programming de novo for low-resolution data. Until then the denovoOnly export
is produced but not concordant.
