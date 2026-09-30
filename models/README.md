# Trained models

- `imitation_jurkat_v1.pkl`: LightGBM classifier predicting whether PEAKS Studio (DL boost on)
  reports a candidate PSM, trained on DrewR_Adam_91526_DDA_OT_4 (Jurkat, ion-trap HCD) with
  the two-pass feature set. Cross-validated on that run: at p >= 0.5 it reproduces 56% of PEAKS's
  PSMs with 77% precision; at p >= 0.3, 63% / 71%. Use with
  `peaksparity run --mode imitate --model models/imitation_jurkat_v1.pkl`. It has seen only one
  run; retrain when PEAKS exports for the HeLa runs exist.
