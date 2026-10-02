# Trained models

- `imitation_jurkat_v1.pkl`: LightGBM classifier predicting whether PEAKS Studio (DL boost on)
  reports a candidate PSM, trained on DrewR_Adam_91526_DDA_OT_4 (Jurkat, ion-trap HCD) with
  the two-pass feature set. Cross-validated on that run: at p >= 0.5 it reproduces 56% of PEAKS's
  PSMs with 77% precision; at p >= 0.3, 63% / 71%. Use with
  `proteoscan run --mode imitate --model models/imitation_jurkat_v1.pkl`. It has seen only one
  run; retrain when PEAKS exports for the HeLa runs exist.
- `casanovo_iontrap_hcd_ft5runs_step15000.ckpt` (not tracked in git, 400+ MB; copy from
  /drive_a/drew/proteoscan-work/denovo_finetune/models/): Casanovo 5.2 fine-tuned on 197k
  ion-trap HCD PSMs from the HeLa runs. Held-out HeLa run: 42% full-sequence accuracy (I=L) vs 17%
  for the released model; Jurkat (turbo-scan acquisition): 20% vs 18% with 5 beams and
  precursor-consistent candidate selection. See docs/DENOVO.md.
