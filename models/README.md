# Trained models

The pipeline's default scoring is mokapot's semi-supervised target-decoy rescoring, learned from each run's own decoys; it uses none of the models below. The PEAKS-trained models are a legacy option (`--mode imitate`) from the project's original goal of reproducing PEAKS's PSM lists; they are not used in any accuracy result.

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

## Validity after the 2026-10 corrections (docs/STATUS.md "Corrections")
- `score_prop_v1`, `score_prop_v1b_ions`, `score_prop_v2`: **invalidated**. Trained on candidates of
  the propionamide search, which carried a hidden static +57 on every Cys (P1), so their Cys
  candidates were C+114 / C+128 forms. Retrain from corrected candidates and the decoy-inclusive
  PEAKS export with `scripts/train_score_model.py`.
- `score_cam_unfiltered_v1`, `imitation_jurkat_v1`: trained on fixed-CAM searches (unaffected by P1)
  and still valid, but their feature lists contain `has_cys`, which is no longer a rescoring feature.
  They keep working (`apply_imitation` reindexes missing columns, and `has_cys` is still computed by
  `sequence_features`); retrain without it when next updated.
