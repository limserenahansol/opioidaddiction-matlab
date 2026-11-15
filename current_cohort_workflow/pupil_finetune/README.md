# Fine-tuning the pupil U-Net for a new rig or cohort

On the current rig (1920x1080 side view), the general U-Net (`general_pupil_unet.mat`, trained on the previous rig's 720x540 video) missed or mislabelled most pupils. This folder holds the procedure used on 2026-11-15. It labels frames of the new cohort, fine-tunes the U-Net on them, and checks the result on held-out mice or days.

Use the result with the pupil/TST package:
- `"pupil_detector": "maxp"` keeps the most probable U-Net blob, if its probability is >= `"maxp_threshold"` (0.9).
- `"pupil_model": "pupil/model/<new model>.mat"`

See `pupil_tst_patch/README.md` for the package edits.

## Steps

1. **Check the eye ROIs.**
   - Run `python pupil_roi_check.py --root <project> --day N` (in `python/daily`).
   - Every eye must be inside its box, with some margin. Mice sit differently every day.
2. **Sample frames.**
   - Run `python pupil_diag_extract3.py <project> <diag dir>`.
   - It takes 24 frames per session, spread over the session, with the window = ROI x2.
3. **Make the candidate sheets.**
   - Run `make_candidate_sheets4(<diag dir>, <sheet dir>, <general model .mat>, <package pupil dir>)` in MATLAB.
   - The sheets have 4 frames each. Every frame is shown plain and with numbered candidate contours: the U-Net blobs, plus dark regions from `pupil_hybrid_candidates.m`.
   - It also writes `candidates.csv` and `<session>_cands.mat` (the masks).
4. **Label.**
   - Two independent labellers per frame choose the contour that outlines the pupil, 0 (visible, not outlined) or -1 (not visible). A third labeller decides the disagreements.
   - `label_workflow_prompt.js` is the Claude Code workflow used for this, with its exact instructions. Labelling by hand with the same rules works too.
5. **Consensus and same-pupil sets.**
   - `fit_pupil_rule3.py <sheet dir> <fit dir> <diag dir>` needs `labels_consensus.csv` in `<fit dir>`.
   - It also compares candidate-score rules (logistic / boosting). On this cohort they were worse than fine-tuning.
6. **Fine-tune and check.**
   - Held-out mice: `finetune_unet_cv(..., 'Fold', 1..4)`.
   - Held-out day: `finetune_unet_split(..., 'Split', 'day1'|'day2')`.
   - Then compare the mask post-processing variants with `eval_ft_variants(..., 'folds'|'split')`.
   - `agreed_eval.py` scores only the frames where the two labellers agree. `sum_ft_eval.py` and `compare_methods.py` give totals.
7. **Final model.**
   - Run `finetune_unet_cv(..., 'Fold', 0)`, which trains on all labelled frames.
   - Copy `unet_ft_fold0.mat` to `pupil/model/` in the package.

## Result on the current cohort (2026-11-15)

Frames where both labellers agree; held-out 2 mice or 1 day:

| Detector | Pupil found | Wrong picks |
|---|---|---|
| General U-Net (package) | 44 % | 45 % |
| Fine-tuned U-Net + `maxp` 0.9 | 82-83 % | 22-24 % |

Errors cluster in sessions where the pupil is hard to see even by eye:
- inflamed or narrowed eyes;
- a dark eye without a visible pupil edge;
- a dark spot under the upper lid.

The package QC sheet and `qc.json` flag these sessions.
