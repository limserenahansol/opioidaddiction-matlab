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

**Which blob, when several pass the threshold.** `eval_ft_variants2(..., 'folds')` compares the choice among the blobs with max p >= 0.9 on the held-out folds (all 384 labelled frames, so the numbers are lower than on agreed frames). 133 of 384 frames had more than one such blob.

| Choice among confident blobs | Pupil found | Wrong picks |
|---|---|---|
| highest max p (`maxp`, used) | 75.7 % | 30.9 % |
| closest to the brightest spot (glint) | 73.6 % | 32.9 % |
| largest sum of p | 67.5 % | 38.4 % |
| most circular | 65.7 % | 40.1 % |
| perfect choice (upper bound) | 80.4 % | 26.7 % |

So `maxp` stays. A better tie-break can gain at most about 5 points; the remaining errors are frames where no confident blob is the pupil.

## Per-session grade after the full run

All 16 day-1/2 sessions were graded by eye from their QC sheets (`ok`, `noisy`, `guess`). Two automatic numbers from `hybrid_detections.csv` separated the 7 `ok` sessions from the other 9 without error (`python/daily/pupil_session_quality.py`):
- `far`: how often the detection centre is more than half a pupil diameter from its 1-s running median;
- `dev30`: how often the detected diameter is more than 30 % off its 1-s running median.

`ok` sessions had both <= 0.028; each of the others had `far` >= 0.032 or `dev30` >= 0.049. The automatic grade is `ok` when both are <= 0.04, else `check`. A `check` session still needs a look at its QC sheet, because `noisy` (the pupil, with dips) and `guess` (often not the pupil) are not told apart by these numbers.

## Fully dilated pupils (pupil fills the eye opening)

When the whole eye is dark, the pupil is fully dilated. The pupil is then the whole dark area inside the eyelids, not only its darkest part (rule from the user). On this rig that is about 100-135 px, against about 20-30 px for a normal pupil. Morphine dilates mouse pupils, so drug days can look like this.

The package's tracker cannot measure such a pupil:
- the U-Net, fine-tuned on frames labelled with a small dark part, finds only that part;
- the package's step 3 keeps only the part of the mask within 30 px of the last centre, so no pupil can be wider than 60 px.

**Retraining the U-Net did not help (not used).** The 72 labelled frames of the 3 fully dilated sessions (all from one cage, one day) were relabelled with the whole dark area: `make_fullpupil_labels.m` -> `add_fullpupil_labels.m`; 1 frame was dropped because its mask leaked into fur. The U-Net was then retrained with `finetune_unet_cv` on the same folds, and the folds were scored with the tracker's detector (`eval_fullpupil.m`: maxp >= 0.9, IrisCorrection on and off):
- On held-out mice, the chosen blob had 0.20-0.42 of the full diameter (median 0.29; 0.18 before), and none overlapped well.
- Even on its own training frames, the final model covered only 28-44 % of the full area.
- Normal pupils got slightly worse: correct 80 % vs 84 %, wrong picks 18 % vs 12 % on agreed frames.

Three sessions of one cage are too few to learn this.

**Direct measurement (used).** The rule is measured directly in sessions the user approves (`fullpupil_sessions.csv`):
- `../matlab/pupil_fullpupil_session.m` applies `fullpupil_mask.m` (the whole dark area inside the lids) to every frame.
- It uses the package's step-3 rules (jump rejection, EMA, re-lock, sgolay smoothing), without the 30-px search radius.
- A region that touches the ROI edge on more than 40 px and is more than 5 % wider than the last accepted value is treated as a leak into fur. An eye that is only cut by the ROI edge is measured as it is.
- No values are made up: gaps longer than 1 s (30 frames) and the session ends stay empty instead of being interpolated.
- A trace is not published to Step 02 when more than half of its frames have no value, or when its median is under 10 px.
- `../python/daily/fullpupil_sessions.py` runs it and writes a Step 02 copy of the newest table in which only the pupil columns are replaced.

**Screening (an option, not run daily).** `fullpupil_sessions.py --screen-only` looks for such sessions. Fully dilated pupils are rare, so the measurement is an option for listed sessions, not a rule for all sessions. The screen's test: when the ring just outside the tracked pupil is as dark as the pupil, there is no lighter iris around it.
- On days 1-2 the screen flagged the 3 fully dilated sessions and 3 normal ones; the numbers overlap.
- So a flag only asks for a look; it never changes the measurement.
- `dilation_features.m` gives the same numbers on the labelled frames.
