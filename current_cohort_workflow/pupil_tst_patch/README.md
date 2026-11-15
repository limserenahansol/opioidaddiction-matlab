# Pupil + TST package: edits for the current rig

The pupil + TST package (`auto_pupil_tst_pipeline_20261114_iris_on`, run with `RUN_ALL_18DAYS.bat`) was written for the previous rig:
- DV `.avi` video at 720x540;
- cage folders named with 4 digits (`0910`);
- Saleae Channel 4 named `Channel 4`;
- colours `black / orange / red / white`.

The current rig differs in four ways, and without changes the package finds **0 sessions**:

| | Previous rig | Current rig |
|---|---|---|
| Video | `*.avi`, 720x540 DV | `*.mp4`, 1920x1080 H.264 (IC Capture), 60 fps header, one frame per camera TTL |
| Cage folder | `0910` | `9670f`, `9669m` (4 digits + sex letter) |
| Saleae reward column | `Channel 4` / `injector` | `reward` |
| Mouse colour | letters only | also `green1`, `green2` |

`current_cohort_edits.patch` holds all edits, applied to a copy of the package:
- The minimal edits were approved by the PI on 2026-11-14.
- The pupil detector was added on 2026-11-15.
- Every changed line is marked `CURRENT-COHORT EDIT` with its date. The original zip is not changed.
- **Every new option defaults to the old behaviour**, so the package still runs the previous rig exactly as before.

## 2026-11-14: names and formats (needed to find the sessions)

| File | Change |
|---|---|
| `pupil/run_pupil_one_click.m` | Video candidates are `*.avi` **and** `*.mp4` (new `list_videos` helper). Cage folder `^\d{4}[a-zA-Z]?$`. Session names like `day1 9670f_green1` parse with a second pattern, tried only when the original pattern fails. |
| `pupil/pupil_align_ttl.m` | `reward` / `Reward` added **last** to the injector column names. |
| `RUN_ALL_18DAYS.m` | The same second pattern for `redo` / `redraw_roi` names and for `mice` entries (`9670f_green1`). |
| `tst/scoring_io.py` | TST video names may have a sex letter after the cage and the colours `green1/green2` (`9670fgreen1tst` -> `9670_green1`). |

## 2026-11-15: pupil detection on the current rig

The general U-Net was trained on the previous rig's video. On the current 1080p side view it found the pupil in only 44 % of the visible-pupil frames, and 45 % of its picks were wrong. These numbers are on labelled frames where two independent labellers agree.

The fix:
1. The ROIs were checked and corrected (whole eye opening + margin).
2. 384 frames were labelled.
3. The U-Net was fine-tuned on them. See `../pupil_finetune/README.md`.

On held-out mice or days, the fine-tuned U-Net with the `maxp` blob choice found 82–83 % of pupils, and 22–24 % of its picks were wrong.

| File | Change (all options off by default) |
|---|---|
| `pupil/pupil_mask_frame.m` | Optional 7th output: the pupil probability map, resized to the crop. It is computed only when asked for. |
| `pupil/pupil_maxp_mask.m` (new) | Keeps the one blob with the highest pupil probability, if that probability is >= the threshold (0.9). Otherwise there is no pupil in that frame. |
| `pupil/pupil_track_session.m` | `Detector` `unet` (default) / `hybrid` / `maxp`. `ReacquireFrames` N: after N frames without an accepted pupil, the tracker forgets its last centre, diameter **and** warm-up anchor. `MaxPThr`. With `hybrid` / `maxp` it writes `hybrid_detections.csv`: the raw detection of every frame with its score, and the best guess even below the threshold. `qc.json` records the detector, threshold and model. |
| `pupil/pupil_qc_sheet.m` | Uses the same detector as the tracker (read from `qc.json`), so the replay check still holds. |
| `pupil/run_pupil_one_click.m` | Passes `Detector`, `ReacquireFrames` and `MaxPThr` through. |
| `RUN_ALL_18DAYS.m` | New config keys: `pupil_detector`, `reacquire_frames`, `maxp_threshold`, `pupil_model`. A relative `pupil_model` is resolved against the package folder. |
| `pupil/pupil_hybrid_mask.m` (new, superseded) | First rule-based attempt (U-Net blob with an edge, else a round dark region). Kept for comparison only. |
| `pupil/model/README.md` | Describes `current_rig_pupil_unet.mat`. The model file itself is not in git (`*.mat`). |

Current-rig settings in `run_all_config.json`:
```json
"pupil_detector": "maxp",
"maxp_threshold": 0.9,
"reacquire_frames": 15,
"pupil_model": "pupil/model/current_rig_pupil_unet.mat",
"expected_frame_size": "1920x1080"
```

## Daily use

`DAILY_summary.bat` (`tools/daily_summary.py`) runs the package for the day. It copies these keys from `run_all_config.json`, then makes `analysis/days/dayN/figs/01_pupil_roi_check.png` (`pupil_roi_check.py`).
- The package re-uses the previous day's ROI for a new day. If an eye is not inside its yellow box, set `"redraw_roi": ["dayN cage_colour"]` and run `RUN_ALL_18DAYS.bat`. It opens the ROI window and redoes that session.
- A narrowed or inflamed eye cannot be fixed with the box. With `maxp` it is often not "no pupil" but a confident wrong blob (the slit, a cloudy patch), so the "no pupil" fraction alone does not flag it.
- **Session grade.** `pupil_session_quality.py` writes `pupil_tst_final/pupil_session_quality.csv`: a manual grade from `pupil_tst_final/pupil_grades.csv` when there is one, else an automatic stability grade (`ok` or `check`; see `../pupil_finetune/README.md`). Grades: `ok` and `noisy` (the pupil, with dips) are used by the daily analysis; `guess` (often not the pupil) and `no_pupil` are not. For a `check` session, look at its `qc_sheet.png` and add a row to `pupil_grades.csv`.
- For `guess` sessions, `pupil_guess_trace.py` makes a separate best-guess trace from `hybrid_detections.csv`. The guess is labelled as such and never replaces the tracked pupil.

## Apply to a fresh copy of the package

```bash
cd <copy of auto_pupil_tst_pipeline_20261114_iris_on>
```
```bash
git apply --whitespace=nowarn <repo>/current_cohort_workflow/pupil_tst_patch/current_cohort_edits.patch
```

Then set `run_all_config.json` for the current rig:
- the keys above;
- `"video_root"`: the cohort's `video` folder;
- `"session_root"`: the `salae` folder;
- `"tst_video_root"`: the `tst_video` folder.

Then copy `current_rig_pupil_unet.mat` into `pupil/model/`. It is not in git; keep it with the cohort data.

## Known limits

- **Dark images.** The eye ROI uses about 10–45 of 255 grey levels. More IR light or a higher camera exposure / gain would help any detector.
- **Hard sessions.** Inflamed or narrowed eyes, a very dark eye without a visible pupil edge, and a dark spot under the upper lid remain unreliable. Use the QC sheets and `qc.json` `fracNoBlob`. The list of affected sessions is kept with the cohort data, not here.
- **TST.** The TST model expects whole 3–4 min recordings at 30 fps. The current camera writes a 60 fps header, so check the first TST videos (day 5).
- **Video vs camera pulses.** When a video has more or fewer frames than camera pulses, `pupil_align_ttl` (the original Step 01a logic, unchanged) matches them at the END. On the current rig the video normally starts with the first pulse: the recording stopped early, or test frames were added at the end. A motion-vs-lick cross-correlation confirmed this for both such sessions so far. `pupil_session_quality.py` compares the package's camera times with `sync_sessions.py` (`align_shift_s`, `align_check`) and lists every session where they differ. The package's output is not changed.
  - For a session the user approves (`alignment_fixes.csv`), `fix_concat_alignment.py` writes a corrected copy next to it, `concat_out_<stamp>_syncaligned`. It uses `../matlab/pupil_align_ttl_sync.m`, which is `pupil_align_ttl.m` with only the frame-to-pulse match taken from `sync_sessions.py`.
  - Step 02 then reads the corrected copy (the newest folder).
