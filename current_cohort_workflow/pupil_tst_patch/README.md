# Pupil + TST package: edits for the current rig

The pupil + TST package (`auto_pupil_tst_pipeline_20261114_iris_on`, run with `RUN_ALL_18DAYS.bat`) was written for the previous rig. That rig had DV `.avi` video at 720x540, cage folders named with 4 digits (`0910`), Saleae Channel 4 named `Channel 4`, and colours `black / orange / red / white`.

The current rig differs in four ways, and without changes the package finds **0 sessions**:

| | Previous rig | Current rig |
|---|---|---|
| Video | `*.avi`, 720x540 DV | `*.mp4`, 1920x1080 H.264 (IC Capture), 60 fps header, one frame per camera TTL |
| Cage folder | `0910` | `9670f`, `9669m` (4 digits + sex letter) |
| Saleae reward column | `Channel 4` / `injector` | `reward` |
| Mouse colour | letters only | also `green1`, `green2` |

`current_cohort_edits.patch` holds the minimal edits, applied to a copy of the package and approved by the PI on 2026-11-14. Every changed line is marked `CURRENT-COHORT EDIT 2026-11-14`. The original zip is not changed.

| File | Change |
|---|---|
| `pupil/run_pupil_one_click.m` | Video candidates are `*.avi` **and** `*.mp4` (new `list_videos` helper). Cage folder `^\d{4}[a-zA-Z]?$`. Session names like `day1 9670f_green1` parse with a second pattern, tried only when the original pattern fails. |
| `pupil/pupil_align_ttl.m` | `reward` / `Reward` added **last** to the injector column names. |
| `RUN_ALL_18DAYS.m` | The same second pattern for `redo` / `redraw_roi` names and for `mice` entries (`9670f_green1`). |
| `tst/scoring_io.py` | TST video names may have a sex letter after the cage and the colours `green1/green2` (`9670fgreen1tst` -> `9670_green1`). |

Checked on 2026-11-14:
- Every previous-cohort name (`day12 0910_red`, `0910red`, `day12/0910/red`, `6100_Orange`, `0910blackhot`, ...) parses exactly as before.
- The dry run finds all 16 current-cohort sessions of day 1–2.
- 14 of the 16 have a video frame count equal to the camera TTL count. The other two are 30- and 26-minute sessions recorded without the time limit: +150 and −212 frames.

## Apply to a fresh copy of the package

```bash
cd <copy of auto_pupil_tst_pipeline_20261114_iris_on>
```
```bash
git apply --whitespace=nowarn <repo>/current_cohort_workflow/pupil_tst_patch/current_cohort_edits.patch
```

Then set `run_all_config.json` for the current rig:
- `"expected_frame_size": "1920x1080"`
- `"video_root"`: the cohort's `video` folder
- `"session_root"`: the `salae` folder
- `"tst_video_root"`: the `tst_video` folder

## Known limitation (to fix last, as agreed)

The general U-Net was trained on the previous rig's video. A 1500-frame test of day 1–2 on 2026-11-14 tracked one session well: `9669m_green2` day 1, with no pupil found in 0.9 % of frames. Most other sessions found no pupil in 48–99 % of frames; the eye is small and dark in the 1080p side view. The package's QC sheet and `qc.json` (`fracNoBlob`) flag these sessions. The daily analysis drops sessions with `fracNoBlob > 0.5` from pupil figures.
