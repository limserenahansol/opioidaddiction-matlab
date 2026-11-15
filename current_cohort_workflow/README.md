# current_cohort_workflow

Daily handling, pupil/TST, Step 02 and cohort-level figures for a **new cohort recorded on the current rig**. The output is in exactly the format of the previous cohort, so the two can be compared and combined.

**Nothing of the previous pipeline is changed or removed:**
- `step01`–`step07` (and `step08`–`step10` where present) are used as they are.
- Step 02 is run through a wrapper that only sets its `BASE` folder.

Run order, logic diagram, definitions and the protocol differences between cohorts are in **[PIPELINE.md](PIPELINE.md)**. The analysis methods of the previous cohort, and how each is reproduced, are in **[docs/ANALYSIS_METHODS.md](docs/ANALYSIS_METHODS.md)**.

## Contents

| Path | What it does |
|---|---|
| `python/daily/daily_summary.py` | End of every day: find the day's task logs, Saleae captures and videos, pair them by start time, arrange `python|salae|video\dayN\<cage>\<mouse>\`, then sync and run the daily analysis. By default it also runs the pupil/TST package for that day (`--no-pupil` skips it), then makes the eye-box check sheet and updates `step02_base`. `--longitudinal` also rebuilds `ALL_mice_longitudinal.csv` with the original Step 02. |
| `python/daily/sync_sessions.py` | Task log + Saleae + video on one clock; per-session QC (rewards, licks, camera pulses vs video frames). |
| `python/daily/quick_analysis.py` | Day dashboards, rasters, cumulative and per-minute licking, ILI, peri-reward licking, PR breakpoint; across-day master tables. Uses the package's pupil output (`pupil_tst_final`) and skips sessions whose `qc.json` reports no pupil in more than 50 % of frames. |
| `python/arrange_assay_videos.py` | Assay days: TST / HP videos to `tst_video|hp_video\dayN\<cage>\<mouse>\`. |
| `python/prepare_step02_base.py` | After day 18: builds the folder Step 02 reads (newest `concat_out_*` of each mouse-day + its task log + assay files). |
| `matlab/run_step02_on_base.m` | Runs the **original** Step 02 script on any `BASE`; only the `BASE = '...'` line is replaced in a temporary copy. Verified: identical to the previous `run_009` for 14/14 day-1 sessions. |
| `python/plot_session_sync.py` | One-session check: pupil, Lick_TTL and Injector_TTL on the camera clock, reward-locked pupil and licking, and the numbers (frames vs camera TTLs, frame-sampled vs raw Saleae vs task-log licks, rewards). |
| `python/cohort_plots.py` | Cohort figures in the style of the previous analysis (PR, licking, bouts, trajectories, pupil, assays, multivariate), from `ALL_mice_longitudinal.csv`. |
| `python/harmonize_longitudinal.py` | Combine cohorts: union of columns, `cohort` column, duplicate sessions removed and listed, per-session table. |
| `pupil_tst_patch/` | Edits that let the pupil/TST package read the current rig (mp4, `9670f` cage folders, `reward` column, `green1` names), plus the current-rig pupil detector (`maxp` + fine-tuned U-Net, re-lock, per-frame detections). All new options default to the old behaviour. |
| `pupil_finetune/` | How the U-Net was fine-tuned for the current rig: sample frames, candidate sheets, two independent labels per frame, leave-mice-out / leave-day-out checks. Repeat for a new cohort or rig. |
| `python/daily/pupil_roi_check.py` | Daily eye-box check sheet (3 frames per mouse, ROI drawn, "no pupil" fraction). |
| `python/daily/pupil_session_quality.py` | Per-session pupil grade (`pupil_session_quality.csv`): the manual grade of `pupil_grades.csv`, else an automatic stability grade. `ok` / `noisy` sessions are used by the daily analysis, `guess` / `no_pupil` are not. |
| `python/daily/fix_concat_alignment.py` + `matlab/pupil_align_ttl_sync.m` | For sessions the user approved (`alignment_fixes.csv`): a new `concat_out_<stamp>_syncaligned` in which video frame k is matched to the camera pulse `sync_sessions.py` gives it. The package matches at the end when the counts differ. The package's output and code are not changed. Each copy is checked row by row and logged. |
| `python/daily/fullpupil_sessions.py` + `matlab/pupil_fullpupil_session.m` | Fully dilated pupils, i.e. the whole eye is dark: the pupil is then the whole dark area inside the lids. This is an option for rare sessions: only sessions listed in `fullpupil_sessions.csv` get a full-pupil trace (the package's filtering and smoothing) and a Step 02 copy with only the pupil columns replaced. Mice listed in `fullpupil_mice.csv` are handled automatically: each of their sessions is measured this way only when the black pupil region is larger than 60 px, the package's limit. An optional screen (`--screen-only`) looks for such sessions; a flag only asks for a look. |
| `python/daily/pupil_guess_trace.py` | Best-guess pupil trace for hard sessions (narrowed or inflamed eye) from `hybrid_detections.csv`. It is labelled as a guess and never replaces the tracked pupil. |

## Requirements

- **Python 3.11+** with `python/daily/requirements.txt` (numpy, pandas, scipy, matplotlib, imageio-ffmpeg).
- **The pupil/TST package** needs MATLAB (Deep Learning, Image Processing, Computer Vision toolboxes) and Python 3.13 with scikit-learn 1.8. See the package README.
- **Steps 03–10** of the previous pipeline need MATLAB with the Statistics and Machine Learning Toolbox.

## Cohort folder layout (copy this one folder to move the analysis to another PC)

```
<cohort>\
  python\dayN\<cage>\<mouse>\*.jsonl        task logs (copied)
  salae\dayN\<cage>\<mouse>\*.sal, digital.csv
  video\dayN\<cage>\<mouse>\*.mp4 (+ roi.mat written by the pupil package)
  tst_video\dayN\..., hp_video\dayN\...    assay videos
  assays\dayN\<cage>\<mouse>\               HP manual scoring .mat, tail-immersion latency file
  sync\, analysis\                         daily outputs
  pupil_tst_final\                         pupil/TST package output (+ QC_index.html)
  step02_base\ -> longitudinal_outputs\run_###\ALL_mice_longitudinal.csv
  tools\                                   python/daily/*, the patched pupil/TST package, venvs
```

No absolute paths are stored in this repository. Paths are given on the command line or in the package's `run_all_config.json`, which stays local.
