# current_cohort_workflow

Daily handling, pupil/TST, Step 02 and cohort-level figures for a **new cohort recorded on the current rig**. The output is in exactly the format of the previous cohort, so the two can be compared and combined.

**Nothing of the previous pipeline is changed or removed:**
- `step01`–`step07` (and `step08`–`step10` where present) are used as they are.
- Step 02 is run through a wrapper that only sets its `BASE` folder.

Run order, logic diagram, definitions and the protocol differences between cohorts are in **[PIPELINE.md](PIPELINE.md)**. The analysis methods of the previous cohort, and how each is reproduced, are in **[docs/ANALYSIS_METHODS.md](docs/ANALYSIS_METHODS.md)**.

## Contents

| Path | What it does |
|---|---|
| `python/daily/daily_summary.py` | End of every day: find the day's task logs, Saleae captures and videos, pair them by start time, arrange `python|salae|video\dayN\<cage>\<mouse>\`, then sync and run the daily analysis. By default it also runs the pupil package for that day (`--no-pupil` skips it; TST is scored by `assays_day.py`, step 5b), the whole-eye step, the session and PR checks, the assay step (`--no-assays` skips it), updates `step02_base`, runs the original Step 02 + the previous cohort's figures (`--no-longitudinal` skips them), the QC summary, the data index and the routine PowerPoint (`--no-ppt`). |
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
| `python/daily/session_check.py` + `SESSION_CHECK.bat` | Read-only check right after each session, against the day's expected values in `protocol_expected.csv`. It prints PASS / CHECK for:<br/>- the file name, including the same mouse used twice on one date;<br/>- the protocol in the task log, including the measured lockout;<br/>- the Saleae valve open time and reward count;<br/>- Python-only licks that the task counted;<br/>- camera frames vs pulses, and whether the video was stopped after the task ended;<br/>- the Saleae settings.<br/>On a passive (replay) log it also checks the replayed log: same date, another mouse, finished first, with the same reward count and times.<br/>On a `test_...` run it checks only the valve time and the Saleae settings. Reports go to `analysis\session_checks\`. |
| `python/daily/pr_daily_checks.py` | Evening table, one row per mouse, report only (`analysis\days\dayN\pr_checks_dayN.*`, `figs\12_pr_checks.png`):<br/>- valve open time;<br/>- measured vs logged lockout;<br/>- Python-only licks that counted toward a ratio;<br/>- session length, last reward and final gap;<br/>- PR for the whole session and in the first 910 s, next to the previous cohort's days 3–5 (reference only). |
| `python/daily/qc_summary.py` | The daily report: one QC picture + CSV of all days so far (rewards Python = Saleae, PR, valve time, session length, licks, frames vs pulses, pupil grade / whole-eye off-eye %). Every CHECK item is explained (yellow) or red; the user's decisions on red flags (`analysis\days\qc_decisions.csv`) turn them yellow. |
| `python/daily/routine_ppt.py`, `make_summary_ppt.py` | The routine PowerPoint (QC per session, problems and fixes, folders, basic plots, next steps) and the multi-day summary deck. |
| `python/daily/data_index.py` | Step 10 of the routine: `analysis\data_index.csv` + `analysis\days\dayN\data_index_dayN.csv` = every saved data file (relative path, rows, size, which one is final), so the analysis after day 18 reads them and never re-runs the daily steps. Loader: `DI.paths(kind, day=, cage=, mouse=)`. |
| `python/daily/assays_day.py` | Step 5b of the routine (TST + hot plate on days 5, 10, 13, 16, 18): copies the day's assay videos to `tst_video` / `hp_video`; scores TST automatically on the **first 3 min** (`tst_mouse_lock.py`) and writes `assays\dayN\<cage>\<mouse>\ManualScoringResults_<video>_AUTO.mat` = the file name the **original Step 02** reads, so the `TST_*` columns are made as for the previous cohort (a hand score of the same mouse-day wins); makes first-30-s copies of the HOT videos (`hp_video\dayN_first30s\`) for the lab's hand scorer (`step01_manual_assays/manual_scoring_video_multibehavior_batch2.m`) and files its `.mat`; writes `analysis\days\dayN\assays_dayN.csv` (TST, HOT, tail immersion). |
| `python/daily/tst_mouse_lock.py` | The package's TST scorer, **not changed**, run on the first 180 s with a per-frame check that the tracked blob is the mouse (area, jumps, position in the upper middle). Only when the untouched run loses the mouse (> 1 % of frames, or most of any second; e.g. the camera was moved during a recording) are the moved frames aligned to the main camera view (ORB + RANSAC) and scored again; that result is used only if its check passes. |
| `python/daily/tst_examples.py`, `tst_clips.py` | Example sheets (3 longest moving / non-moving stretches per video, mouse outline now and 0.5 s before) and short annotated clips of the automatic TST calls, for the user to confirm. |
| `python/daily/tail_immersion_entry.py` | The tail-immersion latencies the user reports each day -> `assays\dayN\<cage>\<mouse>\<cage><mouse>_immersion_dayN.csv` (value in B1, read by the original Step 02) + `assays\tail_immersion_all.csv`. A different value for the same mouse-day is never overwritten silently. |
| `python/daily/wholeeye_*.py`, `pupil_finetune/eye_*.m` | Mice whose pupil cannot be seen (eye damage): a one-time U-Net of the whole eye between the lids from the user's outlines (drawing pages: `wholeeye_teach_export.py` + `pupil_teach_server.py`), tracking (`eye_track_all.m`), the CLEAN trace (`wholeeye_clean.py`) and a Step 02 copy with only the pupil columns replaced (`wholeeye_sessions.py`); `wholeeye_compare_models.py` compares two eye models on the same sessions. Opt-in per mouse (`wholeeye_mice.csv`). |
| `python/daily/besteffort_sessions.py`, `behavioronly_sessions.py` | Options for listed sessions only: a best-effort pupil copy for hard eyes; a behaviour-only Step 02 table for a session without video. |
| `python/daily/claude_jobs_guard.ps1` | Keeps long analysis jobs (MATLAB eye / pupil training and tracking) out of the way of the experiment: freezes them while the task runs, never kills them; resume after the day is over. |

## Requirements

- **Python 3.11+** with `python/daily/requirements.txt` (numpy, pandas, scipy, matplotlib, imageio-ffmpeg).
- **The pupil/TST package** needs MATLAB (Deep Learning, Image Processing, Computer Vision toolboxes) and Python 3.13 with scikit-learn 1.8. See the package README. `assays_day.py`, `tst_mouse_lock.py`, `tst_examples.py` and `tst_clips.py` run with that same Python (they import the package's `tst` folder from `<cohort>\tools\`).
- **Steps 03–10** of the previous pipeline need MATLAB with the Statistics and Machine Learning Toolbox.

## Cohort folder layout (copy this one folder to move the analysis to another PC)

```
<cohort>\
  python\dayN\<cage>\<mouse>\*.jsonl        task logs (copied)
  salae\dayN\<cage>\<mouse>\*.sal, digital.csv
  video\dayN\<cage>\<mouse>\*.mp4 (+ roi.mat written by the pupil package)
  tst_video\dayN\..., hp_video\dayN\...    assay videos (copies; TST first 3 min, HOT first 30 s are scored)
  hp_video\dayN_first30s\                  HOT scoring copies (first 30 s) + output\ManualScoringResults_*.mat
  assays\dayN\<cage>\<mouse>\               files the original Step 02 reads: tail-immersion latency (B1),
                                           ManualScoringResults_*_TST_*_AUTO.mat (or a hand TST score), HOT hand score
  sync\, analysis\                         daily outputs (analysis\tst\dayN\: automatic TST runs and checks;
                                           analysis\data_index.csv: every saved data file)
  pupil_tst_final\, pupil_tst_teach\       pupil/TST package output (pupil_active_root.txt names the one in use)
  step02_base\ -> longitudinal_outputs\run_###\ALL_mice_longitudinal.csv
  tools\                                   python/daily/*, the patched pupil/TST package, venvs, local_paths.json
```

No absolute paths are stored in this repository. Paths are given on the command line, in the package's `run_all_config.json`, or in `<cohort>\tools\local_paths.json` (this PC's raw-data folder, old project copies never searched, assay-video folders), which all stay local.
