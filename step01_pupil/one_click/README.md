# Pupil one-click — one general U-Net + one call for every mouse x day (2026-11-12)

Korean version: [README_KR.md](README_KR.md)

The original Step 01a (`step1_2dec29_nove222222222.m` → `step3_dec291123342.m` → `alignment_pupil_salae_aug20.m`) needs, **for every
video**: editing paths in three scripts, clicking seeds on 2–3 frames, and training a U-Net from scratch (120 epochs) — 394 trainings for the
Dec-2025 cohort. This folder trains **one U-Net for all mice and days** from the ~200,000 pseudo-labels those trainings left behind, and runs
discovery → tracking → TTL alignment in one call. The three original scripts are unchanged; everything here is new.

## Run (MATLAB)

```matlab
% once: copy one_click_paths.example.json to one_click_paths.json (git-ignored) and set your folders
run_pupil_one_click('DryRun', true)   % lists sessions, TTL files, chosen videos, ROIs - writes nothing
run_pupil_one_click                   % everything; sessions already done are skipped
run_pupil_one_click('VideoRoot', V, 'SessionRoot', S, 'OutRoot', O)   % or pass the folders
run_pupil_one_click('Days', [12 13], 'Mice', {'0910_red'})          % a subset
run_pupil_one_click('ROIMode', 'previous')                          % never ask: reuse the mouse's previous ROI
```

The only manual step is the ROI, and only for sessions without `roi.mat`. They are asked **all at the start**: the previous day's ROI is
shown as a movable box — adjust and double-click it, or close the window to keep it. After that the run is unattended.

## Expected layout (as in Dec-2025)

```
VideoRoot\dayN\<cage>\<color>\*.avi  (+ roi.mat)     two cameras per mouse is fine
SessionRoot\dayN\<cage>\<color>\digital.csv           Saleae export: Ch6 camera, Ch7 lick, Ch4 injector
```
If `digital.csv` is not in `SessionRoot`, it is searched anywhere under `VideoRoot\dayN\` as `...\<cage>\<color>\digital.csv`
(the Dec-2025 TTL folders were called `ttl`, `tt`, `day4Jan2`, `Dec29day1`, ...).

## What it does per session

1. Counts camera TTL rising edges in `digital.csv`.
2. **Pupil video** = the `.avi` whose frame count is closest to the TTL count (right in 243/244 Dec-2025 sessions; the other had two
   videos of identical length, which is flagged). Frame counts use `Duration x FrameRate` (same as `NumFrames` within 1 frame, < 1 s
   instead of ~10 s per 3-GB file).
3. ROI from `roi.mat` (see above if missing).
4. `pupil_track_session.m` — general U-Net. Preprocessing as step 1 made its training images (ROI crop → gray → glint suppression →
   CLAHE), then **step 3's tracking rules unchanged** (30-frame warm-up, 30-px search radius, reject centre jumps > 18 px or diameter
   changes > 18 % and hold the last value, EMA alpha 0.18 / 0.12). Same output files and columns as step 3
   (`pupil_diameter_FINAL_U-Net_results2.csv`, `..._CLEAN.csv` with `SmoothedDiameter_px`, `Final_U-Net_TimeSeries.png`) plus
   `qc_frames.png` (12 frames with the mask outline) and `qc.json`.
   **Difference:** old step 3 skipped video frame 1 (pupil rows = frames − 1); this processes every frame, so rows = frames = camera TTLs
   and the alignment is `start_1to1`.
5. `pupil_align_ttl.m` — the original `alignment_pupil_salae_aug20.m` body **copied unchanged** into a function (only the inputs and the
   output folder are arguments); writes `concat_out_<stamp>\` exactly as before.
6. Report `OutRoot\one_click_report_<stamp>.csv`: one row per session, QC numbers, and a `flags` column (no TTL, video vs TTL off by
   > 2 frames, ambiguous camera, ROI reused, no pupil in > 5 % frames, > 10 % frames rejected, diameter > 35 % from the mouse's median,
   failed).

**Where outputs go:** Step 02 counts **every** `concat_out*` folder in a session folder as a separate session, so keep `OutRoot`
separate. `OutRoot = SessionRoot` is only for new sessions without a `concat_out`; the runner refuses to add a second one.

## The general model

`model\general_pupil_unet.mat` — not in git (`*.mat` ignored; trained on unpublished video). See [model/README.md](model/README.md).
128×128 input, U-Net depth 4, 32 filters, generalized Dice loss, 12 epochs. Training data: per-video runs whose own validation Dice was
≥ 0.90 with ≥ 100 pairs, up to 120 evenly spaced pairs per run; augmentation (rotation, scale, shift, brightness / contrast / gamma, blur,
noise). Rebuild:
1. `python -X utf8 inventory_existing_pupil.py` — sessions + runs (read-only)
2. `python -X utf8 make_training_list.py` — training list without 3 held-out mice (`--all` for all 14)
3. MATLAB `pupil_train_general(train_pairs.csv, out.mat, 'Epochs', 12, 'TestCsv', test_pairs.csv)` — ~40 min on a GPU

## Quick test with toy data

```
python -X utf8 make_toy_dataset.py     % 3 real sessions, first 900 frames + matching TTLs
```
```matlab
T = fullfile(<work_dir>, 'toy');
run_pupil_one_click('VideoRoot', fullfile(T,'videos'), 'SessionRoot', fullfile(T,'sessions'), 'OutRoot', fullfile(T,'out'), 'ROIMode', 'previous')
```
```
python -X utf8 test_toy_compare.py     % diameter vs the existing per-video U-Net result
```
The toy set contains two cameras per mouse, a day without `roi.mat`, and a TTL file in a `ttl\` sub-folder on purpose.

## Checks done (2026-11-12, toy only — not run on the full data)

- **Held-out mice** (0910_red, 0911_white, 6100_orange — 64 runs, 1,920 labels, never used in training): Dice median **0.952**
  (IQR 0.872–0.969); per-video models reached 0.962 on frames of their own video (optimistic). The 3-epoch model gave 0.943.
- **Toy one-click** (3 sessions x 900 frames, two of the held-out mice): 3/3 ok; re-run skipped 3/3; video pick, TTL search, ROI reuse
  and `start_1to1` alignment (900 rows) all as expected; 50–125 frames/s (≈ 3 min per full session). Lick onsets 34/34 and 24/24 TTL lick
  edges; 44/56 for 6100_orange — the original alignment samples the lick line at camera frames, so licks within one frame merge.
- **Versus the old per-video results:** the general model is lower by 6–12 % (medians 34.1 / 34.7 / 37.7 vs 39.5 / 38.4 / 40.5 px).
  Zoomed frames show why: when the pupil constricts, the darker pupil disc is visible inside the grey iris; the general model follows it
  (26–30 px) while the old models keep segmenting the whole dark eye including the iris (38–40 px). When the pupil is dilated both agree
  (40.0 vs 40.5 px). So the old pupil data may under-report constriction — seen in 3 short toy sessions only (frame checks made with the
  3-epoch model).
- **Known:** glint suppression can leave a notch in the mask edge, which lowers the diameter a little; the diameter definition
  (EquivDiameter of the mask) was kept as in step 3. A convex hull of the mask would remove the notch.
