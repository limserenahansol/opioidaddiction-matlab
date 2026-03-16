# Step 02: Lick pipeline (run second)

Run after **Step 01**. Full lick analyses: core, microstructure, rhythm, PSTH, learning, prediction, unsupervised, sequence.

## Scripts

- **`run_lick_mega_pipeline.m`** / **`run_lick_mega_pipeline_new.m`** — Main lick mega-pipeline. Writes to `run_*/figs/lick_MEGA/` and can produce `per_session_features.csv` under `run_*/figs/lick_patterns_MASTER/`.
- **`run_lick_MEGA_followup2.m`** / **`run_lick_MEGA_followup2_new.m`** — Follow-up analyses.
- **`analyze_lick_patterns_MASTER.m`** — PCA + k-means on per-session lick features; reads `run_*/figs/lick_patterns_MASTER/per_session_features.csv`; writes to `.../plots_matlab/`.
- **`rebuild_lick_patterns_with_ID.m`** — Rebuilds lick-pattern plots with mouse IDs.

## Order

1. Run **Step 01** first.
2. Run **`run_lick_mega_pipeline`** (or `_new`).
3. Optionally **`analyze_lick_patterns_MASTER`** or **`rebuild_lick_patterns_with_ID`**.

## Next

Then run **Step 03** (motivation), then **Step 04** (build longitudinal).
