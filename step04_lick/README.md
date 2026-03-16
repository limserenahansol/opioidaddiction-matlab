# Step 04: Lick pipeline — downstream

Run **after Step 02**. Lick mega-pipeline and session-level lick patterns (PCA/k-means). Uses latest `run_*/ALL_mice_longitudinal.csv`.

## Scripts

- **`run_lick_mega_pipeline.m`** / **`run_lick_mega_pipeline_new.m`** — Full lick analyses. Writes to `run_*/figs/lick_MEGA/` and `per_session_features.csv` under `run_*/figs/lick_patterns_MASTER/`.
- **`run_lick_MEGA_followup2.m`** / **`run_lick_MEGA_followup2_new.m`** — Follow-up analyses.
- **`analyze_lick_patterns_MASTER.m`** — PCA + k-means on per-session lick features; reads `run_*/figs/lick_patterns_MASTER/per_session_features.csv`.
- **`rebuild_lick_patterns_with_ID.m`** — Rebuild lick-pattern plots with mouse IDs.

## Prerequisite

- **Step 02** (longitudinal CSV).
