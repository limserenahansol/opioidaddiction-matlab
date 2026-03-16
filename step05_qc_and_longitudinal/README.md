# Step 05: QC and longitudinal stats/plots — downstream

Run **after Step 02**. Quality control, requested analyses, and longitudinal plots/stats from the latest `run_*/ALL_mice_longitudinal.csv`.

## Scripts

- **`make_longitudinal_QC_and_requested_analyses_cursor.m`** / **`make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v5.m`** — QC + requested analyses. Writes to `run_*/QC_AND_REQUESTED_ANALYSES_<timestamp>/`.
- **`make_longitudinal_plotsall_statistic.m`** / **`make_longitudinal_plotsall_statistic_NEWCOHORT_2026.m`** — Longitudinal plots by week/group + stats; `run_*/figs/`, `figs/weeks`, `figs/stats`.
- **`make_stats_only_from_csv.m`** / **`make_stats_only_from_csv_dec.m`** — Group and period stats; `run_*/figs/stats_only/`.
- **`make_stats_behavior_simple.m`** — Simple behavior stats.

## Prerequisite

- **Step 02** (longitudinal CSV).

## Advanced variant

- **`make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203_cursor.m`** — Newer longitudinal QC (add to this folder when available).
