# Step 05: QC and longitudinal stats/plots — downstream

Run **after Step 02**. Quality control, requested analyses, and longitudinal plots/stats from the latest `run_*/ALL_mice_longitudinal.csv`.

## Scripts (run one per purpose)

| Purpose | Run **one of** |
|---------|-----------------|
| QC + requested analyses | `make_longitudinal_QC_and_requested_analyses_cursor.m`, `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_v5.m` |
| Longitudinal plots by week/group + stats | `make_longitudinal_plotsall_statistic.m`, `make_longitudinal_plotsall_statistic_NEWCOHORT_2026.m` |
| Group and period stats from CSV | `make_stats_only_from_csv.m`, `make_stats_only_from_csv_dec.m` |
| Simple behavior stats | `make_stats_behavior_simple.m` |

**Prerequisite:** Step 02 (longitudinal CSV). For newer longitudinal QC variants see **step07_advanced** (e.g. `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203_cursor.m`).
