# Step 04: Lick pipeline — downstream

Run **after Step 02**. Lick mega-pipeline and session-level lick patterns (PCA/k-means). Uses latest `run_*/ALL_mice_longitudinal.csv`.

## Scripts (run one per purpose)

| Purpose | Run **one of** |
|---------|-----------------|
| Full lick pipeline | `run_lick_mega_pipeline.m`, `run_lick_mega_pipeline_new.m` |
| Lick follow-up analyses | `run_lick_MEGA_followup2.m`, `run_lick_MEGA_followup2_new.m` |
| PCA/k-means on session lick features | `analyze_lick_patterns_MASTER.m` |
| Rebuild lick-pattern plots with mouse IDs | `rebuild_lick_patterns_with_ID.m` |

## Prerequisite

- **Step 02** (longitudinal CSV).
