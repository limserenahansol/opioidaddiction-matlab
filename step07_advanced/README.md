# Step 07: Advanced pipeline (EFA, modules 5–12, Straub, predictive)

Run **after Step 02**. This folder holds optimized variants of advanced analyses (many scripts are duplicates with small changes — **run one per purpose**).

## Scripts (run one per purpose)

| Purpose | Run **one of** |
|---------|-----------------|
| **Modules 5–12** (QC, GLMM, PCA/EFA, event-locked, cumulative fit, cross-modal, RL, predictive) | `analyze_modules_5_to_11.m`, `analyze_modules_5_to_11_new.m` — or `analyze_advanced_pipeline.m` |
| **Straub tail** | `compute_straub_tail_only_v1.m`, `compute_straub_nonmoving_only_v1.m` |
| **Dashboard / preprocessing** | `analyze_passive_active_dashboard_dec2.m`, `analyze_passive_active_dashboard_dec.m` |
| **Addiction score / EFA** | `analyze_addiction_score_efa_dec4.m`, `analyze_addiction_score_efa_dec3.m`, `analyze_addiction_score_efa_dec2.m`, `analyze_addiction_score_efa_both.m`, `analyze_addiction_score_efa_mousefit.m` — or `compute_addiction_index_EFA_mousefit_projectDays*.m` |
| **Longitudinal QC** | `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203_cursor.m`, `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203.m`, or other `make_longitudinal_QC_*` variants in this folder |
| **Event-locked / rasters** | `analyze_pupil_event_locked.m`, `analyze_pupil_event_locked_REVISED.m`, `analyze_pupil_event_lockednew.m` — or `analyze_PR_pupil_pairs_raster.m` / `analyze_PR_pupil_pairs_rasternew.m` |

**Prerequisite:** Step 02 (longitudinal CSV). Run from this folder or add it to the MATLAB path. Complementary EFA/decoder in the Python repo.
