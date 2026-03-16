# Step 07: Advanced pipeline (EFA, modules 5–12, Straub, predictive)

Run **after Step 02**. Advanced analyses: feature QC, GLMM/LME, PCA/clustering/EFA, event-locked, cumulative fit, cross-modal, RL model, predictive/decoder, Straub tail, addiction score, longitudinal QC, rasters.

## Plan / module → MATLAB implementation

| Plan / Module | MATLAB implementation |
|---------------|------------------------|
| **Module 5** (Feature QC) | `analyze_modules_5_to_11` |
| **Module 6** (GLMM/LME) | `analyze_modules_5_to_11` |
| **Module 7** (PCA, clustering, EFA) | `analyze_modules_5_to_11` |
| **Module 8** (Event-locked) | `analyze_modules_5_to_11` |
| **Module 9** (Cumulative fit) | `analyze_modules_5_to_11` |
| **Module 10** (Cross-modal) | `analyze_modules_5_to_11` |
| **Module 11** (RL model) | `analyze_modules_5_to_11` |
| **Module 12** (Predictive) | `analyze_modules_5_to_11` (or separate `analyze_modules_5_to_12`) |
| **Straub tail** | `compute_straub_tail_only_v1` |
| **Dashboard / preprocessing** | `analyze_passive_active_dashboard_dec2` |
| **Addiction score / EFA** | `analyze_addiction_score_efa_*.m` |
| **Longitudinal QC** | `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203_cursor.m` |
| **Rasters** | `plotLickAndBoutRasters_SelectedDays.m` |

Add the corresponding `.m` files to this folder as they are implemented. The same logic is available in the Python repo (decoder, EFA, cross-generalization) — see that repo for complementary analyses.

## Prerequisite

- **Step 02** (longitudinal CSV and run folder).
