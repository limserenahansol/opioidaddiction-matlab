# Step 06: Dashboards and event-locked visualization — downstream

Run **after Step 02**. Passive/active dashboards, PR+pupil rasters, event-locked pupil, arranged figures. Uses latest `run_*/ALL_mice_longitudinal.csv` (and aligned pupil from Step 01a where needed).

## Scripts (run one per purpose)

| Purpose | Run **one of** |
|---------|-----------------|
| Passive vs active dashboard | `analyze_passive_active_dashboard1.m`, `analyze_passive_active_dashboardnew.m`, `analyze_passive_active_dashboardnew_dec.m` |
| PR + pupil rasters | `analyze_PR_pupil_pairs_raster.m`, `analyze_PR_pupil_pairs_rasternew.m` |
| Event-locked pupil | `analyze_pupil_event_locked.m` |
| Multi-panel arranged figures | `arranged_plots_all_ids.m`, `arranged_plots_all_ids_new.m`, `arranged_plots_all_ids_new2.m` |

**Prerequisite:** Step 02 (longitudinal CSV). For pupil scripts, Step 01a done for sessions of interest. Dashboard/rasters variants also in **step07_advanced**.
