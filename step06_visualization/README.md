# Step 06: Dashboards and event-locked visualization — downstream

Run **after Step 02**. Passive/active dashboards, PR+pupil rasters, event-locked pupil, arranged figures. Uses latest `run_*/ALL_mice_longitudinal.csv` (and aligned pupil from Step 01a where needed).

## Scripts

- **`analyze_passive_active_dashboard1.m`** / **`analyze_passive_active_dashboardnew.m`** / **`analyze_passive_active_dashboardnew_dec.m`** — Passive vs active dashboard.
- **`analyze_PR_pupil_pairs_raster.m`** / **`analyze_PR_pupil_pairs_rasternew.m`** — PR + pupil rasters.
- **`analyze_pupil_event_locked.m`** — Event-locked pupil (aligned pupil + behavioral data).
- **`arranged_plots_all_ids.m`** / **`arranged_plots_all_ids_new.m`** / **`arranged_plots_all_ids_new2.m`** — Multi-panel arranged figures.

## Prerequisite

- **Step 02** (longitudinal CSV). For pupil scripts, **Step 01a** should be done for the sessions of interest.

## Related

- **`analyze_passive_active_dashboard_dec2`** — Dashboard/preprocessing (add to this folder or Step 07 when available).
- **`plotLickAndBoutRasters_SelectedDays.m`** — Rasters (see Step 07 advanced or add here).
