# Step 01b: Manual scoring — pharmacological assays (per-session)

**Per-session** work — run before building the longitudinal file. Manual scoring for hotplate (HP), tail suspension test (TST), Straub tail, etc. Results are then fed into **Step 02** (build longitudinal) so that one folder has all mice x all days including assay data.

## Scripts (run one per purpose)

| Purpose | Run **one of** |
|---------|-----------------|
| TST/HOT/Straub/tail summaries | `make_tst_hot_tail_simple.m`, `make_tst_hot_tail_simple_new.m` |
| Batch manual video scoring (TST, HOT, etc.) | `manual_scoring_video_multibehavior_batch2.m` (set folder path inside) |

## Next

After **Step 01a (pupil)** and **Step 01b (manual assays)** are done for all sessions, run **Step 02** to build **longitudinal_outputs** (one folder, all mice x all days).
