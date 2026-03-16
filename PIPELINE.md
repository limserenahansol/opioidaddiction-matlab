# Pipeline: Run order and roadmap

**Run order:** Do all **per-session** work (Step 01a + 01b), then **Step 02** (build longitudinal). Then run **downstream** steps 03–07 in any order.

---

## Flow

```mermaid
flowchart TB
    S01a[Step 01a: Pupil per session]
    S01b[Step 01b: Manual assays per session]
    S02[Step 02: Build longitudinal]
    CSV[("ALL_mice_longitudinal.csv")]
    S03[Step 03: Motivation]
    S04[Step 04: Lick]
    S05[Step 05: QC and stats]
    S06[Step 06: Visualization]
    S07[Step 07: Advanced pipeline]
    S01a --> S02
    S01b --> S02
    S02 --> CSV
    CSV --> S03
    CSV --> S04
    CSV --> S05
    CSV --> S06
    CSV --> S07
```

| Step | Description |
|------|-------------|
| 01a | Pupil tracking + alignment (per session/video) |
| 01b | Manual scoring: TST, HOT, Straub, etc. (per session); results feed into Step 02 |
| 02 | Build one folder: `longitudinal_outputs/run_###/`, **ALL_mice_longitudinal.csv** (all mice × all days) |
| 03 | Motivation (PR) — downstream |
| 04 | Lick pipeline — downstream |
| 05 | QC, longitudinal stats/plots — downstream |
| 06 | Dashboards, rasters, event-locked pupil — downstream |
| 07 | Advanced: EFA, modules 5–12, Straub, addiction score, predictive — downstream |

---

## Dependencies

- **01a, 01b:** Per-session; set paths in each script. Do for all sessions before Step 02.
- **02:** Set `BASE` in `step02_build_longitudinal/Longitudinal_final_trialrequire_HOTTST_passive_final_handle_nomatchingstrabu.m`.
- **03–07:** Use the **latest** `run_*` under `longitudinal_outputs` (script logic or path at top).

---

## Advanced pipeline (Step 07)

Scripts in **step07_advanced/** are often triplicated/optimized variants — **run one per purpose.** See step07_advanced/README.md for the full table.

| Purpose | Run **one of** |
|---------|-----------------|
| Modules 5–12 | `analyze_modules_5_to_11.m`, `analyze_modules_5_to_11_new.m`, or `analyze_advanced_pipeline.m` |
| Straub tail | `compute_straub_tail_only_v1.m`, `compute_straub_nonmoving_only_v1.m` |
| Dashboard | `analyze_passive_active_dashboard_dec2.m`, `analyze_passive_active_dashboard_dec.m` |
| Addiction score / EFA | one of `analyze_addiction_score_efa_dec*.m`, `_both.m`, `_mousefit.m` |
| Longitudinal QC | `make_longitudinal_QC_and_requested_analyses_NEWCOHORT_20260203_cursor.m` or other QC variant |

Complementary EFA/decoder/cross-generalization in the Python repo.
