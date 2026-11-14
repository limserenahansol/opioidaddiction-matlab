# Step 08 — PR trial licking trajectory & trajectory phenotype (downstream)

**Primary entry point (recommended):** `run_pr_licking_full_analysis.m`  
- **`USE_DEMO = false` (default):** loads the latest `longitudinal_outputs/run_*/figs/lick_MEGA/per_trial_features.csv`, enriches **between-trial** lick rate from `ALL_mice_longitudinal.csv`, and plots **within-** and **between-trial** lick rate vs **PR trial index** for:
  - **All** data pooled, and
  - Each **phase** (Pre, During, Post, Withdrawal, Reexposure),
  with **individual mouse** traces (grey) and **group** mean ± SEM (all groups in the table: Active/Passive, or e.g. G1/G2 in a toy CSV).  
- Writes **`PR_lick_stats_trialwise.csv`** (group comparison per PR level) and **`PR_lick_slopes_by_mouse.csv`**, plus slope comparison boxplots when two groups are present.  
- Set **`OVERRIDE_RUN`** to a specific `run_*` folder if you do not want “latest by date”.

**Legacy / minimal:** `run_pr_trial_licking_trajectory.m` (single demo or one custom CSV, fewer facets).

---

This step adds **downstream** analyses inspired by *longitudinal trajectory phenotyping* in the fentanyl escalation literature, adapted to a **PR + licking** (no lever-press) design.

It does **not** replace Step 04 (lick extraction) or Step 07 (EFA, LME, etc.). It **consumes** trial-level summaries and produces:

1. **Within-trial vs between-trial lick rate** as a function of PR trial index (PR1, PR2, …).
2. (Optional) **K-means on z-scored trial-by-trial curves** as a simple stand-in for group-based trajectory models (e.g. PROC TRAJ–style) when sample size allows.

## When to use

- After you have (or can export) a table with one row per **(mouse, day, session, PR trial)** and lick rates.
- To compare with **EFA trajectories** from Step 07: trajectory clusters here are based on **raw lick/PR structure**; EFA is based on **latent factors** first—complementary, not identical.

## Input table (CSV)

**Required columns**

| Column | Meaning |
|--------|--------|
| `mouse_key` | Mouse ID (string) |
| `day_index` | Day / run index (numeric) |
| `session_idx` | Session within day if needed (numeric) |
| `Trial` | PR trial index: 1 = first ratio step, 2 = second, … |
| `WithinTrialLickRate` | Licks per second during the trial reward window (definition must match your extraction) |

**Optional**

| Column | Meaning |
|--------|--------|
| `BetweenTrialLickRate` | Licks per second in the inter-trial interval after the previous trial (NaN/empty for `Trial==1` is OK) |
| `Group` | e.g. `A` / `B` for coloring |

If `BetweenTrialLickRate` is missing, the script only plots **within-trial** curves and warns.

## How to run

### A) Demo (no lab paths)

```matlab
cd('.../opioidaddiction-matlab/step08_pr_trajectory_downstream');
run_pr_trial_licking_trajectory
```

`USE_DEMO` is `true` by default; outputs go to `output_demo/`.

### B) Your data

1. In `run_pr_trial_licking_trajectory.m`, set `USE_DEMO = false`.
2. Set `CUSTOM_CSV` to your exported trial table (same columns as above).
3. Run the script.

### C) Optional clustering

After you have a CSV (demo or real):

```matlab
run_pr_trajectory_phenotype_cluster
```

Tweak `K` and `RATE_KIND` (`within`, `between`, or `both`) in the file.

## Outputs

- `GroupMean_LickRate_by_PRLevel.csv` (and per-mouse summary)
- `WithinTrial_LickRate_by_PRLevel.png` (+ `.fig` if `savefig` available)
- `BetweenTrial_LickRate_by_PRLevel.png` (if between column present)
- `pr_trajectory_cluster_assignments.csv` (clustering script)

## Linking to Step 04

If your lick pipeline writes `per_trial_features.csv` or `per_trial_rates.csv` under `figs/lick_patterns_MASTER/`, map columns to the schema above. You may need a small **adapter script** to convert `lick_rate` naming or trial windows into **within** vs **between** using your event timestamps—Step 08 assumes those rates are **already** defined per trial.

## References (conceptual)

- Neuropsychopharmacology-style escalation designs often use **trajectory groups** on repeated drug self-administration; you adapt the same *idea* to **PR lick trajectories** when IV SA is not in the task.
- This repo’s Step 07 EFA scripts remain the place for **factor-score** longitudinal modeling.
