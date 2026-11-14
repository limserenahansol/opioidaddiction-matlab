# Step 9 — Paper-oriented trajectory & benchmark (hub)

This is an **organizational** step: it does **not** re-implement EFA or the full lick MEGA pipeline. It links the analyses you listed for a methods/benchmark section and notes what is **TBD** (e.g. neural / imaging / transcriptomic prediction).

**Step 10** is not a separate code folder. Optional extensions (e.g. outcome prediction) are listed as **TBD** below and belong in the same private analysis stream once outcome tables exist.

---

## Map: your 4 “benchmark” blocks → this repo

| Your step | What to run / where | Status in repo |
|-----------|----------------------|----------------|
| **1) EFA addiction index → trajectory / clustering** | `step07_advanced/`: e.g. `compute_addiction_index_EFA_mousefit_projectDays*.m`, `analyze_addiction_score_efa_*.m` | Implemented in **Step 07** (not duplicated here). |
| **2) “Paper-like” trajectory: 3 (or 3–5) behavioral streams** | Minimal set aligned with *your* data: **PR breakpoint**, **within-trial lick rate**, **between-trial lick** (and optionally ITI / total lick). | **Step 08** `run_pr_licking_full_analysis.m` + `run_pr_trajectory_phenotype_cluster.m` (k-means on z-scored trajectories) — *not* lever-press SA; PR/lick analog as discussed. |
| **3) Compare the two clusterings (EFA vs raw-trajectory)** | Join on `mouse_key`, then overlap metrics | **`compare_trajectory_clusterings.m`** in this folder: confusion matrix, **ARI**, **NMI**; optional CSV/figure export. |
| **4) Clusters → neural / imaging / transcriptomic outcomes** | Project-specific; needs outcome tables and governance | **TBD** — not automated here. After you have e.g. region response or DE scores per mouse, add a small GLM or correlation script under this folder or `step07_advanced/`. |

---

## Variable mapping (no lever-press; PR / lick design)

| Role | Possible operational measures in *your* pipeline |
|------|--------------------------------------------------|
| Motivation / effort | PR breakpoint, PR-derived score, trial count in session |
| Consummatory / “within” engagement | Lick rate within trial (`lick_rate` in MEGA) |
| Persistent / inter-trial seeking | Between-trial lick rate (Step 08 enrichment) |
| Extra (if not redundant) | ITI lick frequency, total lick, etc. (watch r > 0.8) |

**Naming in prose:** e.g. “PR-based licking trajectory phenotyping” (not a drug-escalation label unless your design matches).

---

## One-click (recommended)

From MATLAB, after **Step 7** (EFA / `AI_mouse` table exists) and **Step 8** (`run_pr_licking_full_analysis` with `USE_DEMO = false` so `output_full/MouseMean_LickRate_by_PRLevel.csv` exists):

```matlab
cd('<path to opioidaddiction-matlab>\step09_paper_trajectory_benchmark');
run_step9_one_click
```

This script:

1. Resolves the **same latest `run_*`** as Step 8 (`pr_step08_locate_data`).
2. Loads EFA / addiction summary from the first path that exists, e.g.  
   `figs/addiction_index_EFA/AI_mouse_summary.csv` or  
   `figs/addiction_score_efa_mousefit/addiction_score_mousefit_by_mouse.csv`.
3. Runs **k-means on `AI_mouse`** (or `AI_mouse_z`) → `out_step9_auto/efa_kmeans_for_step9.csv`.
4. Builds **PR-trajectory clusters** from Step 8’s `MouseMean...` (same k-means as `run_pr_trajectory_phenotype_cluster`) → `out_step9_auto/pr_trajectory_for_step9.csv`.
5. Calls **`compare_trajectory_clusterings`** → `out_step9_auto/compare/` (ARI, NMI, confusion).

Edit `K_EFA`, `K_PR`, and `RATE_PR` at the top of `run_step9_one_click.m` if needed.

---

## How to run the comparison (step 3) — manual two CSVs

1. Export or save two **mouse-level** label tables (same `mouse_key`):
   - **EFA-based** cluster: from Step 07 outputs (you may need a one-line export to CSV).
   - **Raw-trajectory** cluster: e.g. `pr_trajectory_cluster_assignments.csv` from Step 08 or a variant with your 3–5 feature stack.
2. From this directory:

   ```matlab
   compare_trajectory_clusterings('path/to/efa_clusters.csv', 'path/to/pr_clusters.csv', ...
       'outDir', fullfile(pwd,'out_compare'));
   ```

   Each CSV should include `mouse_key` and a `cluster` (or `trajectory_cluster`) column; see the script’s header for column names it accepts.

---

## Relation to the public GitHub “benchmark” repo

The **synthetic** benchmark materials live in a **separate** public repo (toy data + generic ML). Step 9 here is the **in-lab, paper-structured** counterpart and should stay **private** until you choose what to share.

---

## Files

- `compare_trajectory_clusterings.m` — ARI, NMI, confusion matrix, matched mice
- (Optional) you can add your own `predict_outcomes_from_clusters.m` when outcome vectors are ready.
