# Step 10 — Methods: PR licking trajectory phenotyping

This file is copied to `out_step10/METHODS_step10.md` when you run `run_step10_pr_lick_trajectory_phenotype.m`. Edit the source in `step10_pr_licking_trajectory_phenotype/docs/METHODS_step10.md` if you change the pipeline.

---

## 1. Purpose

Unsupervised **trajectory phenotyping** of mice using only **PR + licking** measures from Step 8’s **per-trial enriched** table (no EFA). The goal is to find subgroups whose **multivariate longitudinal pattern** (across calendar days) differs, in analogy to lever-press trajectory work, but adapted to your task.

---

## 2. Input

- **Default:** `step08_pr_trajectory_downstream/output_full/per_trial_enriched_step08.csv`
- Required columns (at minimum): `mouse_key`, `day_index`, `Trial`, lick-related fields as used in Step 8 (`lick_count`, `WithinTrialLickRate` [from sanitize], `BetweenTrialLickRate`, `iei_median`, `Group`, and for trial-level masks: `Period` or `PeriodClean`).

---

## 3. Day window

- Design: **PR starts on calendar day 3**; analysis uses **`day_index >= DAY_INDEX_MIN`** (default **4**), so **days 1–3** are excluded from Step 10.
- `day_index` is coerced to numeric so filters cannot fail silently on text/categorical imports.

---

## 4. Mouse × day table (`dsum`)

For each **(mouse_key, day_index)** group (all trials that day):

| Variable | Definition |
|----------|------------|
| `pr_breakpoint` | **Max `Trial`** = highest PR ratio step reached that day. |
| `total_lick_count` | **Sum of `lick_count`** over trials. If every `lick_count` is NaN for that day, the sum is **NaN** (not 0). |
| `mean_within_hz` | Mean of **within-trial lick rate** (`WithinTrialLickRate` after Step 10 sanitize caps). |
| `mean_between_hz` | Mean of **`BetweenTrialLickRate`** when present. |
| `mean_iti_lick_freq_hz` | **Mean(1 / `iei_median`)** over trials with valid IEI (seconds), a **within-trial** inter-lick frequency (Hz-like). Not a separate “between-trial ITI lick count” unless you add that in Step 8. |

Written: `out_step10/pr_licking_mouse_by_day.csv`.

---

## 5. Passive / “During” / no-lick days (optional)

**Script constants** (edit in `run_step10_pr_lick_trajectory_phenotype.m`):

- `PASSIVE_DURING_LICK_DAYS` (e.g. **6:10**), `PASSIVE_DURING_PERIOD_NAME` (e.g. **`During`**).
- Set `PASSIVE_DURING_LICK_DAYS = []` to **disable** all of this and use full lick data (e.g. to test whether clusters still separate without missingness rules).

**Trial table `E`:** For rows with **Passive** + **Period** matching + **day in range**, set to **NaN:** `WithinTrialLickRate`, `BetweenTrialLickRate`, `lick_rate`, `lick_count`, **`iei_median`**.

**`dsum`:** For **Passive** + **calendar day in range**, set **`mean_within_hz`**, **`mean_between_hz`**, **`total_lick_count`**, **`mean_iti_lick_freq_hz`** to **NaN**. **PR (`pr_breakpoint`) is not nulled** here.

**Imputation (below):** For **day-trajectory** mode, **Passive** × **those days** × lick feature blocks in the design matrix, missing values are **not** filled with the column mean (they stay NaN through imputation, then become 0 in `Z` after z-scoring). **PR** columns are imputed as usual if missing.

---

## 6. Two analysis modes

| Mode | When | Design matrix `X` |
|------|------|---------------------|
| **daytraj** | At least **`MIN_DAYS`** distinct `day_index` in `dsum` (default 2) and enough mice | One row per mouse. Columns = **stacked blocks** in order: for each of the **five** features, **all calendar days** (≥ `DAY_INDEX_MIN`) in order. So length = **5 × (number of distinct days in data)**. |
| **mouse3** | Too few distinct days | One row per mouse, **5 pooled** variables over all remaining trials: max PR, total licks, mean within, mean between, mean 1/IEI (see code). |

Feature list for trajectory: `STEP10_TRAJ_FNS` in the main script (PR, total licks, within, between, ITI-frequency proxy).

---

## 7. Missing data, imputation, z-score, clustering

1. **Column-wise imputation:** Remaining NaNs (except where **omitted** for Passive lick window in daytraj) are filled with the **column mean** of that feature column across mice.
2. **Z-score** each column of `X` (mean/SD across mice; `omitnan`). Columns with **near-zero variance** are dropped. Residual non-finite values in `Z` are set to **0**.
3. **k-means** on `Z` with **`K_TRAJ`** (default 2), replicates, `kmeans` “plus” starts if available. **Cluster labels 1…K are arbitrary** (not “high vs low” by definition).

**Output:** `out_step10/pr_licking_trajectory_clusters.csv` (`mouse_key`, `trajectory_phenotype_cluster`, `Group` from manifest).

---

## 8. Manifest Group vs cluster

- `Group` is read from the enriched file (e.g. Active / Passive) **per mouse**.
- **Crosstabs:** `pr_cluster_x_group_matrix.csv`, `pr_cluster_x_group_long.csv`, figure `pr_group_vs_trajectory_cluster_heatmap.png`.
- **Trajectory cluster ≠ experimental group** a priori; interpret overlap from these tables.

---

## 9. Figures and tables (typical `out_step10` contents)

| File | Content |
|------|--------|
| `pr_licking_mouse_by_day.csv` | Raw mouse×day features. |
| `pr_feature_corr_mouse_day.csv` / `pr_feature_correlation_heatmap.png` | Pairwise **Pearson** correlations on mouse-day vectors (up to 5 features). |
| `pr_cluster_mean_trajectory_shape.png` | **Raw** (pre–z-score) feature means ± SE by cluster; x = stacked feature×day (or 5 pooled in mouse3). |
| `pr_multifeature_cluster_means_by_day.png` | Mean ± SEM by cluster vs **calendar day** (one panel per feature). |
| `pr_breakpoint_spaghetti_by_cluster.png` | Individual **PR only** vs day, color = trajectory cluster, line style = manifest Active/Passive. |
| `pr_breakpoint_spaghetti_by_cluster_split_Active_Passive.png` | Same, split by manifest group. |
| `METHODS_step10.md` | This document (copy of `docs/METHODS_step10.md`). |

---

## 10. Limitations (brief)

- k-means assumes **spherical** groups in **z-scored** space; sample size and **K** should be justified.
- **ITI** here = **1/median IEI** aggregated within-trial, not a dedicated between-trial ITI lick count unless you add it upstream.
- **Total lick count** is **sum of trial licks**, not normalized by session duration unless you add duration to `dsum`.
- If `PASSIVE_DURING_LICK_DAYS` is empty, passive animals still contribute full lick metrics for all days in range.

---

## 11. 한국어 요약 (짧게)

- **데이터:** Step 8 `per_trial_enriched` 기준, **day 4~**만 사용(실험일 1–3 제외).
- **쥐×날 5가지:** (1) **PR breakpoint**(그날 max trial) (2) **총 핥기**(`lick_count` 합) (3) **trial 내 핥기율** (4) **trial 사이 핥기율** (5) **1/IEI 평균**(trial 내 inter-lick 빈도에 가깝게).
- **Passive + During + 지정일(6–10 등):** 핥기·IEI는 **NaN** 처리( PR은 유지); 군집에 쓰는 lick 쪽은 **impute로 채우지 않고** z 후 0.
- **군집:** 5×날(또는 5요약) 벡터를 **열 z-score** 한 뒤 **k-means**로 2(또는 K)집.
- **Cluster 번호 1/2**는 **임의**이며, **Active/Passive**와 다를 수 있음 → 교차표/heatmap으로 비교.

---

*Document version: aligned with `run_step10_pr_lick_trajectory_phenotype.m` (five-feature trajectory, passive window, copy to out_step10).*
