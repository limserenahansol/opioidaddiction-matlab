# Step 02: Build longitudinal (all mice × all days)

Run **after all per-session work (Step 01a pupil + Step 01b manual assays)**. Aggregates raw session data into one folder: **longitudinal_outputs/run_###/** with **ALL_mice_longitudinal.csv**.

## Script

- **`Longitudinal_final_trialrequire_HOTTST_passive_final_handle_nomatchingstrabu.m`** — Scans `BASE\day1..dayN\<cage>\<mouse>\concat_out_*\` (combined_pupil_digital.csv/.xlsx, JSONL, assay inputs) and writes:
  - `BASE\longitudinal_outputs\run_###\ALL_mice_longitudinal.csv`
  - Related outputs in the same `run_###` folder.

## Before running

1. Set **`BASE`** at the top (e.g. `'K:\addiction_concate_Dec_2025'`).
2. Ensure layout: `BASE\day1..dayN\<cage>\<mouse>\concat_out_*\` with `combined_pupil_digital.csv`/`.xlsx` and `*.jsonl` (and assay data as expected by the script).

## Output

- **One folder per run:** `longitudinal_outputs\run_###\` with **ALL_mice_longitudinal.csv** (one row per session: mouse × day). All **downstream** steps (03–07) use this file (latest `run_*` by default).

## Next

Run **Steps 03–07** (downstream): motivation, licking, reward, pupil, delta, statistics, plots, pharmacology, and advanced pipeline (EFA, modules 5–12, etc.).
