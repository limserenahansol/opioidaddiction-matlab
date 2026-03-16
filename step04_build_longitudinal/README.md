# Step 04: Build longitudinal_outputs (run after 01–03)

Run **after Step 01, 02, 03**. This step scans raw data and builds the master longitudinal CSV used by Steps 05–07.

## Script

- **`Longitudinal_final_trialrequire_HOTTST_passive_final_handle_nomatchingstrabu.m`** — Aggregates per-session files (`combined_pupil_digital.csv`/`.xlsx`, JSONL) under `BASE\day1..dayN\<cage>\<mouse>\concat_out_*\` and writes:
  - `BASE\longitudinal_outputs\run_###\ALL_mice_longitudinal.csv`
  - Related outputs in the same `run_###` folder.

## Before running

1. Set **`BASE`** at the top of the script (e.g. `'K:\addiction_concate_Dec_2025'`).
2. Ensure directory layout: `BASE\day1..dayN\<cage>\<mouse>\concat_out_*\` with `combined_pupil_digital.csv` or `.xlsx` and `*.jsonl`.

## Output

- `longitudinal_outputs\run_###\ALL_mice_longitudinal.csv` — one row per session (mouse × day). **Steps 05–07** use this file (latest `run_*` by default).

## Next

Then run **Step 05** (QC/longitudinal), **Step 06** (assays), **Step 07** (visualization) in any order.
