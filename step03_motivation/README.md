# Step 03: Motivation (PR) — downstream

Run **after Step 02**. Trial- and session-level motivation (Progressive Ratio) from the latest `run_*/ALL_mice_longitudinal.csv`.

## Scripts

- **`run_motivation_analysis.m`** / **`run_motivation_analysis_new.m`** — Trial/session motivation tables and plots. Writes to `run_*/figs/motivation/` (e.g. `TrialsMotiv.csv`, `SessionsMotiv.csv`, PNGs).
- **`motivation_extras_independent_new.m`** — Extra motivation analyses; outputs under `run_*/figs/motivation/`.

## Prerequisite

- **Step 02** must be run first (latest `run_*/ALL_mice_longitudinal.csv`).
