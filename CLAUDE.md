# Claude Code — Morphine PR MATLAB pipeline

> OMCR: [omcr/INSTALL.md](./omcr/INSTALL.md) · **AMM spec:** [opioid-model_wiki](https://github.com/limserenahansol/opioid-model_wiki)

## Project context

- **Working title:** Morphine progressive-ratio (PR) behavioral pipeline + longitudinal analysis
- **Field:** Behavioral neuroscience (head-fixed mouse operant)
- **Target venue:** [TBD — linked to AMM modeling paper / grant]
- **Central hypothesis:** Same as AMM wiki — latent motivational state explains PR across active/passive and phases; this repo implements **data processing and phenotypes** feeding that model.
- **Research topic:** Pipeline steps 01a–07: pupil, manual assays, longitudinal CSV, motivation, lick, QC, visualization, advanced (EFA, Straub, addiction score).
- **Datasets:** Local `BASE` path (e.g. `K:\addiction_concate_Dec_2025` or run_009); per-session `combined_pupil_digital.csv`, `*.jsonl`; output `ALL_mice_longitudinal.csv`.
- **Narrative spine:** Phase 1 per-session → Phase 2 longitudinal → Phase 3 downstream (any order). See `PIPELINE.md`, `README.md`.

## Research stack

- **Pipeline docs:** `README.md`, `PIPELINE.md`, each `step*/README.md`
- **Computational model (canonical):** https://github.com/limserenahansol/opioid-model_wiki — `model/07`, `model/08`, `model/05`
- **Logic figures:** `pipeline_schematic.png`, `addiction_rl_logic_flowchart.png`
- **Wiki stub:** `wiki/README.md` → TRAP_analysis research hub
- **BibTeX:** `paper/references.bib` · **Summary:** `references.csv`
- **Python modeling (future):** may live in opioid-model_wiki `src/` or Eli's fit repo

## Language preference

- **Code comments / README:** English
- **Internal notes:** Korean OK
- **Do not** commit absolute paths like `K:\...` — use config at script top only locally

## OMCR agent routing

| Task | Agent |
|------|--------|
| Which step to run, pipeline order | `@supervisor` |
| MATLAB scripts, longitudinal build, lick/motivation | `@analysis-implementer` |
| Methods describing pipeline | `@paper-writer` |
| Dashboards, rasters, flowcharts | `@figure-descriptor` |
| QC / passive lockout in data | `@reviewer` |
| Citations | `@literature-curator` |

## Analysis rules (align with wiki model/07)

- Analyze `day_index >= 4` in downstream stats
- Passive During days 6–10: lick TTL invalid / masked in model fits; rewards yoked — still valid
- `PairID` for paired active/passive inference
- Prefer **second day of each phase** for stable summaries (4, 7, 12, 15, 18) when phase means are noisy

## Hard rules

1. Run Phase 1 (01a, 01b) for all sessions before step 02
2. One script per purpose in step07 — see README "run one of" tables
3. Model fitting rules live in **opioid-model_wiki**, not invented here
4. No PII in committed paths or manifests
