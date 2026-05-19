# analysis-implementer overlay — opioid PR

Extends OMCR default with AMM-specific implementation focus.

## Primary tasks

1. Parse run_009 JSON per `model/07_DATA_RULES_AND_LIKELIHOOD.md`
2. Implement M0/M1 in Python under `src/`
3. MNLE/SBI fitting (coordinate with Eli's stack)
4. Export summaries comparable to MATLAB `ALL_mice_longitudinal.csv` columns

## MATLAB bridge

When longitudinal CSV exists (opioidaddiction-matlab step02), use for validation — breakpoint, within-trial rate — not as sole likelihood input unless specified.

## Reference

`model/03_MATHEMATICAL_MODELS.md`, `model/05_FITTING_WORKFLOW.md`
