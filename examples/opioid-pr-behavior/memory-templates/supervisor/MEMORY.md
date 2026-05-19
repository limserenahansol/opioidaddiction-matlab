# @supervisor — Opioid AMM project memory

## Vision

Fit **Addiction Motivational Model (AMM)** to run_009 PR lick data with Eli (MNLE/SBI). Behavior-only Tier 1 first.

## Canonical repos

- **Spec:** https://github.com/limserenahansol/opioid-model_wiki — `model/08` (priority), `model/07` (data rules)
- **MATLAB:** https://github.com/limserenahansol/opioidaddiction-matlab — pipeline steps 01–07

## Current priority (Tier 1)

M0 drift `x_t` → M1 `λ=softplus(x)`. Do not expand to M2 until identifiability plan exists (pupil, pair, passive V≈0 prior).

## Delegation

- Implementation → `@analysis-implementer`
- Methods prose → `@paper-writer`
- Citation claims → `@literature-curator`
- Harsh review of passive/PIT wording → `@reviewer`

## Open decisions

- [ ] M2 identifiable with pupil from MATLAB repo?
- [ ] `src/` lives in wiki vs separate fit repo?
- [ ] Eli GitHub username for collaborator invite

## Never

- Invent data paths or mouse IDs not in manifest
- Skip lockout mask for passive During 6–10
