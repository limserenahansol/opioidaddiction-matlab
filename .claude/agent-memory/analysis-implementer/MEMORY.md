# @analysis-implementer — Opioid PR / AMM

## Stack (this repo = MATLAB primary)

- **MATLAB:** step01_pupil, step01_manual_assays, step02_build_longitudinal, step03_motivation, step04_lick, step05–07
- **Outputs:** `ALL_mice_longitudinal.csv`, figs under `run_###/`
- **Python / SBI:** https://github.com/limserenahansol/opioid-model_wiki `src/` (when built)
- **Fitting:** MNLE / SBI per Eli — spec in wiki `model/08`

## Data rules (mandatory — model/07)

```python
# day_index >= 4
# passive During 6-10: lick_observed = False (mask), reward_observed = True
# PairID: random effect in hierarchical fit
```

## Tier 1 equations

```
dx/dt = ((α + r(t)) / τ) dt + (σ/τ) dW
λ_t = softplus(x_t)
r = +R rewarded lick, -L unrewarded
disengage if x_t < θ_stop
```

Also implement smoothed: `E[r|T] = R/T - L*(T-1)/T` for Tier 2 comparison.

## Starter sessions

- Active: 6099_orange, day 12 or 13 (Post)
- Passive: 6099_red, day 6–10 (During)

## JSON source

Drive: https://drive.google.com/drive/folders/1tmojU4ahssZEvAdNGa5w6BfPAYERuYsV?usp=drive_link

## Validation

- Within-trial lick slowing reproduced by drift model
- PPC on breakpoint (requirementLast)
- Lockout-off sanity: worse fit if licks forced to zero

## Do not

- Treat x as lick rate
- Fit days 1–3
- Implement full PR schedule POMDP as v1
