"""Same frames, same criteria: fine-tuned U-Net / original U-Net (eval_fold*.csv) vs the candidate-score rules
(candidates_scored.csv: H3 rule, LR and GBM leave-one-mouse-out). Label-0 frames (pupil visible, not on a contour)
are reported apart, as in the U-Net evaluation.
usage: python -I compare_methods.py <unet_ft dir> <fit dir>"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
U, F = Path(sys.argv[1]), Path(sys.argv[2])
E = pd.concat([pd.read_csv(f) for f in sorted(U.glob("eval_fold*.csv"))], ignore_index=True)
D = pd.read_csv(F / "candidates_scored.csv")
same = {tuple(k.split("|")): set(v) for k, v in json.load(open(F / "same_pupil_sets.json")).items()} if (F / "same_pupil_sets.json").exists() else None


def rule_outcome(col, thr):
    out = []
    for r in E.itertuples():
        g = D[(D.session == r.session) & (D.frame == r.frame)]
        pick = 0
        if len(g) and np.isfinite(g[col].max()) and g[col].max() >= thr: pick = int(g.loc[g[col].idxmax(), "cand"])
        if r.label >= 1:
            ok = g.loc[g.cand == pick, "pos"].max() if pick else 0
            out.append("miss" if pick == 0 else ("correct" if ok == 1 else "wrong"))
        elif r.label == -1: out.append("false" if pick else "none_ok")
        else: out.append("blob_unknown" if pick else "none_unknown")
    return out


rows = {"fine-tuned U-Net": E.ft, "original U-Net": E.orig}
for col, thr, nm in [("s_h3", 0.05, "current hybrid rule"), ("s_lr", 0.6, "LR score >= 0.6"), ("s_lr", 0.8, "LR score >= 0.8"),
                     ("s_gbm", 0.6, "GBM score >= 0.6")]:
    if col in D: rows[nm] = pd.Series(rule_outcome(col, thr))
lab = E.label.to_numpy()
print(f"{'method':22s} {'correct':>8s} {'wrong':>6s} {'miss':>5s} {'false':>6s} {'none-ok':>8s} | recall  bad-pick rate (wrong+false)/(picks with a known answer)")
for nm, v in rows.items():
    v = pd.Series(np.asarray(v)); c = v.value_counts()
    picks = c.get("correct", 0) + c.get("wrong", 0) + c.get("false", 0)
    print(f"{nm:22s} {c.get('correct', 0):8d} {c.get('wrong', 0):6d} {c.get('miss', 0):5d} {c.get('false', 0):6d} {c.get('none_ok', 0):8d} | "
          f"{(v[lab >= 1] == 'correct').mean():6.1%}  {(c.get('wrong', 0) + c.get('false', 0)) / max(1, picks):6.1%}")
