"""Methods compared only on frames where the two labellers agree (same pupil, or both 'not visible')."""
import sys
import numpy as np, pandas as pd, scipy.io as sio
from pathlib import Path
CAND, FIT, VAR = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
L = pd.read_csv(FIT / "labels_long.csv")
masks = {}
for f in sorted(CAND.glob("*_cands.mat")):
    key = f.name.replace("_cands.mat", ""); a = sio.loadmat(f, squeeze_me=True, struct_as_record=False)["allc"]
    for j, c in enumerate(a):
        cs = c if isinstance(c, np.ndarray) else [c]
        for i, x in enumerate(cs):
            if hasattr(x, "m"): masks[(key, j, i + 1)] = x.m.astype(bool)
def same(s, f, a, b):
    if a == b: return True
    if a < 1 or b < 1: return False
    A, B = masks.get((s, f, a)), masks.get((s, f, b))
    if A is None or B is None: return False
    iou = (A & B).sum() / max(1, (A | B).sum()); ya, xa = np.nonzero(A); yb, xb = np.nonzero(B)
    return iou >= 0.3 or B[int(round(ya.mean())), int(round(xa.mean()))] or A[int(round(yb.mean())), int(round(xb.mean()))]
agreed = set()
for (s, f), g in L.groupby(["session", "frame"]):
    a, b = g.choice.iloc[0], g.choice.iloc[-1]
    if same(s, f, a, b) or (a == -1 and b == -1): agreed.add((s, int(f)))
V = pd.read_csv(VAR)
V["agreed"] = [(s, int(f)) in agreed for s, f in zip(V.session, V.frame)]
D = pd.read_csv(FIT / "candidates_scored.csv")
def rule(col, thr):
    out = []
    for r in V.itertuples():
        g = D[(D.session == r.session) & (D.frame == r.frame)]; pick = 0
        if len(g) and np.isfinite(g[col].max()) and g[col].max() >= thr: pick = int(g.loc[g[col].idxmax(), "cand"])
        if r.label >= 1: out.append("miss" if pick == 0 else ("correct" if g.loc[g.cand == pick, "pos"].max() == 1 else "wrong"))
        elif r.label == -1: out.append("false" if pick else "none_ok")
        else: out.append("blob_unknown" if pick else "none_unknown")
    return out
for col, thr in [("s_lr", 0.6), ("s_lr", 0.8), ("s_gbm", 0.6), ("s_h3", 0.05)]:
    V[f"{col}>={thr}"] = rule(col, thr)
A = V[V.agreed]
print(f"agreed frames {len(A)} of {len(V)} (pupil on contour {int((A.label >= 1).sum())}, not visible {int((A.label == -1).sum())})")
cols = [c for c in V.columns if c not in ("session", "frame", "label", "agreed")]
print(f"{'method':16s} {'correct':>7s} {'wrong':>6s} {'miss':>5s} {'false':>6s} {'none-ok':>8s} | recall  bad-pick rate")
for c in cols:
    v = A[c].value_counts(); nc, nw, nm, nf, no = (v.get(k, 0) for k in ("correct", "wrong", "miss", "false", "none_ok"))
    print(f"{c:16s} {nc:7d} {nw:6d} {nm:5d} {nf:6d} {no:8d} | {nc / max(1, (A.label >= 1).sum()):6.1%}  {(nw + nf) / max(1, nc + nw + nf):6.1%}")
