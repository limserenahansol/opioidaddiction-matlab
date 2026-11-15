"""Overlap-aware fit/evaluation of the pupil candidate score (leave-one-MOUSE-out).
A pick is correct when it is the labelled candidate or the same pupil with another boundary (IoU >= 0.3 with the
labelled contour, or centroid inside it and diameter ratio 0.6-1.6). Training positives = that same-pupil set.
usage: python -I fit_pupil_rule3.py <cand dir> <fit dir (labels_consensus.csv)> <diag dir>"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd, scipy.io as sio
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.preprocessing import StandardScaler

CD, FD, DIAG = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
Cn = pd.read_csv(FD / "labels_consensus.csv")
D = pd.read_csv(CD / "candidates.csv")
KEYS = [m["key"] for m in json.load(open(DIAG / "meta.json"))]
masks = {}
for f in sorted(CD.glob("*_cands.mat")):
    key = f.name.replace("_cands.mat", "")
    a = sio.loadmat(f, squeeze_me=True, struct_as_record=False)["allc"]
    for j, c in enumerate(a):
        cs = c if isinstance(c, np.ndarray) else [c]
        for i, x in enumerate(cs):
            if hasattr(x, "m"): masks[(key, j, i + 1)] = x.m.astype(bool)
same = {}
for r in Cn.dropna(subset=["label"]).itertuples():
    lab = int(r.label)
    if lab < 1 or (r.session, int(r.frame), lab) not in masks: continue
    L = masks[(r.session, int(r.frame), lab)]; ld = 2 * np.sqrt(L.sum() / np.pi); s = {lab}; i = 1
    while (r.session, int(r.frame), i) in masks:
        if i != lab:
            M = masks[(r.session, int(r.frame), i)]
            iou = (M & L).sum() / max(1, (M | L).sum())
            my, mx = np.nonzero(M); md = 2 * np.sqrt(M.sum() / np.pi)
            if iou >= 0.3 or (L[int(round(my.mean())), int(round(mx.mean()))] and 0.6 <= md / ld <= 1.6): s.add(i)
        i += 1
    same[(r.session, int(r.frame))] = s
print(f"labelled pupil frames {len(same)}; with other same-pupil candidates: {sum(len(v) > 1 for v in same.values())}")
print("fit of labelled contour:", Cn[Cn.label >= 1].fit.value_counts().to_dict())

sz = {k: tuple(sio.loadmat(DIAG / f"{k}.mat")["roi_in_crop"].ravel()[2:4].astype(float)) for k in KEYS}
D["dfrac"] = D.diam / D.session.map(lambda s: min(sz[s]))
D["unet"] = (D.src == "unet").astype(float); D["stab_f"] = D.stab.fillna(0); D["has_glint"] = D.dglint.notna().astype(float)
DG_FILL = float(D.dglint.median()); D["ldg"] = np.log1p(D.dglint.fillna(DG_FILL))
g = D.groupby(["session", "frame"])
D["con_rel"] = D.con / g.con.transform("max").clip(lower=1e-6); D["dark_rel"] = D.dark - g.dark.transform("min")
D["n_cand"] = g.cand.transform("size"); D["area_rank"] = g.area.rank(pct=True); D["con_rank"] = g.con.rank(pct=True)
D["mouse"] = D.session.str.replace(r"^day\d+_", "", regex=True)
D = D.merge(Cn[["session", "frame", "label"]], on=["session", "frame"], how="left")
D["pos"] = [int(c in same.get((s, f), set())) for s, f, c in zip(D.session, D.frame, D.cand)]
FEAT = ["unet", "con", "sec1", "sec2", "sec3", "sec4", "circ", "sol", "ecc", "dark", "stab_f", "ldg", "has_glint", "dctr",
        "dfrac", "unet_ov", "con_rel", "dark_rel", "con_rank", "area_rank"]
D[FEAT] = D[FEAT].fillna(0)
frames = [(r.session, int(r.frame), int(r.label)) for r in Cn.dropna(subset=["label"]).itertuples()]


def evaluate(col, thr, name, quiet=False):
    out = []
    for s, f, lab in frames:
        gg = D[(D.session == s) & (D.frame == f)]
        pick = 0
        if len(gg) and np.isfinite(gg[col].max()) and gg[col].max() >= thr: pick = int(gg.loc[gg[col].idxmax(), "cand"])
        if lab >= 1: o = "miss" if pick == 0 else ("correct" if pick in same.get((s, f), {lab}) else "wrong")
        else: o = "correct_none" if pick == 0 else "false_pick"
        out.append(dict(session=s, frame=f, label=lab, pick=pick, outcome=o))
    E = pd.DataFrame(out); c = E.outcome.value_counts(); npk = int((E.pick > 0).sum()); bad = c.get("wrong", 0) + c.get("false_pick", 0)
    if not quiet:
        print(f"{name:36s} correct {c.get('correct', 0):3d} none-ok {c.get('correct_none', 0):3d} WRONG {c.get('wrong', 0):3d} "
              f"FALSE {c.get('false_pick', 0):3d} miss {c.get('miss', 0):3d} | recall {(E[E.label >= 1].outcome == 'correct').mean():.1%} "
              f"| bad {bad}/{npk} = {bad / max(1, npk):.1%}")
    return E, bad / max(1, npk), int(c.get("correct", 0))


def h3(gg):
    sc = pd.Series(-np.inf, index=gg.index)
    u = gg[(gg.src == "unet") & (gg.circ >= 0.6) & (gg.con >= 0.08)]
    if len(u): sc.loc[u.index] = 10 + u.con * u.circ * u.sol; return sc
    dk = gg[(gg.src == "dark") & (gg.circ >= 0.6) & (gg.sol >= 0.8) & (gg.con >= 0.12)]
    prox = np.where(dk.dglint.notna(), np.maximum(0.35, np.exp(-np.maximum(0, dk.dglint.fillna(0) - 1.5) / 3)), 1)
    sc.loc[dk.index] = dk.con * dk.circ * dk.sol * np.minimum(1, dk.stab / 3) * prox
    return sc


D["s_h3"] = D.groupby(["session", "frame"], group_keys=False)[D.columns.tolist()].apply(h3)
evaluate("s_h3", 0.05, "current hybrid rule (H3)")
D["s_unet"] = np.where(D.src == "unet", D.area, -np.inf); evaluate("s_unet", 0, "U-Net only (largest blob)")
oracle = sum(1 for s, f, lab in frames if lab >= 1)
print(f"best possible: {oracle} frames have the pupil on a candidate, {sum(1 for _, _, l in frames if l < 1)} frames have none")


def lomo(model_fn, col):
    D[col] = np.nan
    for m in sorted(D.mouse.unique()):
        tr, te = D.mouse != m, D.mouse == m
        mdl = model_fn(); mdl.fit(D.loc[tr, FEAT], D.loc[tr, "pos"])
        D.loc[te, col] = mdl.predict_proba(D.loc[te, FEAT])[:, 1]


class LR:
    def __init__(self, C=0.5): self.sc = StandardScaler(); self.m = LogisticRegression(C=C, max_iter=5000, class_weight="balanced")
    def fit(self, X, y): self.m.fit(self.sc.fit_transform(X), y); return self
    def predict_proba(self, X): return self.m.predict_proba(self.sc.transform(X))


for name, fn, col in [("LR", lambda: LR(0.5), "s_lr"),
                      ("GBM (upper bound)", lambda: HistGradientBoostingClassifier(max_depth=3, learning_rate=0.08, max_iter=200, class_weight="balanced"), "s_gbm")]:
    lomo(fn, col)
    print(f"--- {name}, leave-one-mouse-out")
    for thr in [0.3, 0.5, 0.6, 0.7, 0.8, 0.9]:
        evaluate(col, thr, f"{name} thr {thr}")
E, _, _ = evaluate("s_lr", 0.6, "LR thr 0.6 (per session below)", quiet=True)
print(E.groupby("session").outcome.value_counts().unstack(fill_value=0).to_string())
D.to_csv(FD / "candidates_scored.csv", index=False)
