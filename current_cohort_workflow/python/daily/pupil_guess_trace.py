"""Best-guess pupil trace for hard sessions (narrowed / inflamed eye), next to the tracked pupil. Never replaces it.

The hybrid detector writes hybrid_detections.csv for every frame: the best pupil candidate (diameter, centre) and its
score (0-1). The tracked pupil (pupil_diameter_FINAL_U-Net_results2.csv) only uses detections with score >= the
model threshold (qc.json "hybridThreshold"). Here the bar is lower (--min-score, default 0.25): candidates between
--min-score and the threshold are kept as a GUESS, then median-filtered (--smooth-s) and gaps up to --max-gap-s are
interpolated. Output per session (in its pupil_outputs_general folder):
  pupil_guess.csv   frame, t_s, guess_diam_px, score, kind (tracked-quality / guess / none)
  pupil_guess.png   tracked diameter vs guess, and the score over time
Use the guess only with care (label it as such); for analyses keep the tracked pupil.

usage: python pupil_guess_trace.py --pupil-root <pupil_tst_final> [--days 1 2] [--min-score 0.25]
"""
import argparse, json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def med_filter(x, k):
    if k < 3: return x
    s = pd.Series(x)
    return s.rolling(k, center=True, min_periods=max(1, k // 3)).median().to_numpy()


def session_grade(sd: Path) -> str:
    """grade of this session in <pupil root>/pupil_session_quality.csv ('' if none)."""
    root = sd.parent.parent.parent.parent
    pq = root / "pupil_session_quality.csv"
    try:
        day = int(sd.parent.parent.parent.name[3:])
        Q = pd.read_csv(pq, dtype={"cage": str, "mouse": str})
        r = Q[(Q["day"].astype(int) == day) & (Q["cage"].astype(str) == sd.parent.parent.name) &
              (Q["mouse"].astype(str).str.lower() == sd.parent.name.lower())]
        return str(r.iloc[-1]["grade"]) if len(r) else ""
    except (FileNotFoundError, ValueError, KeyError):
        return ""


def one(sd: Path, min_score, smooth_s, max_gap_s):
    det = sd / "hybrid_detections.csv"
    if not det.exists(): return None
    D = pd.read_csv(det)
    q = json.loads((sd / "qc.json").read_text(encoding="utf-8")) if (sd / "qc.json").exists() else {}
    thr = float(q.get("hybridThreshold", np.nan))
    t = D["t_s"].to_numpy(); sc = D["score"].to_numpy(); d = D["best_diam_px"].to_numpy(dtype=float)
    dt = np.nanmedian(np.diff(t)) if len(t) > 2 else 1 / 30
    keep = np.isfinite(d) & (sc >= min_score)
    g = np.where(keep, d, np.nan)
    g = med_filter(g, int(round(smooth_s / max(dt, 1e-3))) | 1)
    s = pd.Series(g)
    gap = int(round(max_gap_s / max(dt, 1e-3)))
    g = s.interpolate(limit=gap, limit_area="inside").to_numpy()
    kind = np.where(np.isfinite(thr) & (sc >= thr), "tracked-quality", np.where(keep, "guess", "none"))
    # session grade (pupil_session_quality.csv, 2026-11-15): with grade guess / no_pupil the detector is often not on the
    # pupil even when its score is high (narrowed or inflamed eye), so every kept frame of that session is a guess
    grade = session_grade(sd)
    if grade in ("guess", "no_pupil"):
        kind = np.where(keep, "guess", "none")
    out = pd.DataFrame(dict(frame=D["frame"], t_s=t, guess_diam_px=g, score=sc, kind=kind))
    out.to_csv(sd / "pupil_guess.csv", index=False)
    tr = None
    r2 = sd / "pupil_diameter_FINAL_U-Net_results2.csv"
    if r2.exists():
        T = pd.read_csv(r2); col = next((c for c in T.columns if "iam" in c), None)
        if col is not None: tr = T[col].to_numpy(dtype=float)
    fig, ax = plt.subplots(2, 1, figsize=(11, 4.6), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
    if tr is not None: ax[0].plot(t[: len(tr)] / 60, tr[: len(t)], color="0.3", lw=0.8, label="tracked pupil (results2.csv)")
    ax[0].plot(t / 60, g, color="tab:orange", lw=0.8, alpha=0.9, label=f"guess (score >= {min_score})")
    ax[0].set_ylabel("diameter (ROI px)"); ax[0].legend(fontsize=8, loc="upper right")
    ax[1].plot(t / 60, sc, color="tab:blue", lw=0.5)
    if np.isfinite(thr): ax[1].axhline(thr, color="k", lw=0.8, ls="--")
    ax[1].axhline(min_score, color="tab:orange", lw=0.8, ls=":")
    ax[1].set_ylabel("score"); ax[1].set_xlabel("session time (container min)"); ax[1].set_ylim(0, 1)
    frac_t = float(np.mean(kind == "tracked-quality")); frac_g = float(np.mean(kind == "guess"))
    gtxt = f"session grade '{grade}': " if grade else ""
    fig.suptitle(f"{sd.parent.parent.parent.name} {sd.parent.parent.name} {sd.parent.name}: {gtxt}tracked-quality "
                 f"{100 * frac_t:.0f} % of frames, guess {100 * frac_g:.0f} % (use the guess with care)", fontsize=10)
    fig.tight_layout(); fig.savefig(sd / "pupil_guess.png", dpi=90); plt.close(fig)
    return dict(session=f"{sd.parent.parent.parent.name} {sd.parent.parent.name}_{sd.parent.name}", grade=grade,
                tracked_quality=frac_t, extra_guess=frac_g, threshold=thr)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pupil-root", required=True); ap.add_argument("--days", type=int, nargs="*")
    ap.add_argument("--min-score", type=float, default=0.25); ap.add_argument("--smooth-s", type=float, default=0.5)
    ap.add_argument("--max-gap-s", type=float, default=1.0)
    a = ap.parse_args(argv)
    root = Path(a.pupil_root); rows = []
    for sd in sorted(root.glob("day*/*/*/pupil_outputs_general")):
        if a.days and int(sd.parent.parent.parent.name[3:]) not in a.days: continue
        r = one(sd, a.min_score, a.smooth_s, a.max_gap_s)
        if r: rows.append(r); print(f"{r['session']:22s} tracked-quality {100 * r['tracked_quality']:5.1f} %  extra guess {100 * r['extra_guess']:5.1f} %")
    if rows:
        new = pd.DataFrame(rows); sp = root / "pupil_guess_summary.csv"
        if sp.exists():  # keep the rows of sessions not redone in this call (2026-11-15)
            old = pd.read_csv(sp); new = pd.concat([old[~old["session"].isin(new["session"])], new]).sort_values("session")
        new.to_csv(sp, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
