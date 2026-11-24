"""straub_auto_v9.py, v9.1 (2026-11-23): the v9 measure (straub_angle_v9.py) and the v9 decision over time (decide());
frozen after the independent review (user: use the current version from now on). History:
straub_auto_v8.py (2026-11-22) = straub_auto.py with the v8 measure (straub_angle_v8.py) and frames with a
vertical body (rearing / grooming) kept out of the counts (user 2026-11-22). Original docstring:
Automatic Straub tail from a side-view video (added 2026-11-21; user: "could you be able to detect straub tail from
video? if the tail is more than 25 degree from body", "you can train with this" (2026-11-21_09-18-25.mp4), "you can
just crop only mouse part and then also increase contrast").

Per frame: straub_angle.measure() (background subtraction, crop around the mouse at full resolution, contrast scaled to
the mouse's own darkness, body = thick part, tail = longest thin piece from the body's rear end; angle of the tail's
first part above the body axis, or above the horizontal when the body is round or upright). Straub = angle > 25 deg and
the tail rises > 10 deg above the horizontal. Then over time: majority of the measured frames within +-7 frames
(0.25 s); frames with no measurement take the neighbours' state for gaps up to 0.5 s.
Outputs in --out: <video>_frames.csv (per frame), <video>_summary.json, <video>_trace.png.
Not validated against hand scores yet (previous cohort on drive K:, or this cohort's day 10 hand scores).
usage: python straub_auto.py <video> [--out <folder>] [--seconds 300]
"""
import argparse
import json
import sys
import time
from pathlib import Path

import cv2
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import straub_angle_v9 as SA  # noqa: E402  (v9 measure, 2026-11-23)

WIN, GAP = 15, 15
FLIP = 0.75   # v9: a measured frame keeps its own call unless >= 75 % of the measured frames around it disagree
VERT, FLOOR = -1, -2   # state codes: 1 Straub, 0 no, -1 not counted (vertical body), -2 not counted (tail on the floor)


def state(r) -> float:
    if r.get("excluded"):
        return FLOOR if "floor" in str(r.get("why", "")) else VERT
    if r.get("ok"):
        return 1.0 if r["straub"] else 0.0
    return np.nan


def decide(st: np.ndarray, win: int = WIN, flip: float = FLIP, gap: int = GAP) -> np.ndarray:
    """v9 decision over time (user answers 2026-11-23: the v8 0.5-s majority turned 2 correct Straub frames into 'no').
    Measured Straub / no frames keep their own call unless >= flip of the Straub / no frames within +-win//2 disagree;
    not-counted frames stay not counted; an unmeasured frame takes 'not counted' when not-counted frames are at least
    as many as the Straub / no frames around it, else the Straub / no majority; still empty -> the nearest decided
    state within gap frames."""
    st = np.asarray(st, float)
    k = np.ones(win)

    def wsum(x):     # centred +-win//2 sums that keep the input length (review 2026-11-23: 'same' broke clips < 15 frames)
        return np.convolve(x, k, "full")[win // 2: win // 2 + len(x)]
    cS = wsum((st == 1).astype(float))
    cN = wsum((st == 0).astype(float))
    cV = wsum((st == VERT).astype(float))
    cF = wsum((st == FLOOR).astype(float))
    out = np.full(len(st), np.nan)
    tot = np.maximum(cS + cN, 1e-9)
    isS, isN = st == 1, st == 0
    out[isS] = np.where(cN[isS] / tot[isS] >= flip, 0.0, 1.0)
    out[isN] = np.where(cS[isN] / tot[isN] >= flip, 1.0, 0.0)
    out[st == VERT], out[st == FLOOR] = VERT, FLOOR
    un = np.isnan(st)
    cX = cV + cF
    # v9.1 (2026-11-23, independent review, confirmed): unmeasured frames vote with the DECIDED calls of the measured
    # frames, so a call that the 75 % rule overruled cannot still fill its neighbours
    dS = wsum((out == 1).astype(float))
    dN = wsum((out == 0).astype(float))
    xs = un & (cX > 0) & (cX >= dS + dN)
    out[xs] = np.where(cF[xs] >= cV[xs], FLOOR, VERT)
    sn = un & ~xs & (dS + dN > 0)
    out[sn] = np.where(dS[sn] >= dN[sn], 1.0, 0.0)
    return pd.Series(out).ffill(limit=gap).bfill(limit=gap).to_numpy()


def run(video: Path, out: Path, seconds: float | None):
    cap = cv2.VideoCapture(str(video))
    n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    fps = cap.get(cv2.CAP_PROP_FPS) or 30.0
    if seconds:
        n = min(n, int(round(seconds * fps)))
    bg = SA.background(cap, n)
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
    rows, t0 = [], time.time()
    for k in range(n):
        ok, f = cap.read()
        if not ok:
            break
        r = SA.measure(SA.gray_full(f), bg)
        row = dict(frame=k + 1, t_s=round(k / fps, 3), measured=bool(r.get("ok")), why=r.get("why", ""),
                   excluded=bool(r.get("excluded", False)), state_frame=state(r))
        if r.get("ok"):
            row.update(angle_deg=round(r["angle"], 1), rise_deg=round(r["rise"], 1), reference=r["ref"],
                       straub_frame=bool(r["straub"]), tail_len_px=round(r["tail_len"], 1))
        if "floor_gap" in r:
            row.update(floor_gap_H=round(r["floor_gap"] / r["Hb"], 3), lift=round(r["lift"], 3))
        if r.get("excluded") and "rise" in r:      # v9.1: not-counted rows keep the tail's rise (for checking)
            row.update(rise_deg=round(r["rise"], 1))
        rows.append(row)
    D = pd.DataFrame(rows)
    s = D.get("straub_frame", pd.Series(np.nan, index=D.index)).map({True: 1.0, False: 0.0})
    D["state"] = decide(D.state_frame.to_numpy())
    D["straub"] = D.state.where(D.state >= 0)        # 1 / 0 counted; vertical and tail-on-floor frames not counted
    ex = D.state < 0
    D["excluded"] = ex
    D["not_counted_why"] = np.select([D.state == VERT, D.state == FLOOR], ["vertical body", "tail on / near the floor"],
                                     "")
    known = D.state.isin([1.0, 0.0, FLOOR])          # the tail state is known (vertical-body frames are not)
    stem = video.stem
    out.mkdir(parents=True, exist_ok=True)
    D.to_csv(out / f"{stem}_frames.csv", index=False)
    cov = D.straub.notna()
    S = dict(video=str(video), frames=len(D), fps=fps, seconds=round(len(D) / fps, 1),
             measured_pct=round(100 * D.measured.mean(), 1),
             excluded_vertical_pct=round(100 * float((D.state == VERT).mean()), 1),
             excluded_floor_pct=round(100 * float((D.state == FLOOR).mean()), 1),
             counted_frames=int(cov.sum()),
             covered_after_smoothing_pct=round(100 * float(D.state.notna().mean()), 1),
             straub_pct=round(100 * float(D.straub[cov].mean()), 1) if cov.any() else None,
             # the previous cohort's way: the hand score counted every frame (key held while the tail is up), so a
             # tail on the floor counts as Straub-free; vertical-body frames stay out (tail not measured there)
             straub_pct_all_frames=round(100 * float((D.state[known] == 1).mean()), 1) if known.any() else None,
             straub_pct_frame_only=round(100 * float(s.dropna().mean()), 1) if s.notna().any() else None,
             straub_pct_30deg=round(100 * float(((D.angle_deg > 30) & (D.rise_deg > 10)).sum() / max(1, D.measured.sum())), 1),
             not_measured=D.loc[~D.measured, "why"].value_counts().to_dict(),
             rule="v9.1 (2026-11-23): v8h + tail on / near the floor (lowest point within 0.30 body heights of the "
                  "feet, unless clearly raised: > 25 deg and >= 0.3 body heights above its root; legs / forward "
                  "pieces left to v8) = not counted; long tail over the top of the body = Straub; body line tilted "
                  "> 25 deg = floor reference; a measured frame keeps its call unless >= 75 % of the frames within "
                  "0.25 s disagree; unmeasured frames take the decided calls around them",
             version="v9.1",
             seconds_run=round(time.time() - t0))
    (out / f"{stem}_summary.json").write_text(json.dumps(S, indent=1), encoding="utf-8")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(12, 3.4))
    on = D.straub == 1
    for a, b in zip(np.flatnonzero(np.diff(np.r_[0, on.astype(int)]) == 1), np.flatnonzero(np.diff(np.r_[on.astype(int), 0]) == -1)):
        ax.axvspan(D.t_s[a], D.t_s[b], color="#f4c7c3", lw=0)
    ax.plot(D.t_s, D.get("angle_deg"), lw=0.4, color="0.55")
    med = D.set_index("t_s").get("angle_deg").rolling(30, center=True, min_periods=10).median()
    ax.plot(med.index, med.values, lw=1.4, color="#1f2933", label="tail angle, 1-s median")
    ax.axhline(25, color="#c0392b", ls="--", lw=1, label="25 deg")
    ax.set_ylim(-90, 120)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("tail angle above body (deg)")
    ax.set_title(f"{stem}: Straub tail {S['straub_pct']} % of the time (pink = Straub), angle measured in "
                 f"{S['measured_pct']} % of frames", loc="left", fontsize=10)
    ax.legend(loc="lower right", fontsize=8, frameon=False)
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    fig.tight_layout()
    fig.savefig(out / f"{stem}_trace.png", dpi=120)
    plt.close(fig)
    return S


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("video", type=Path)
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "results")
    ap.add_argument("--seconds", type=float, default=None, help="only the first N s (e.g. 300 = the 5-min observation)")
    a = ap.parse_args(argv)
    print(json.dumps(run(a.video, a.out, a.seconds), indent=1))


if __name__ == "__main__":
    main()
