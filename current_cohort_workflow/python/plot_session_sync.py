"""One-session sync check: pupil, licks and rewards on the camera-frame time axis of Step 02's
ALL_mice_longitudinal.csv (or one combined_pupil_digital.csv), plus count checks against the raw Saleae export and the
task log.

    python plot_session_sync.py --longitudinal <run_###>\\ALL_mice_longitudinal.csv --mouse 9669_green2 --day 1
        [--digital <salae>\\day1\\9669m\\green2\\digital.csv] [--jsonl <python>\\day1\\9669m\\green2\\x.jsonl] --out fig.png

Panels: A whole session (pupil + lick raster + reward ticks), B a 60-s zoom at frame resolution, C reward-locked pupil
(baseline-subtracted, mean +- SEM) and lick rate, D the numbers (frames vs camera TTLs; frame-sampled licks vs raw
Saleae contacts vs task-log licks; rewards Injector_TTL vs Saleae pulses vs task log).
Definitions follow the previous pipeline: a lick / reward is a 0->1 change of Lick_TTL / Injector_TTL between frames,
time = CamTime_rel_s.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

COLS = ["mouse_key", "day_index", "session_idx", "Frame", "CamTime_rel_s", "Diameter_px", "Lick_TTL", "Injector_TTL",
        "Trial", "TrialRequirement", "Session_Paradigm", "isPassive", "RequirementLast"]


def onsets(t, x):
    x = (np.nan_to_num(np.asarray(x, float)) > 0.5).astype(int)
    return t[np.nonzero(np.diff(np.r_[0, x]) == 1)[0]]


def load(a):
    if a.combined:
        T = pd.read_csv(a.combined)
        T["mouse_key"], T["day_index"], T["session_idx"] = a.mouse, a.day, 1
        return T
    keep = []
    for ch in pd.read_csv(a.longitudinal, usecols=lambda c: c in COLS, chunksize=1_000_000, low_memory=False):
        m = (ch.mouse_key.astype(str) == a.mouse) & (ch.day_index == a.day)
        if m.any():
            keep.append(ch[m])
        elif keep:
            break
    if not keep:
        raise SystemExit(f"{a.mouse} day {a.day} not in {a.longitudinal}")
    T = pd.concat(keep)
    return T[T.session_idx == T.session_idx.min()]


def raw_counts(a):
    out = {}
    if a.digital and Path(a.digital).exists():
        D = pd.read_csv(a.digital)
        t = D.iloc[:, 0].to_numpy(float)
        col = {c.strip().lower(): c for c in D.columns}
        cam = col.get("camera", col.get("channel 6")); lick = col.get("lick", col.get("channel 7"))
        rew = col.get("reward", col.get("channel 4", col.get("injector")))
        up = lambda c: np.nonzero(np.diff(np.r_[0, D[c].to_numpy(int)]) == 1)[0]
        out["camera TTLs (Saleae)"] = len(up(cam))
        cu = t[up(cam)]
        lu = t[up(lick)]
        out["lick contacts (Saleae raw, camera on)"] = int(np.sum((lu >= cu[0]) & (lu <= cu[-1])))
        if rew:
            xr = D[rew].to_numpy(int)
            ui = np.nonzero(np.diff(np.r_[0, xr]) == 1)[0]; di = np.nonzero(np.diff(np.r_[xr, 0]) == -1)[0]
            ru, wid = t[ui], (t[np.minimum(di, len(t) - 1)] - t[ui])[:len(ui)]
            real = wid >= 0.005            # < 5 ms = line glitch (no Saleae glitch filter on this channel)
            on = (ru >= cu[0]) & (ru <= cu[-1])
            out["reward pulses (Saleae, >= 5 ms, camera on)"] = int(np.sum(real & on))
            if np.sum(~real):
                out["reward-line glitches < 5 ms (not rewards)"] = int(np.sum(~real))
            if np.sum(real & ~on):
                out["reward pulses outside camera time (e.g. valve test)"] = int(np.sum(real & ~on))
    if a.jsonl and Path(a.jsonl).exists():
        L = [json.loads(s) for s in open(a.jsonl, encoding="utf-8-sig") if s.strip()]
        lk = [r for r in L if r.get("tag") in ("lick", "observed_lick")]
        out["licks (task log, all)"] = len(lk)
        out["licks (task log, valid)"] = sum(1 for r in lk if r.get("valid", True))
        out["rewards (task log)"] = sum(1 for r in L if r.get("tag") == "reward_cmd")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--longitudinal", type=Path); ap.add_argument("--combined", type=Path)
    ap.add_argument("--mouse", required=True); ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--digital"); ap.add_argument("--jsonl"); ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--zoom-start", type=float, default=None, help="s; default = 30 s before the 5th reward")
    a = ap.parse_args(argv)
    T = load(a).sort_values("Frame")
    t = T.CamTime_rel_s.to_numpy(float)
    pup = T.Diameter_px.to_numpy(float)
    lick_t, rew_t = onsets(t, T.Lick_TTL), onsets(t, T.Injector_TTL)
    good = np.isfinite(pup)
    smooth = pd.Series(pup).rolling(15, center=True, min_periods=5).median().to_numpy()

    fig = plt.figure(figsize=(15, 10))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.1, 1.0, 1.1])
    ax = fig.add_subplot(gs[0, :])
    ax.plot(t / 60, smooth, color="#3060a0", lw=0.8, label="pupil diameter (px, 0.5-s median)")
    lo, hi = np.nanpercentile(pup, [1, 99]) if good.any() else (0, 1)
    ax.vlines(lick_t / 60, hi + 0.04 * (hi - lo), hi + 0.12 * (hi - lo), color="k", lw=0.3, label="lick (Lick_TTL 0->1)")
    ax.vlines(rew_t / 60, hi + 0.14 * (hi - lo), hi + 0.24 * (hi - lo), color="#d04020", lw=1.0, label="reward (Injector_TTL 0->1)")
    ax.set_xlabel("session time (min, CamTime_rel_s)"); ax.set_ylabel("pupil (px)")
    ax.set_title(f"A  {a.mouse} day {a.day}: whole session ({T.Session_Paradigm.iloc[0]}, isPassive={int(T.isPassive.iloc[0])})", loc="left")
    ax.legend(loc="lower right", fontsize=8, ncol=3)

    z0 = a.zoom_start if a.zoom_start is not None else (max(0, rew_t[4] - 30) if len(rew_t) > 4 else 0)
    zm = (t >= z0) & (t <= z0 + 60)
    ax = fig.add_subplot(gs[1, :])
    ax.plot(t[zm], pup[zm], ".", ms=2, color="#3060a0", label="pupil per frame")
    ax2 = ax.twinx()
    ax2.step(t[zm], T.Lick_TTL.to_numpy(float)[zm] * 0.8, where="post", color="k", lw=0.6, label="Lick_TTL")
    ax2.step(t[zm], T.Injector_TTL.to_numpy(float)[zm] * 1.0, where="post", color="#d04020", lw=1.0, label="Injector_TTL")
    ax2.set_ylim(-0.05, 3); ax2.set_yticks([])
    ax.set_xlabel("session time (s)"); ax.set_ylabel("pupil (px)")
    ax.set_title(f"B  60-s zoom from {z0:.0f} s at frame resolution (one point = one camera frame)", loc="left")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    win = np.arange(-5, 10.01, 1 / 30)
    snips = []
    for r in rew_t:
        if r + win[0] < t[0] or r + win[-1] > t[-1]:
            continue
        y = np.interp(r + win, t[good], pup[good]) if good.sum() > 10 else np.full(win.size, np.nan)
        base = np.nanmean(y[(win >= -2) & (win < 0)])
        snips.append(y - base)
    ax = fig.add_subplot(gs[2, 0])
    if snips:
        S = np.array(snips); mu = np.nanmean(S, 0); se = np.nanstd(S, 0) / np.sqrt(len(S))
        ax.fill_between(win, mu - se, mu + se, color="#3060a0", alpha=0.3); ax.plot(win, mu, color="#3060a0")
    ax.axvline(0, color="#d04020", lw=1); ax.axhline(0, color="0.6", lw=0.5)
    ax.set_xlabel("time from reward (s)"); ax.set_ylabel("pupil change (px, -2..0 s baseline)")
    ax.set_title(f"C  reward-locked pupil (n = {len(snips)} rewards)", loc="left")
    ax = fig.add_subplot(gs[2, 1])
    edges = np.arange(-5, 10.01, 0.25)
    cnt = np.zeros(edges.size - 1)
    for r in rew_t:
        cnt += np.histogram(lick_t - r, edges)[0]
    ax.bar(edges[:-1], cnt / max(1, len(rew_t)) / 0.25, width=0.25, align="edge", color="0.3")
    ax.axvline(0, color="#d04020", lw=1)
    ax.set_xlabel("time from reward (s)"); ax.set_ylabel("licks / s (mean over rewards)")
    ax.set_title("C  reward-locked licking", loc="left")

    ax = fig.add_subplot(gs[2, 2]); ax.axis("off")
    rc = raw_counts(a)
    lines = [f"frames in Step 02 table: {len(T)}",
             f"pupil tracked frames: {good.sum()} ({100 * good.mean():.1f} %)",
             f"licks, frame-sampled (old pipeline): {len(lick_t)}",
             f"rewards, Injector_TTL onsets: {len(rew_t)}",
             f"RequirementLast: {T.RequirementLast.max():g}   max Trial: {T.Trial.max():g}"]
    lines += [f"{k}: {v}" for k, v in rc.items()]
    if "lick contacts (Saleae raw, camera on)" in rc and rc["lick contacts (Saleae raw, camera on)"]:
        lines.append(f"frame-sampling loss: {100 * (1 - len(lick_t) / rc['lick contacts (Saleae raw, camera on)']):.1f} %")
    ax.text(0, 1, "D  checks\n\n" + "\n".join(lines), va="top", family="monospace", fontsize=8.5)
    fig.tight_layout()
    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=130)
    print("saved", a.out); print("\n".join(lines))


if __name__ == "__main__":
    main()
