#!/usr/bin/env python3
"""
wholeeye_compare_models.py - old vs new whole-eye model on the same sessions (added 2026-11-18 for the retrain on 21
user outlines: the 18 of 2026-11-17 plus 3 late day-1 9670f black frames that include the white part; user
2026-11-18: the white/cloudy part IS eye). Reads only eye_track_all.m outputs; changes nothing.

For every session in --new (<root>/day<N>/<cage>/<mouse>/wholeeye_frames.csv) that also exists in --old:
  * eye size = equivalent diameter 2*sqrt(area/pi) px (what the tables carry), per-minute median, old vs new;
  * found %, median size, first vs last 10 min, off-eye % and its worst tenth (same rule as qc_summary.py: mask centre
    more than 20 px from the 1800-frame rolling median centre);
  * the 24 qc frames that eye_track_session.m saves (same frames for both runs): old outline orange, new outline cyan,
    and the overlap (IoU) of the two masks.
Outputs in --out: <day>_<cage>_<mouse>.png per session, compare_traces.png (all sessions), compare_summary.csv.
Time axis: frame k = camera pulse k at 30.26 Hz (not the mp4's 1/60 s stamps).
usage: python wholeeye_compare_models.py [--old analysis/wholeeye] [--new analysis/wholeeye/black3] [--out ...]
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image

P = Path(__file__).resolve().parent.parent
HZ = 30.26
OLD_C, NEW_C = "#ff8c1a", "#00d5ff"


def eqd(area):
    return 2.0 * np.sqrt(np.asarray(area, float) / np.pi)


def offeye(W):
    """same rule as qc_summary.wholeeye_offeye"""
    fnd = W.found.to_numpy() == 1
    cx, cy = np.where(fnd, W.cx, np.nan), np.where(fnd, W.cy, np.nan)
    med = lambda x: pd.Series(x).rolling(1800, center=True, min_periods=100).median().to_numpy()  # noqa: E731
    with np.errstate(invalid="ignore"):
        off = fnd & (np.hypot(cx - med(cx), cy - med(cy)) > 20)
    return 100 * off.mean(), 100 * max(x.mean() for x in np.array_split(off, 10))


def per_minute(W):
    fnd = W.found.to_numpy() == 1
    d = np.where(fnd, eqd(W.area_px.fillna(0)), np.nan)
    m = ((W.frame.to_numpy() - 1) / HZ / 60).astype(int)
    s = pd.Series(d).groupby(m).median()
    return s.index.to_numpy() + 0.5, s.to_numpy()


def stats(W):
    fnd = W.found.to_numpy() == 1
    d = np.where(fnd, eqd(W.area_px.fillna(0)), np.nan)
    t = (W.frame.to_numpy() - 1) / HZ / 60
    oe = offeye(W)
    return dict(found_pct=100 * fnd.mean(), size_med=np.nanmedian(d), size_first10=np.nanmedian(d[t < 10]),
                size_last10=np.nanmedian(d[t >= t.max() - 10]), offeye_pct=oe[0], offeye_worst10=oe[1])


def stretch(img):
    a = img.astype(float)
    lo, hi = np.percentile(a, [2, 98])
    return np.clip((a - lo) / max(1.0, hi - lo), 0, 1)


def qc_frames(folder):
    out = {}
    for f in (folder / "qc").glob("f*_mask.png"):
        k = re.match(r"f(\d+)_mask", f.stem)
        if k:
            out[int(k.group(1))] = f
    return out


def one(key, old_dir, new_dir, out_dir):
    Wo, Wn = pd.read_csv(old_dir / "wholeeye_frames.csv"), pd.read_csv(new_dir / "wholeeye_frames.csv")
    so, sn = stats(Wo), stats(Wn)
    qo, qn = qc_frames(old_dir), qc_frames(new_dir)
    ks = sorted(set(qo) & set(qn))
    ious = []
    fig = plt.figure(figsize=(13.2, 11.2))
    gs = fig.add_gridspec(6, 6, height_ratios=[1.25, 0.18, 1, 1, 1, 1], hspace=0.35, wspace=0.05)
    ax = fig.add_subplot(gs[0, :])
    for W, c, lab in ((Wo, OLD_C, "old model (18 outlines)"), (Wn, NEW_C, "new model (21 outlines)")):
        x, y = per_minute(W)
        ax.plot(x, y, "-o", ms=3, color=c, label=lab)
    ax.set_xlabel("minutes (camera pulses at 30.26 Hz)")
    ax.set_ylabel("eye size, equivalent diameter (px)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_title(f"{key}: per-minute median eye size | old: found {so['found_pct']:.1f} %, median {so['size_med']:.1f} px, "
                 f"off-eye {so['offeye_pct']:.1f} % (worst tenth {so['offeye_worst10']:.1f} %) | new: found "
                 f"{sn['found_pct']:.1f} %, median {sn['size_med']:.1f} px, off-eye {sn['offeye_pct']:.1f} % "
                 f"(worst tenth {sn['offeye_worst10']:.1f} %)", fontsize=8.5)
    for i, k in enumerate(ks[:24]):
        a = fig.add_subplot(gs[2 + i // 6, i % 6])
        img = np.array(Image.open(new_dir / "qc" / f"f{k}_img.png").convert("L"))
        mo = np.array(Image.open(qo[k]).convert("L")) > 0
        mn = np.array(Image.open(qn[k]).convert("L")) > 0
        iou = (mo & mn).sum() / max(1, (mo | mn).sum()) if (mo.any() or mn.any()) else np.nan
        ious.append(iou)
        a.imshow(stretch(img), cmap="gray", vmin=0, vmax=1)
        for m, c in ((mo, OLD_C), (mn, NEW_C)):
            if m.any():
                a.contour(m.astype(float), levels=[0.5], colors=c, linewidths=1.1)
        a.set_title(f"f{k}  {(k - 1) / HZ / 60:.1f} min  IoU {iou:.2f}", fontsize=7.5)
        a.set_xticks([]), a.set_yticks([])
    fig.suptitle(f"{key}: orange = old model, cyan = new model (same 24 frames; contrast stretched)", fontsize=11)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{key}.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    row = dict(session=key, frames_old=len(Wo), frames_new=len(Wn))
    row.update({f"{k}_old": round(v, 2) for k, v in so.items()})
    row.update({f"{k}_new": round(v, 2) for k, v in sn.items()})
    ious = np.array(ious, float)
    row.update(qc_frames=len(ks), qc_iou_mean=round(float(np.nanmean(ious)), 3) if ks else np.nan,
               qc_iou_min=round(float(np.nanmin(ious)), 3) if ks else np.nan)
    return row, (Wo, Wn)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", default=str(P / "analysis" / "wholeeye"))
    ap.add_argument("--new", default=str(P / "analysis" / "wholeeye" / "black3"))
    ap.add_argument("--out", default=str(P / "analysis" / "wholeeye" / "qc" / "compare_black2_black3"))
    a = ap.parse_args(argv)
    old, new, out = Path(a.old), Path(a.new), Path(a.out)
    rows, traces = [], {}
    for f in sorted(new.glob("day*/*/*/wholeeye_frames.csv")):
        day, cage, mouse = f.parts[-4], f.parts[-3], f.parts[-2]
        o = old / day / cage / mouse
        key = f"{day}_{cage}_{mouse}"
        if not (o / "wholeeye_frames.csv").exists():
            print(f"  {key}: no old output - skipped")
            continue
        row, tr = one(key, o, f.parent, out)
        rows.append(row)
        traces[key] = tr
        print(f"  {key}: size {row['size_med_old']:.1f} -> {row['size_med_new']:.1f} px, last 10 min "
              f"{row['size_last10_old']:.1f} -> {row['size_last10_new']:.1f} px, off-eye {row['offeye_pct_old']:.1f} -> "
              f"{row['offeye_pct_new']:.1f} % (worst tenth {row['offeye_worst10_old']:.1f} -> "
              f"{row['offeye_worst10_new']:.1f} %), qc IoU mean {row['qc_iou_mean']}")
    if not rows:
        print("no session to compare")
        return 1
    pd.DataFrame(rows).to_csv(out / "compare_summary.csv", index=False)
    n = len(traces)
    fig, axs = plt.subplots(n, 1, figsize=(11, 2.1 * n), sharex=True, squeeze=False)
    for ax, (key, (Wo, Wn)) in zip(axs[:, 0], traces.items()):
        for W, c in ((Wo, OLD_C), (Wn, NEW_C)):
            x, y = per_minute(W)
            ax.plot(x, y, "-o", ms=2.5, color=c)
        ax.set_ylabel("px", fontsize=8)
        ax.set_title(key, fontsize=9, loc="left")
        ax.grid(alpha=0.3)
    axs[-1, 0].set_xlabel("minutes")
    fig.suptitle("eye size (equivalent diameter, per-minute median): orange = old model, cyan = new model", fontsize=10)
    fig.tight_layout()
    fig.savefig(out / "compare_traces.png", dpi=110)
    plt.close(fig)
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
