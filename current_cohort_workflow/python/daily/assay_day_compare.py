#!/usr/bin/env python3
"""
assay_day_compare.py - TST, hot plate and Straub tail of two (or more) assay days side by side, per mouse (added
2026-11-23; user: "TST and HoT and STraub comparison with day 5 vs day 10"). Report only.
  TST        non-moving %, automatic, first 180 s (analysis/days/dayN/assays_dayN.csv, tst_pct_nonmoving)
  hot plate  non-moving %, the user's hand score, first 30 s (assays_dayN.csv, hot_pct_Non-moving) - the two measures
             the previous pipeline's dashboard shows; the other hot-plate behaviours are in the csv
  Straub     % of the first 5 min, final value (analysis/days/dayN/straub_dayN.csv straub_final_pct: automatic v9.1, or
             the user's hand score over the counting frames), on every day that has one (first Straub day = 9)
One line per mouse (Active red, Passive blue, as cohort_plots.py), thick = group median.
usage: python assay_day_compare.py --days 5 10   -> analysis/days/assay_compare_days5-10.png + .csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
DAYS_DIR = P / "analysis" / "days"
GCOL = {"Active": "#c0392b", "Passive": "#2c6fbb"}
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]


def groups():
    G = pd.read_csv(P / "cohort_manifest.csv", dtype=str).fillna("")
    return {(r.cage_folder, r.mouse): r.group for r in G.itertuples()}


def collect(days):
    G, rows = groups(), []
    for d in sorted(set(days) | set(range(min(days), max(days) + 1))):
        f = DAYS_DIR / f"day{d}" / f"assays_day{d}.csv"
        if d in days and f.exists():
            A = pd.read_csv(f, dtype={"cage": str})
            for _, r in A.iterrows():
                for meas, col in (("TST non-moving (%)", "tst_pct_nonmoving"), ("Hot plate non-moving (%)", "hot_pct_Non-moving")):
                    v = pd.to_numeric(r.get(col), errors="coerce")
                    if pd.notna(v):
                        rows.append(dict(day=d, cage=r.cage, mouse=r.mouse, group=G.get((r.cage, r.mouse), ""), measure=meas,
                                         value=float(v), source=f"analysis/days/day{d}/assays_day{d}.csv"))
        s = DAYS_DIR / f"day{d}" / f"straub_day{d}.csv"
        if s.exists():
            S = pd.read_csv(s, dtype={"cage": str})
            for _, r in S.iterrows():
                v = pd.to_numeric(r.get("straub_final_pct"), errors="coerce")
                if pd.isna(v):
                    v = pd.to_numeric(r.get("auto_straub_pct"), errors="coerce")
                if pd.notna(v):
                    rows.append(dict(day=d, cage=r.cage, mouse=r.mouse, group=G.get((r.cage, r.mouse), ""),
                                     measure="Straub tail (%)", value=float(v), source=f"analysis/days/day{d}/straub_day{d}.csv"))
    return pd.DataFrame(rows)


def panel(ax, D, meas, days, first):
    g = D[D.measure == meas]
    dd = sorted(int(x) for x in g.day.unique()) if len(g) else []
    ax.set_title(meas.replace(" (%)", ""), loc="left", fontsize=12, color="#1f3a5f", fontweight="bold", pad=10)
    if not dd:
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", color="#7b8794")
        ax.set_xticks([]); ax.set_yticks([])
        return
    xs = {d: i for i, d in enumerate(dd)}
    ends = []                                     # (label x, y, text, colour): spread below so they do not overlap
    for c, m in MICE:
        h = g[(g.cage == c) & (g.mouse == m)].sort_values("day")
        if h.empty:
            continue
        col = GCOL.get(h.group.iloc[0], "#7b8794")
        x = [xs[d] for d in h.day] if len(dd) > 1 else [xs[d] + (0.06 if h.group.iloc[0] == "Passive" else -0.06) for d in h.day]
        ax.plot(x, h.value, "-o", color=col, alpha=0.45, lw=1.2, ms=4, zorder=2)
        lx = x[-1] + 0.06 if len(dd) > 1 else xs[h.day.iloc[-1]] + (0.16 if h.group.iloc[0] == "Passive" else -0.32)
        ends.append([lx, float(h.value.iloc[-1]), f"{c[-1]}{m}", col])
    for lx in sorted({round(e[0], 3) for e in ends}):
        L = sorted((e for e in ends if round(e[0], 3) == lx), key=lambda e: e[1])
        for k in range(1, len(L)):
            L[k][1] = max(L[k][1], L[k - 1][1] + 3.6)    # at least 3.6 percentage points apart
        for x_, y_, t_, c_ in L:
            ax.text(x_, y_, t_, fontsize=7, color=c_, va="center",
                    bbox=dict(boxstyle="square,pad=0.05", facecolor="white", edgecolor="none", alpha=0.75))
    for grp, col in GCOL.items():
        M = g[g.group == grp].groupby("day").value.median()
        if len(M):
            ax.plot([xs[d] for d in M.index], M.values, "-o", color=col, lw=3, ms=7, zorder=3, label=f"{grp} median")
    ax.set_xticks(range(len(dd)))
    ax.set_xticklabels([f"day {d}" for d in dd])
    ax.set_xlim(-0.4, len(dd) - 0.4 + 0.35)
    ax.set_ylim(0, 105)
    ax.set_ylabel("% of time")
    ax.grid(axis="y", color="#e4e7eb", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    if first:
        ax.legend(fontsize=8, frameon=False, loc="lower right")
    if meas.startswith("Straub") and set(dd) != set(days):
        ax.set_xlabel("Straub was recorded from day 9 (no day-5 baseline)", fontsize=8.5, color="#52606d")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+", required=True)
    a = ap.parse_args(argv)
    days = sorted(a.days)
    D = collect(days)
    tag = f"days{'-'.join(map(str, days))}"
    out_csv, out_png = DAYS_DIR / f"assay_compare_{tag}.csv", DAYS_DIR / f"assay_compare_{tag}.png"
    D.to_csv(out_csv, index=False, encoding="utf-8-sig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.4))
    for k, meas in enumerate(("TST non-moving (%)", "Hot plate non-moving (%)", "Straub tail (%)")):
        panel(axes[k], D, meas, days, k == 0)
    hand = D[(D.measure == "Straub tail (%)")].copy()
    hand_txt = ""                                 # 2026-11-23 (review): the Straub source is read from the tables
    srcs = []
    for d in sorted(hand.day.unique()) if len(hand) else []:
        S = pd.read_csv(DAYS_DIR / f"day{d}" / f"straub_day{d}.csv", dtype={"cage": str})
        if "straub_final_source" in S:
            srcs += [f"{r.cage} {r.mouse} d{d}" for r in S.itertuples() if str(r.straub_final_source).startswith("hand")]
    if srcs:
        hand_txt = ", your hand score: " + ", ".join(srcs)
    fig.suptitle(f"Assays: day {' vs day '.join(map(str, days))} (TST automatic 3 min · hot plate your hand score 30 s · "
                 f"Straub automatic 5 min{hand_txt})", fontsize=11, color="#1f3a5f", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(out_png, dpi=120)
    plt.close(fig)
    print(out_png)
    print(out_csv)
    for meas in ("TST non-moving (%)", "Hot plate non-moving (%)", "Straub tail (%)"):
        g = D[D.measure == meas]
        if len(g):
            med = g.groupby(["day", "group"]).value.median().round(1)
            print(f"  {meas}: " + "; ".join(f"day {d} {grp} {v:g}" for (d, grp), v in med.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
