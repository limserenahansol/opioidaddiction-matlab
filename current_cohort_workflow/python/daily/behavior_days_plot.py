#!/usr/bin/env python3
"""
behavior_days_plot.py - behaviour of every mouse across days, NO pupil (added 2026-11-20, user: "just except for pupil
show me the data day 34567 data"). Report only; reads the routine's tables, changes nothing.

Panels (first 15 min 10 s of each session, the same window as the daily summary):
  PR (last ratio requirement) · rewards · licks per min · licks per reward · time to reward (median s) · tail immersion
Passive mice on yoked days: rewards = drops delivered (the partner's replay); their lick-line signals are NOT shown
(the passive spout is inside the mouth, so those signals are not licks). PR / licks / time to reward on yoked days are
the active mice only. Days 3-5 every mouse did PR (Baseline, natural reward).
Inputs: analysis/session_metrics.csv (window first_15m10s), assays/tail_immersion_all.csv, cohort_manifest groups via
pr_table_days*.csv (group per mouse). Outputs: analysis/days/behavior_days<first>-<last>.png + .csv (long table).
usage: python behavior_days_plot.py [--days 3 7]
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
GCOL = {"Active": "#c0392b", "Passive": "#2c6fbb"}          # as cohort_plots.py / ti_table_days.py (validated pair)
INK, MUTED, GRID, SHADE = "#1f2933", "#7b8794", "#e4e7eb", "#fdf0e2"
PHASES = [(1, 5, "Baseline · natural reward"), (6, 10, "Yoked · morphine"), (11, 13, "Equal access"),
          (14, 16, "Withdrawal"), (17, 18, "Re-exposure")]


def groups() -> dict:
    """mouse key 'cage_mouse' -> Active / Passive from the newest PR table (the cohort manifest groups)"""
    fs = [f for f in (P / "analysis" / "days").glob("pr_table_days1-*.csv") if f.stem.split("-")[-1].isdigit()]
    fs = sorted(fs, key=lambda f: int(f.stem.split("-")[-1]))       # skips kept copies like ..._before_<date>_...
    G = pd.read_csv(fs[-1], encoding="utf-8-sig")
    return {f"{r.cage}_{r.mouse}": r.group for r in G.drop_duplicates(["cage", "mouse"], keep="last").itertuples()}


def table(d0: int, d1: int) -> pd.DataFrame:
    M = pd.read_csv(P / "analysis" / "session_metrics.csv")
    M = M[(M["window"] == "first_15m10s") & M["day"].between(d0, d1)].copy()
    gr = groups()
    M["group"] = M["mouse"].map(gr)
    passive = M["task"].str.upper().eq("PASSIVE")
    out = pd.DataFrame(dict(day=M["day"], mouse=M["mouse"], group=M["group"], sex=M["sex"], task=M["task"],
                            rewards=M["rewards"],
                            pr=M["pr_breakpoint"].where(~passive),
                            licks_per_min=M["licks_per_min"].where(~passive),
                            licks_per_reward=M["licks_per_reward"].where(~passive),
                            time_to_reward_s=M["time_to_reward_median_s"].where(~passive)))
    TI = pd.read_csv(P / "assays" / "tail_immersion_all.csv", encoding="utf-8-sig")
    TI["mouse"] = TI["cage"].astype(str) + "_" + TI["mouse"].astype(str)
    TI = TI.drop_duplicates(["day", "mouse"], keep="last")[["day", "mouse", "latency_s"]]
    out = out.merge(TI.rename(columns={"latency_s": "tail_immersion_s"}), on=["day", "mouse"], how="left")
    return out.sort_values(["group", "mouse", "day"]).reset_index(drop=True)


def short(m: str) -> str:
    cage, col = m.split("_")
    return f"{'f' if cage.endswith('f') else 'm'} {col}"


def panel(ax, D, col, title, unit, d0, d1, label_ends=False):
    for lo, hi, name in PHASES:
        a, b = max(lo, d0) - 0.5, min(hi, d1) + 0.5
        if a < b:
            if lo >= 6:
                ax.axvspan(a, b, color=SHADE, zorder=0, lw=0)
            ax.text((a + b) / 2, 1.0, name, transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                    fontsize=7.5, color=MUTED)
    ends = []
    for m, g in D.groupby("mouse"):
        g = g.dropna(subset=[col]).sort_values("day")
        if g.empty:
            continue
        grp = g["group"].iloc[0]
        ls = "-" if str(g["sex"].iloc[0]).upper().startswith("F") else "--"
        ax.plot(g["day"], g[col], ls, marker="o", ms=3.5, lw=1.0, color=GCOL.get(grp, MUTED), alpha=0.45, zorder=2)
        ends.append((g["day"].iloc[-1], g[col].iloc[-1], short(m), GCOL.get(grp, MUTED)))
    for grp, c in GCOL.items():
        s = D[D["group"] == grp].groupby("day")[col].median().dropna()
        n = D[D["group"] == grp].groupby("day")[col].count()
        s = s[n.reindex(s.index) >= 2]
        if len(s):
            ax.plot(s.index, s.values, "-o", color=c, lw=2.4, ms=6, zorder=3, label=f"{grp} median",
                    markeredgecolor="white", markeredgewidth=1.2)
    if label_ends and ends:
        # direct labels at each line's last day, pushed apart vertically so they do not overlap
        lo_y, hi_y = ax.get_ylim()
        gap = (hi_y - lo_y) * 0.055
        for x in sorted({e[0] for e in ends}):
            grp_ends = sorted([e for e in ends if e[0] == x], key=lambda e: e[1])
            ys = []
            for e in grp_ends:
                y = e[1] if not ys else max(e[1], ys[-1] + gap)
                ys.append(y)
            for e, y in zip(grp_ends, ys):
                ax.text(x + 0.12, y, e[2], fontsize=7, color=INK, va="center")
                ax.plot([x + 0.02, x + 0.1], [e[1], y], color=e[3], lw=0.6, alpha=0.6)
    ax.set_title(f"{title}" + (f"  ({unit})" if unit else ""), fontsize=10, color=INK, loc="left", pad=14)
    ax.set_xticks(range(d0, d1 + 1))
    ax.set_xlim(d0 - 0.4, d1 + 0.75)
    ax.grid(axis="y", color=GRID, lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(MUTED)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.set_xlabel("day", fontsize=8, color=MUTED)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", nargs=2, type=int, default=[3, 7])
    a = ap.parse_args(argv)
    d0, d1 = a.days
    D = table(d0, d1)
    out = P / "analysis" / "days" / f"behavior_days{d0}-{d1}"
    D.to_csv(out.with_suffix(".csv"), index=False)
    fig, axs = plt.subplots(2, 3, figsize=(15.5, 8.6))
    spec = [("pr", "PR = last ratio requirement", "", True), ("rewards", "Rewards", "first 15 min 10 s", False),
            ("licks_per_min", "Licks per min", "active / PR sessions", False),
            ("licks_per_reward", "Licks per reward", "", False),
            ("time_to_reward_s", "Time to reward, median", "s", False),
            ("tail_immersion_s", "Tail immersion latency", "s", False)]
    for ax, (col, title, unit, lab) in zip(axs.ravel(), spec):
        panel(ax, D, col, title, unit, d0, d1, label_ends=lab)
    h, l = axs[0, 0].get_legend_handles_labels()
    h += [plt.Line2D([], [], color=MUTED, lw=1, ls="-", marker="o", ms=3.5, alpha=0.6),
          plt.Line2D([], [], color=MUTED, lw=1, ls="--", marker="o", ms=3.5, alpha=0.6)]
    l += ["one mouse, female (9670f)", "one mouse, male (9669m)"]
    fig.legend(h, l, loc="upper right", ncol=4, fontsize=8.5, frameon=False, bbox_to_anchor=(0.99, 0.995))
    fig.suptitle(f"Behaviour per mouse, days {d0}-{d1} (no pupil)", x=0.01, ha="left", fontsize=13, color=INK,
                 fontweight="bold", y=0.995)
    fig.text(0.01, 0.005,
             "Active = self-administering mice (9669m black, 9669m green2, 9670f green2, 9670f white); Passive = their "
             "yoked partners (days 6+: rewards = drops delivered by the replay of the partner's session). Days 3-5 every "
             "mouse did PR.\nPassive lick-line signals on yoked days are not shown (spout inside the mouth = not licks). "
             "Thin = one mouse, thick = group median (n >= 2).\n"
             "Sources: analysis\\session_metrics.csv (first 15 min 10 s), assays\\tail_immersion_all.csv (the values you "
             "sent) · table: " + str(out.with_suffix('.csv').relative_to(P)) + " · tools\\behavior_days_plot.py",
             fontsize=7.5, color=MUTED, va="bottom")
    fig.tight_layout(rect=(0, 0.06, 1, 0.955))
    fig.savefig(out.with_suffix(".png"), dpi=130)
    plt.close(fig)
    print(out.with_suffix(".png"))
    print(out.with_suffix(".csv"))


if __name__ == "__main__":
    main()
