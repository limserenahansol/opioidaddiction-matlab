#!/usr/bin/env python3
"""
results_plots.py - one picture of the RESULTS so far for the routine (added 2026-11-23; user: "in routine you should
show me the plot of pupil and straub and TI, TST, HOT plate always! the result"). Report only: reads the routine's own
tables and changes nothing.

Six panels, one thin line per mouse (Active red, Passive blue, as cohort_plots.py), thick line = group median, shaded =
morphine days - the same drawing as the tail-immersion plot (ti_table_days.draw_plot):
  pupil size       session median pupil diameter (px), first 15 min 10 s - analysis/master_day_level.csv
                   (quick_analysis.py, the pupil the tables use); the 2 whole-eye mice (eye damage) and sessions graded
                   guess / no_pupil are left out, as in the previous cohort's PU06
  pupil at reward  pupil change after reward (% of the -1..0 s baseline), same table and the same sessions left out
  tail immersion   latency (s) - assays/tail_immersion_all.csv via ti_table_days.table() (the values the user sends)
  Straub tail      % of the first 5 min, automatic score - analysis/days/dayN/straub_dayN.csv (straub_day.py)
  TST              non-moving % - analysis/days/dayN/assays_dayN.csv, tst_pct_nonmoving (automatic, first 180 s)
  hot plate        non-moving % of the first 30 s - assays_dayN.csv, hot_pct_Non-moving (the user's hand score)
TST and hot plate use the non-moving %, the two measures the previous pipeline's dashboard shows
(step06 analyze_passive_active_dashboard1.m: 'TST nonmoving (%)', 'HOT nonmoving (%)').
Outputs: analysis/days/results_days<a>-<b>.png and results_days<a>-<b>.csv (long table: day, cage, mouse, group,
measure, value, source).
usage: python results_plots.py [--days 1 2 ... N]   (default: day 1 up to the last day with a session)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import ti_table_days as TI

T = Path(__file__).resolve().parent
P = T.parent
DAYS_DIR = P / "analysis" / "days"
MASTER = P / "analysis" / "master_day_level.csv"
LEAVE_OUT_GRADES = {"whole_eye", "guess", "no_pupil"}
COLS = ["day", "cage", "mouse", "group", "measure", "value", "source"]
# (measure, title, y label, note under the panel, y range or None)
# 2026-11-23 (user: "put all in our routine ppt with licking PR and pupil all together"): PR and licking added in front
PANELS = [("pr", "PR (last requirement)", "PR = rewards + 1", "active sessions only (passive = yoked, no PR)", None),
          ("licks_per_min", "Licking", "licks / min", "first 15 min 10 s · passive sessions left out (spout in the mouth)",
           None),
          ("pupil_px", "Pupil size", "session median pupil (px)",
           "first 15 min 10 s · whole-eye mice left out (eye damage)", None),
          ("pupil_reward_pct", "Pupil at reward", "pupil change after reward (%)",
           "vs the -1..0 s baseline · whole-eye mice left out", None),
          ("ti_s", "Tail immersion", "latency (s)", "the latencies you send", None),
          ("straub_pct", "Straub tail", "Straub (% of counted time)",
           "first 5 min · automatic v9.1, or your hand score where you scored", (0, 105)),
          ("tst_nonmoving_pct", "Tail suspension (TST)", "non-moving (%)", "automatic, first 180 s", (0, 105)),
          ("hot_nonmoving_pct", "Hot plate", "non-moving (%)", "your hand score, first 30 s", (0, 105))]


def last_day() -> int:
    filed = [int(p.name[3:]) for p in (P / "python").glob("day*") if p.name[3:].isdigit()]
    return max(filed) if filed else 1


def num(v) -> float:
    v = pd.to_numeric(v, errors="coerce")
    return float(v) if pd.notna(v) else float("nan")


def collect(days) -> pd.DataFrame:
    G = {k: v[0] for k, v in TI.groups().items()}
    rows = []

    def add(d, c, m, measure, v, src):
        if np.isfinite(v):
            rows.append(dict(day=int(d), cage=str(c), mouse=str(m), group=G.get((str(c), str(m)), ""),
                             measure=measure, value=float(v), source=src))

    prf = DAYS_DIR / f"pr_table_days{days[0]}-{days[-1]}.csv"   # PR = pr (rewards + 1), as the PR slide (pr_table_days.py)
    if not prf.exists():
        cand = sorted(DAYS_DIR.glob("pr_table_days*.csv"), key=lambda p: p.stat().st_mtime)
        prf = cand[-1] if cand else None
    if prf is not None:
        for _, r in pd.read_csv(prf, dtype={"cage": str}).iterrows():
            if int(r["day"]) in days and str(r.get("task", "")) != "PASSIVE":
                add(r["day"], r["cage"], r["mouse"], "pr", num(r.get("pr")), f"analysis/days/{prf.name}")
    if MASTER.exists():                                   # the tables' pupil (first 15 min 10 s)
        M = pd.read_csv(MASTER, dtype={"cage": str})
        for _, r in M.iterrows():
            if int(r["day"]) not in days:
                continue
            c, m = str(r["cage"]), str(r["mouse"]).split("_", 1)[-1]      # '9669m_black' -> 'black'
            if str(r.get("task", "")) != "PASSIVE":       # passive: the spout is inside the mouth, not licks (11-19)
                add(r["day"], c, m, "licks_per_min", num(r.get("licks_per_min")), "analysis/master_day_level.csv")
            if str(r.get("pupil_grade", "")) in LEAVE_OUT_GRADES:
                continue
            add(r["day"], c, m, "pupil_px", num(r.get("pupil_median_px")), "analysis/master_day_level.csv")
            add(r["day"], c, m, "pupil_reward_pct", num(r.get("pupil_reward_response_pct")),
                "analysis/master_day_level.csv")
    ti_days = [d for d in days if d >= TI.FIRST_TI_DAY]
    if ti_days:
        for _, r in TI.table(ti_days).iterrows():
            add(r["day"], r["cage"], r["mouse"], "ti_s", num(r["latency_s"]), "assays/tail_immersion_all.csv")
    for d in days:
        f = DAYS_DIR / f"day{d}" / f"straub_day{d}.csv"
        if f.exists():
            for _, r in pd.read_csv(f, dtype={"cage": str}).iterrows():
                v = num(r.get("straub_final_pct"))       # 2026-11-23: final = hand (counting frames) or auto v9.1
                add(d, r["cage"], r["mouse"], "straub_pct", v if np.isfinite(v) else num(r.get("auto_straub_pct")),
                    f"analysis/days/day{d}/straub_day{d}.csv")
        f = DAYS_DIR / f"day{d}" / f"assays_day{d}.csv"
        if f.exists():
            for _, r in pd.read_csv(f, dtype={"cage": str}).iterrows():
                add(d, r["cage"], r["mouse"], "tst_nonmoving_pct", num(r.get("tst_pct_nonmoving")),
                    f"analysis/days/day{d}/assays_day{d}.csv")
                add(d, r["cage"], r["mouse"], "hot_nonmoving_pct", num(r.get("hot_pct_Non-moving")),
                    f"analysis/days/day{d}/assays_day{d}.csv")
    return pd.DataFrame(rows, columns=COLS)


def panel(ax, D, key, title, ylabel, note, yr, first):
    g = D[D.measure == key]
    ax.set_title(title, fontsize=12, loc="left", pad=18, color="#1f3a5f", fontweight="bold")
    if g.empty:
        ax.text(0.5, 0.5, "no data yet", transform=ax.transAxes, ha="center", va="center", color="#7b8794",
                fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(note, fontsize=8.5, color="#52606d")
        ax.spines[["top", "right"]].set_visible(False)
        return
    dd = sorted(int(x) for x in g.day.unique())
    TI.draw_plot(ax, g.rename(columns={"value": "latency_s"}), dd)
    ax.set_ylabel(ylabel)
    if yr is not None:
        ax.set_ylim(*yr)
    else:
        lo, hi = float(min(0.0, g.value.min())), float(g.value.max())
        ax.set_ylim(lo - 0.06 * (hi - lo or 1.0), hi + 0.08 * (hi - lo or 1.0))
    n = g.groupby("day").size()
    ax.set_xlabel(f"day · {note} · n {n.min()}-{n.max()} mice per day" if n.min() != n.max()
                  else f"day · {note} · n {n.min()} mice per day", fontsize=8.5, color="#52606d")
    leg = ax.get_legend()
    if leg is not None and not first:
        leg.remove()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+", help="default: day 1 up to the last day with a session")
    a = ap.parse_args(argv)
    days = sorted(a.days) if a.days else list(range(1, last_day() + 1))
    D = collect(days).sort_values(["measure", "day", "cage", "mouse"])
    tag = f"days{days[0]}-{days[-1]}"
    DAYS_DIR.mkdir(parents=True, exist_ok=True)
    out_csv, out_png = DAYS_DIR / f"results_{tag}.csv", DAYS_DIR / f"results_{tag}.png"
    D.to_csv(out_csv, index=False, encoding="utf-8-sig")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 4, figsize=(22, 9.8))
    for k, (ax, spec) in enumerate(zip(axes.ravel(), PANELS)):
        panel(ax, D, *spec, first=(k == 0))
    fig.suptitle(f"Results so far, days {days[0]}-{days[-1]}: PR, licking, pupil, tail immersion, Straub tail, TST, "
                 "hot plate",
                 fontsize=14, color="#1f3a5f", fontweight="bold", x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_png, dpi=115)
    plt.close(fig)
    print(out_png)
    print(out_csv)
    for key, title, *_ in PANELS:
        g = D[D.measure == key]
        if g.empty:
            print(f"  {title}: no data yet")
            continue
        d = int(g.day.max())
        h = g[g.day == d]
        med = h.groupby("group").value.median()
        print(f"  {title}: day {d} " + " · ".join(f"{k or 'no group'} median {v:.3g}" for k, v in med.items())
              + f" (n {len(h)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
