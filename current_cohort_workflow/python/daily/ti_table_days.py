#!/usr/bin/env python3
"""
ti_table_days.py - the tail-immersion latency of every mouse on every day so far (added 2026-11-20; user: "in summary
we should also show TI data aswell ofc in all days"). Report only: reads assays/tail_immersion_all.csv (written by
tail_immersion_entry.py from the latencies the user sends in chat) and changes nothing.

  latency  the value the user sent for that mouse-day, in s (the original Step 02 reads the same value: cell B1 of
           assays/dayN/<cage>/<mouse>/<cage><mouse>_immersion_dayN.csv)
  change   latency minus the mouse's day-5 latency (day 5 = the first TI day, the last natural-reward day before
           morphine)
  phases   Baseline 1-5 (natural reward) · Yoked 6-10 (morphine; Active earns, Passive yoked) · Equal Access 11-13
           (morphine) · Withdrawal 14-16 (natural reward) · Re-exposure 17-18 (morphine)
  empty    "not sent" = the mouse had a session that day but no TI value yet; "-" = no session that day (yet)
Outputs: analysis/days/ti_table_days<a>-<b>.csv (one row per mouse and day), ti_table_days<a>-<b>.png (table + plot)
and ti_plot_days<a>-<b>.png (the plot alone, for the routine deck). Active red, Passive blue, as in cohort_plots.py.

usage: python ti_table_days.py [--days 5 6 7]   (default: day 5 up to the last day with a TI value or a session)
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

T = Path(__file__).resolve().parent
P = T.parent
DAYS_DIR = P / "analysis" / "days"
TI_ALL = P / "assays" / "tail_immersion_all.csv"
MICE = [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")]
FIRST_TI_DAY = 5                                   # tail immersion every day from day 5 (user 2026-11-18)
GCOL = {"Active": "#c0392b", "Passive": "#2c6fbb"}  # as cohort_plots.py
PHASES = [(1, 5, "Baseline", "natural reward"), (6, 10, "Yoked", "morphine"), (11, 13, "Equal Access", "morphine"),
          (14, 16, "Withdrawal", "natural reward"), (17, 18, "Re-exposure", "morphine")]


def phase(d):
    for a, b, name, rew in PHASES:
        if a <= d <= b:
            return name, rew
    return "", ""


def groups():
    f = P / "cohort_manifest.csv"
    if not f.exists():
        return {}
    G = pd.read_csv(f, dtype=str).fillna("")
    return {(r.cage_folder, r.mouse): (r.group, r.pair_id) for r in G.itertuples()}


def has_session(d, c, m):
    f = P / "python" / f"day{d}" / c / m
    return f.is_dir() and any(f.glob("*.jsonl"))


def read_all():
    if not TI_ALL.exists():
        return pd.DataFrame(columns=["day", "cage", "mouse", "latency_s", "entered", "source"])
    return pd.read_csv(TI_ALL, dtype={"cage": str, "mouse": str})


def all_days(A=None):
    A = read_all() if A is None else A
    last_ti = int(A.day.max()) if len(A) else FIRST_TI_DAY
    filed = [int(p.name[3:]) for p in (P / "python").glob("day*") if p.name[3:].isdigit()]
    return list(range(FIRST_TI_DAY, max([last_ti] + filed) + 1))


def table(days=None):
    """one row per mouse and day"""
    A = read_all()
    days = days or all_days(A)
    G = groups()
    rows = []
    for d in days:
        ph, rew = phase(d)
        for c, m in MICE:
            grp, pair = G.get((c, m), ("", ""))
            r = A[(A.day == d) & (A.cage == c) & (A.mouse == m)]
            rows.append(dict(day=d, phase=ph, reward=rew, cage=c, mouse=m, group=grp, pair_id=pair,
                             latency_s=float(r.latency_s.iloc[-1]) if len(r) else float("nan"),
                             session=has_session(d, c, m),
                             entered=str(r.entered.iloc[-1]) if len(r) else "",
                             source=str(r.source.iloc[-1]) if len(r) else ""))
    Tb = pd.DataFrame(rows)
    base = Tb[Tb.day == FIRST_TI_DAY].set_index(["cage", "mouse"]).latency_s
    Tb["day5_s"] = [base.get((c, m), float("nan")) for c, m in zip(Tb.cage, Tb.mouse)]
    Tb["change_vs_day5_s"] = Tb.latency_s - Tb.day5_s
    Tb["ratio_vs_day5"] = Tb.latency_s / Tb.day5_s
    return Tb


def gap_xy(x, y, break_gaps=True):
    """x, y with a NaN between days more than 1 apart, so a line does not run across days without a value"""
    x, y = list(x), list(y)
    if not break_gaps or len(x) < 2:
        return x, y
    xs, ys = [x[0]], [y[0]]
    for a, b, v in zip(x[:-1], x[1:], y[1:]):
        if b - a > 1:
            xs.append((a + b) / 2); ys.append(float("nan"))
        xs.append(b); ys.append(v)
    return xs, ys


def draw_plot(ax, Tb, days, break_gaps=False):
    """break_gaps (2026-11-24 review): no line across missing days (results PR / licks: Passive Yoked days)"""
    import matplotlib.pyplot as plt  # noqa: F401  (the caller set the backend)
    drawn = 0
    for a, b, name, rew in PHASES:
        lo, hi = max(a, days[0]) - 0.5, min(b, days[-1]) + 0.5
        if lo >= hi:
            continue
        if rew == "morphine":
            ax.axvspan(lo, hi, color="#fdf0e2", zorder=0)
        # 2026-11-28 (day-14 check D06): a phase that shows only 1-2 days AFTER another phase (e.g. Withdrawal on day
        # 14) is written a line higher, so it does not run into the label before it (the first phase stays down,
        # away from the plot title)
        lift = hi - lo < 2.6 and drawn > 0
        ax.text((lo + hi) / 2, 1.075 if lift else 1.01, name, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=8, color="#7b8794")
        drawn += 1
    ends = []                                       # (last day, last latency, label, colour) of every mouse
    for c, m in MICE:
        g = Tb[(Tb.cage == c) & (Tb.mouse == m)].sort_values("day")
        g = g[g.latency_s.notna()]
        if g.empty:
            continue
        col = GCOL.get(g.group.iloc[0], "#7b8794")
        ax.plot(*gap_xy(g.day, g.latency_s, break_gaps), "-o", color=col, alpha=0.4, lw=1.1, ms=3.5, zorder=2)
        ends.append([int(g.day.iloc[-1]), float(g.latency_s.iloc[-1]), f"{c[-1]}{m}", col])
    # 2026-11-28 (day-14 check D06): the gap is a share of the whole y range (all days), not of today's values - on
    # day 14 all 8 latencies were 0.8-1.44 s and 6 of the 8 labels sat on top of each other
    yv = pd.to_numeric(Tb.latency_s, errors="coerce")
    gap = 0.052 * max(float(yv.max()) if yv.notna().any() else 1.0, 1.0)
    for d in sorted({e[0] for e in ends}):
        L = sorted((e for e in ends if e[0] == d), key=lambda e: e[1])
        y0 = [e[1] for e in L]
        for k in range(1, len(L)):
            L[k][1] = max(L[k][1], L[k - 1][1] + gap)
        shift = (sum(e[1] for e in L) - sum(y0)) / len(L)          # centre the stack on the points ...
        shift = min(shift, max(0.0, L[0][1] - gap / 2))              # ... but keep the lowest label above 0 s
        for e in L:
            e[1] -= shift
        for (x, y, lab, col), yp in zip(L, y0):
            if abs(y - yp) > gap / 3:                                # moved label: a thin line to its point
                ax.plot([x + 0.02, x + 0.07], [yp, y], color=col, lw=0.5, alpha=0.6, zorder=4)
            ax.text(x + 0.08, y, lab, fontsize=6.5, color=col, va="center", zorder=4,
                    bbox=dict(boxstyle="square,pad=0.05", facecolor="white", edgecolor="none", alpha=0.75))
    for grp, col in GCOL.items():
        M = Tb[(Tb.group == grp) & Tb.latency_s.notna()].groupby("day").latency_s.median()
        if len(M):
            ax.plot(*gap_xy(M.index, M.values, break_gaps), "-o", color=col, lw=3, ms=6, zorder=3,
                    label=f"{grp} median")
    ax.plot([], [], "-o", color="#7b8794", alpha=0.5, lw=1.1, ms=3.5, label="one mouse (f = 9670f, m = 9669m)")
    ax.set_xticks(days)
    ax.set_xlim(days[0] - 0.4, days[-1] + 0.6)
    ax.set_xlabel("day")
    ax.set_ylabel("tail-immersion latency (s)")
    ax.set_ylim(bottom=0)
    ax.grid(axis="y", color="#e4e7eb", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="upper left")


def figure(Tb, days, out_png: Path, out_plot: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    nr, nc = len(MICE) + 3, len(days) + 1            # header + mice + Active / Passive median
    tw = 1.55 * nc
    fig = plt.figure(figsize=(tw + 6.2, 0.62 * (nr + 1.6)), dpi=130)
    ax = fig.add_axes([0.005, 0.02, tw / (tw + 6.2) - 0.01, 0.96])
    ax.set_xlim(-0.05, nc + 0.05)
    ax.set_ylim(nr + 0.75, -0.85)
    ax.axis("off")
    v = Tb.latency_s.dropna()
    lo, hi = (v.min(), v.max()) if len(v) else (0, 1)
    cmap = plt.get_cmap("Oranges")

    def cell(x, y, txt, sub="", face="white", bold=False, col="#1f2933", size=None):
        ax.add_patch(plt.Rectangle((x, y), 1, 1, facecolor=face, edgecolor="#c9ced6", lw=0.8))
        ax.text(x + 0.5, y + (0.42 if sub else 0.52), txt, ha="center", va="center",
                fontsize=size or (12 if sub else 11),
                fontweight="bold" if bold else "normal", color=col)
        if sub:
            ax.text(x + 0.5, y + 0.78, sub, ha="center", va="center", fontsize=7.5, color="#52606d")

    ax.text(0.0, -0.45, "Tail immersion per mouse and day  ·  latency (s)  ·  small line = change from day 5 (before "
            "morphine)", fontsize=11, fontweight="bold", color="#1f3a5f", va="center")
    cell(0, 0, "mouse (group)", face="#1f3a5f", bold=True, col="white")
    for j, d in enumerate(days):
        cell(j + 1, 0, f"day {d}\n{phase(d)[0]}", face="#1f3a5f", bold=True, col="white", size=9.5)   # 2 lines (11-24)
    for i, (c, m) in enumerate(MICE, start=1):
        g = Tb[(Tb.cage == c) & (Tb.mouse == m)]
        grp = g.group.iloc[0] if len(g) else ""
        cell(0, i, f"{c} {m}", sub=grp, bold=True, col=GCOL.get(grp, "#1f2933"))
        for j, d in enumerate(days):
            r = g[g.day == d]
            if r.empty or pd.isna(r.latency_s.iloc[0]):
                sent = bool(len(r)) and bool(r.session.iloc[0])
                cell(j + 1, i, "not sent" if sent else "-", face="#fff6db" if sent else "#f2f2f2",
                     col="#7a5a00" if sent else "#1f2933")
                continue
            r = r.iloc[0]
            sub = (f"{r.change_vs_day5_s:+.2f} s vs day {FIRST_TI_DAY}" if d > FIRST_TI_DAY and
                   pd.notna(r.change_vs_day5_s) else ("before morphine" if d == FIRST_TI_DAY else ""))
            face = matplotlib.colors.to_hex(cmap(0.08 + 0.45 * (r.latency_s - lo) / max(1e-6, hi - lo)))
            cell(j + 1, i, f"{r.latency_s:.2f}", sub=sub, face=face, bold=True)
    for k, grp in enumerate(("Active", "Passive")):
        i = len(MICE) + 1 + k
        cell(0, i, f"{grp} median", sub="(range)", bold=True, face="#fff6db", col=GCOL[grp])
        for j, d in enumerate(days):
            x = Tb[(Tb.day == d) & (Tb.group == grp)].latency_s.dropna()
            if len(x):
                cell(j + 1, i, f"{x.median():.2f}", sub=f"{x.min():.2f}-{x.max():.2f} · n {len(x)}", face="#fff6db",
                     bold=True)
            else:
                cell(j + 1, i, "-", face="#fff6db")
    ax.text(0.0, nr + 0.4, "Source: assays\\tail_immersion_all.csv (the values you sent; Step 02 reads the same, B1 of "
            "assays\\dayN\\<cage>\\<mouse>\\*_immersion_dayN.csv) · tools\\ti_table_days.py",
            fontsize=7.5, color="#7b8794", va="center")
    axp = fig.add_axes([tw / (tw + 6.2) + 0.045, 0.14, 6.2 / (tw + 6.2) - 0.075, 0.74])
    draw_plot(axp, Tb, days)
    fig.savefig(out_png, dpi=130)
    plt.close(fig)
    fig, axp = plt.subplots(figsize=(6.4, 4.6), dpi=150)
    draw_plot(axp, Tb, days)
    axp.set_title("Tail immersion per day", fontsize=11, loc="left", pad=16, color="#1f3a5f", fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_plot, dpi=150)
    plt.close(fig)


def tag_of(days):
    return f"days{days[0]}-{days[-1]}" if len(days) > 1 else f"day{days[0]}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+")
    a = ap.parse_args(argv)
    Tb = table(a.days)
    if Tb.empty or Tb.latency_s.isna().all():
        print("TI table: no tail-immersion value yet (assays/tail_immersion_all.csv)")
        return 0
    days = sorted(int(d) for d in Tb.day.unique())
    tag = tag_of(days)
    out_csv, out_png, out_plot = (DAYS_DIR / f"ti_table_{tag}.csv", DAYS_DIR / f"ti_table_{tag}.png",
                                  DAYS_DIR / f"ti_plot_{tag}.png")
    Tb.to_csv(out_csv, index=False, encoding="utf-8-sig")
    figure(Tb, days, out_png, out_plot)
    print(f"TI table -> {out_csv}\n            {out_png}\n            {out_plot}")
    for d in days:
        parts = []
        for grp in ("Active", "Passive"):
            x = Tb[(Tb.day == d) & (Tb.group == grp)].latency_s.dropna()
            if len(x):
                parts.append(f"{grp} median {x.median():.2f} s (n {len(x)})")
        if parts:
            print(f"  day {d} ({phase(d)[0]}): " + ", ".join(parts))
    ns = Tb[Tb.session & Tb.latency_s.isna()]
    if len(ns):
        print("  not sent yet (session but no TI): " + ", ".join(f"d{r.day} {r.cage} {r.mouse}" for r in ns.itertuples()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
