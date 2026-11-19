#!/usr/bin/env python3
"""
pupil_dark_all.py - run pupil_dark_contrast.py (package pupil repaired where it jumps to another dark blob; strong
dark contrast + glint) on every pupil session of the given days and make one summary (added 2026-11-18; the user:
"try this for day1,2,3,4,5 all"). Nothing is put into the tables here: the summary says per session whether the
repair is usable and how much it changes; the user decides which sessions get it.

Sessions = <pupil tree>/day<N>/<cage>/<mouse> with the package's outputs (results2, hybrid_detections, CLEAN) and a
video. Mice in wholeeye_mice.csv (eye damage: whole eye instead of the pupil) are listed but not run.
Jobs run as separate processes at below-normal priority, --workers at a time; no job starts while the task
(main_hansol.py) is running (user: "실험 중엔 멈춤") - the runner waits.
Output: analysis/pupil_dark/summary_days<a>-<b>.csv (+ .png), run_log.csv; per session the files of
pupil_dark_contrast.py in <session>/pupil_outputs_darkcontrast/.
usage: python pupil_dark_all.py --days 1-5 [--workers 6] [--redo | --rebuild] [--only day5_9670f_green2 ...]
"""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import concat_choice as CC  # noqa: E402

GEN = "pupil_outputs_general"
OUT = P / "analysis" / "pupil_dark"


def task_running() -> bool:
    """True while the behaviour task (main_hansol.py) runs (same test as the jobs guard: its command line)"""
    try:
        r = subprocess.run(["powershell", "-NoProfile", "-Command",
                            "Get-CimInstance Win32_Process -Filter \"Name='python.exe' or Name='pythonw.exe'\" | "
                            "Where-Object { $_.CommandLine -like '*main_hansol*' } | Measure-Object | "
                            "Select-Object -ExpandProperty Count"], capture_output=True, text=True, timeout=60)
        return int((r.stdout or "0").strip() or 0) > 0
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return False


def sessions(root: Path, days):
    whole = {(str(r["cage"]).lower(), str(r["mouse"]).lower())
             for r in CC.read_rows(root / CC.WHOLE_MICE, None, need=("cage", "mouse")) if CC.is_approved(r.get("approved"))}
    out = []
    for d in days:
        for sd in sorted((root / f"day{d}").glob("*/*")):
            if not sd.is_dir():
                continue
            cage, mouse = sd.parent.name, sd.name
            g = sd / GEN
            vid = sorted((P / "video" / f"day{d}" / cage / mouse).glob("*.mp4"))
            note = ""
            if (cage.lower(), mouse.lower()) in whole:
                note = "whole eye (eye damage) - pupil not used"
            elif not vid:
                note = "no video"
            elif not all((g / f).exists() for f in ("pupil_diameter_FINAL_U-Net_results2.csv", "hybrid_detections.csv",
                                                    "pupil_diameter_FINAL_U-Net_results_CLEAN.csv")):
                note = "no package output"
            out.append(dict(day=d, cage=cage, mouse=mouse, sd=sd, note=note, videos=len(vid)))
    return out


def run_one(s, redo: bool, rebuild=False):
    od = s["sd"] / "pupil_outputs_darkcontrast"
    if rebuild:
        if not (od / "pass1.json").exists():
            return s, "no pass 1 yet (run without --rebuild)", 0.0
    elif not redo and (od / "source.json").exists() and (od / "pupil_dark_glint.csv").exists():
        return s, "kept (done before)", 0.0
    while task_running():
        time.sleep(60)
    t0 = time.time()
    cmd = [sys.executable, "-B", str(T / "pupil_dark_contrast.py"), "--day", str(s["day"]), "--cage", s["cage"],
           "--mouse", s["mouse"], "--quiet"] + (["--rebuild", "--no-sheet"] if rebuild else ["--all"])
    flags = getattr(subprocess, "BELOW_NORMAL_PRIORITY_CLASS", 0)
    r = subprocess.run(cmd, capture_output=True, text=True, creationflags=flags)
    msg = "ok" if r.returncode == 0 else "FAILED: " + (r.stderr or r.stdout).strip().splitlines()[-1][:300]
    return s, msg, time.time() - t0


def summary(rows, days, out_csv: Path):
    keys = ["day", "cage", "mouse", "status", "usable", "not_usable_because", "glint_pct", "dark_found_pct",
            "same_blob_pct_of_both", "same_blob_pct_worst_minute", "replaced_pct", "size_gate_none_pct", "package_unchecked_pct", "size_factor",
            "package_median_px", "repaired_median_px", "frames_changed_gt_1px_pct", "port_check_within_0p01_pct",
            "port_check_p99_px", "folder"]
    table = []
    for s in rows:
        sj = s["sd"] / "pupil_outputs_darkcontrast" / "source.json"
        r = dict(day=s["day"], cage=s["cage"], mouse=s["mouse"], status=s.get("status", s["note"]),
                 folder=str(s["sd"] / "pupil_outputs_darkcontrast"))
        if not s["note"] and sj.exists():
            i = json.loads(sj.read_text(encoding="utf-8"))
            a = i.get("actions", {})
            r.update(usable=i.get("usable"), not_usable_because="; ".join(i.get("not_usable_because", [])),
                     glint_pct=i.get("glint_pct"), dark_found_pct=i.get("dark_found_pct"),
                     same_blob_pct_of_both=i.get("same_blob_pct_of_both"),
                     same_blob_pct_worst_minute=i.get("same_blob_pct_worst_minute"), replaced_pct=a.get("replaced_pct"),
                     size_gate_none_pct=a.get("other_blob_size_gate_none_pct"),
                     package_unchecked_pct=a.get("package_unchecked_pct"),
                     size_factor=i.get("size_factor_package_over_dark"),
                     package_median_px=i["median_diam_px"]["package"], repaired_median_px=i["median_diam_px"]["repaired"],
                     frames_changed_gt_1px_pct=i.get("frames_changed_gt_1px_pct"),
                     port_check_within_0p01_pct=i["port_check"].get("frames_diam_within_0p01px_pct"),
                     port_check_p99_px=i["port_check"].get("diam_p99_abs_diff_px"))
        table.append(r)
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(table)
    return table


def figure(table, png: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    run = [r for r in table if r.get("usable") is not None]
    if not run:
        return
    mice = sorted({(r["cage"], r["mouse"]) for r in run})
    days = sorted({r["day"] for r in run})
    fig, ax = plt.subplots(1, 2, figsize=(14, 0.55 * len(mice) + 2.4))
    M = np.full((len(mice), len(days)), np.nan)
    U = np.zeros_like(M, dtype=bool)
    for r in run:
        i, j = mice.index((r["cage"], r["mouse"])), days.index(r["day"])
        U[i, j] = str(r["usable"]) == "True"
        M[i, j] = r["frames_changed_gt_1px_pct"] if U[i, j] else np.nan
    ax[0].imshow(np.where(np.isnan(M), 0.5, np.nan), cmap="Greys", vmin=0, vmax=3, aspect="auto")
    im = ax[0].imshow(M, cmap="Oranges", vmin=0, vmax=max(10, np.nanmax(M) if np.isfinite(M).any() else 10),
                      aspect="auto")
    for i in range(len(mice)):
        for j in range(len(days)):
            ax[0].text(j, i, f"{M[i, j]:.1f}%\nrepair OK" if U[i, j] else "package kept\n(repair not usable)",
                       ha="center", va="center", fontsize=7.5, color="k")
    ax[0].set_xticks(range(len(days)), [f"day {d}" for d in days])
    ax[0].set_yticks(range(len(mice)), [f"{c} {m}" for c, m in mice])
    ax[0].set_title("usable sessions: % frames where the repaired pupil differs from the package by > 1 px",
                    fontsize=10)
    fig.colorbar(im, ax=ax[0], fraction=0.03)
    for (c, m) in mice:
        rr = sorted([r for r in run if (r["cage"], r["mouse"]) == (c, m)], key=lambda r: r["day"])
        x = [r["day"] for r in rr]
        l, = ax[1].plot(x, [r["package_median_px"] for r in rr], "o--", ms=4, lw=1, alpha=0.7)
        ax[1].plot(x, [r["repaired_median_px"] if str(r["usable"]) == "True" else r["package_median_px"] for r in rr],
                   "o-", ms=4, lw=1.6, color=l.get_color(), label=f"{c} {m}")
    ax[1].set_xlabel("day")
    ax[1].set_ylabel("session median pupil diameter (px)")
    ax[1].set_title("dashed = package, solid = with the repair where usable", fontsize=10)
    ax[1].legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", default="1-5")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--redo", action="store_true")
    ap.add_argument("--rebuild", action="store_true", help="pass 2 only (from the saved pass 1), no check sheets")
    ap.add_argument("--only", nargs="*", default=None, help="session names day<N>_<cage>_<mouse>")
    a = ap.parse_args(argv)
    lo, _, hi = a.days.partition("-")
    days = list(range(int(lo), int(hi or lo) + 1))
    root = CC.active_pupil_root(P)
    OUT.mkdir(parents=True, exist_ok=True)
    ss = sessions(root, days)
    todo = [s for s in ss if not s["note"] and (a.only is None or f"day{s['day']}_{s['cage']}_{s['mouse']}" in a.only)]
    for s in ss:
        if s["note"]:
            s["status"] = s["note"]
    print(f"{len(todo)} sessions to run ({len(ss) - len(todo)} not run: "
          f"{sorted({s['note'] for s in ss if s['note']})}), {a.workers} at a time", flush=True)
    log = []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_one, s, a.redo, a.rebuild) for s in todo]
        for f in as_completed(futs):
            s, msg, dt = f.result()
            s["status"] = msg
            print(f"  day{s['day']} {s['cage']} {s['mouse']}: {msg} ({dt / 60:.1f} min)", flush=True)
            log.append([datetime.now().isoformat(timespec="seconds"), s["day"], s["cage"], s["mouse"], msg,
                        round(dt, 1)])
    lp = OUT / "run_log.csv"
    new = not lp.exists()
    with open(lp, "a", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["when", "day", "cage", "mouse", "result", "seconds"])
        w.writerows(log)
    tag = f"days{days[0]}-{days[-1]}"
    table = summary(ss, days, OUT / f"summary_{tag}.csv")
    figure(table, OUT / f"summary_{tag}.png")
    print(f"summary -> {OUT / f'summary_{tag}.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
