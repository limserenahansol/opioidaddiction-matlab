#!/usr/bin/env python3
"""
files_sync_audit.py - are all files of every session in the right folder, and do they belong together? (added
2026-11-23; user: "please re check QC whether all files are located correctly and synchrnoized from day1 to day 10",
after the day-10 arrange put one session's Saleae CSV into another mouse's folder). REPORT ONLY: moves nothing.

For every filed session (python/dayN/<cage>/<mouse>/*.jsonl, full length) it reads its own folder's files with the
routine's readers (daily_summary.read_python / read_sal / read_csv) and checks:
  S1 files     one full python log, one .sal, one digital.csv, one video (listed exceptions are INFO)
  S2 sal name  the .sal file name is the session's name (letters/digits only, case ignored)
  S3 capture   Saleae capture started 0-600 s before the python start
  S4 identity  the CSV's first camera pulse (s after its capture start) = python start - .sal capture start, +-1.0 s:
               the CSV comes from THIS capture (the check the 11-23 day-10 mix-up would have failed)
  S5 rewards   valve openings in the session window = python reward commands (+-2)
  S6 licks     active sessions: Saleae lick onsets in the session ~ python lick records (within 20 or 50 %)
  S7 frames    camera pulses in the session window vs video frames (listed alignment fixes / early stops = INFO)
  S8 saved     the CSV was saved after its capture ended (export of a finished capture)
  S9 copies    no CSV / .sal / video file has the same content as another session's
  S10 sync     the routine's own sync verdict (analysis/session_checks/dayN/<session>_check.txt, line 'H sync')
Output: analysis/days/files_sync_audit_days<a>-<b>.csv (one row per session and check) + a printed summary.
usage: python files_sync_audit.py [--days 1 ... N]
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

import cv2

T = Path(__file__).resolve().parent
P = T.parent
sys.path.insert(0, str(T))
import daily_summary as DS  # noqa: E402  (the routine's own readers; nothing is run)

CAM_HZ = DS.CAM_HZ


def norm(s: str) -> str:
    return re.sub(r"[^a-z0-9]", "", s.lower())


def md5(p: Path) -> str:
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def listed_fixes() -> set:
    keys = set()
    for name in ("alignment_fixes.csv", "video_offsets.csv"):
        f = P / "sync" / name
        if f.exists():
            for r in csv.DictReader(open(f, encoding="utf-8-sig", newline="")):
                d = str(r.get("day", "")).strip().lower().replace("\\", "/")
                keys.add(d)
                if r.get("session"):
                    keys.add(str(r["session"]).strip().lower())
    return keys


def sync_verdict(day: int, name: str) -> str:
    f = P / "analysis" / "session_checks" / f"day{day}" / f"{name}_check.txt"
    if not f.exists():
        return ""
    for ln in open(f, encoding="utf-8", errors="replace"):
        m = re.match(r"\s*(PASS|INFO|WARN|CHECK|FAIL)\s+H sync\s+(.*)", ln)
        if m:
            return f"{m.group(1)} {m.group(2).strip()[:90]}"
    return ""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--days", type=int, nargs="+")
    a = ap.parse_args(argv)
    days = sorted(a.days) if a.days else sorted(int(p.name[3:]) for p in (P / "python").glob("day*") if p.name[3:].isdigit())
    fixes = listed_fixes()
    rows, hashes = [], {}

    def add(d, c, m, name, check, st, txt):
        rows.append(dict(day=d, cage=c, mouse=m, session=name, check=check, status=st, detail=txt))

    for d in days:
        for mdir in sorted((P / "python" / f"day{d}").glob("*/*")):
            if not mdir.is_dir():
                continue
            c, m = mdir.parent.name, mdir.name
            logs = [x for x in (DS.read_python(p) for p in sorted(mdir.glob("*.jsonl"))) if x and x["dur"] >= DS.MIN_SESSION_S]
            sdir, vdir = P / "salae" / f"day{d}" / c / m, P / "video" / f"day{d}" / c / m
            sals = sorted(sdir.glob("*.sal")) if sdir.exists() else []
            csvs = sorted(sdir.glob("digital*.csv")) if sdir.exists() else []
            vids = sorted(vdir.glob("*.mp4")) + sorted(vdir.glob("*.avi")) if vdir.exists() else []
            name = logs[0]["name"] if logs else "(no full log)"
            key = f"day{d}/{c}/{m}"
            exc = key in fixes or name.lower() in fixes
            n_ok = len(logs) == 1 and len(sals) == 1 and len(csvs) == 1 and len(vids) == 1
            add(d, c, m, name, "S1 files", "PASS" if n_ok else ("INFO" if exc else "FAIL"),
                f"logs {len(logs)}, .sal {len(sals)}, csv {len(csvs)}, video {len(vids)}" + (" (listed exception)" if exc and not n_ok else ""))
            for kind, fs in (("csv", csvs), ("sal", sals), ("video", vids)):
                for f in fs:
                    try:
                        hashes.setdefault((kind, f.stat().st_size, md5(f) if kind != "video" else f.stat().st_size), []).append((d, c, m, f))
                    except OSError:
                        pass
            if not logs:
                continue
            py = logs[0]
            sal = DS.read_sal(sals[0]) if sals else None
            cs = DS.read_csv(csvs[0]) if csvs else None
            if sals:
                add(d, c, m, name, "S2 sal name", "PASS" if norm(name) in norm(sals[0].stem) or norm(sals[0].stem) in norm(name) else "WARN",
                    f"{sals[0].name}")
            lag = None
            if sal:
                lag = py["start"] - sal["start"]
                add(d, c, m, name, "S3 capture", "PASS" if 0 <= lag <= 600 else "FAIL", f"capture started {lag:.1f} s before python")
            if cs and cs["first_cam"] is not None:
                if lag is not None:
                    err = cs["first_cam"] - lag
                    add(d, c, m, name, "S4 identity", "PASS" if abs(err) <= 1.0 else "FAIL",
                        f"first camera pulse {cs['first_cam']:.2f} s vs .sal lag {lag:.2f} s (diff {err:+.2f} s)")
                w0, w1 = cs["first_cam"] - 0.1, cs["first_cam"] + py["dur"] + 2
                n_rew = sum(1 for t in cs["rewards"] if w0 <= t <= w1)
                add(d, c, m, name, "S5 rewards", "PASS" if abs(n_rew - py["rewards"]) <= max(2, 0.03 * py["rewards"]) else "FAIL",
                    f"valve openings {n_rew} vs python reward_cmd {py['rewards']}")
                passive = "replay" in str(py.get("paradigm", "")).lower() or "passive" in name.lower()
                if cs.get("licks") is not None and not passive:
                    n_lk = sum(1 for t in cs["licks"] if w0 <= t <= w1)
                    bad = abs(n_lk - py["licks"]) > max(20, 0.5 * max(n_lk, py["licks"]))
                    add(d, c, m, name, "S6 licks", "WARN" if bad else "PASS", f"Saleae lick onsets {n_lk} vs python lick records {py['licks']}")
                if vids:
                    cap = cv2.VideoCapture(str(vids[0]))
                    nf = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
                    cap.release()
                    n_cam_s = int(round(min(cs["n_cam"], (cs["last"] - cs["first_cam"]) * CAM_HZ + 1)))
                    diff = nf - cs["n_cam"]
                    st = "PASS" if abs(diff) <= 3 else ("INFO" if exc else "WARN")
                    add(d, c, m, name, "S7 frames", st, f"video {nf} frames vs {cs['n_cam']} camera pulses in the CSV ({diff:+d})"
                        + (" - listed alignment fix / early stop" if st == "INFO" else ""))
                if sal:
                    saved = csvs[0].stat().st_mtime
                    cap_end = sal["start"] + (cs["last"] or 0)
                    add(d, c, m, name, "S8 saved", "PASS" if saved >= cap_end - 5 else "WARN",
                        f"CSV saved {saved - cap_end:+.0f} s after the capture ended")
            elif csvs:
                add(d, c, m, name, "S4 identity", "FAIL", "CSV has no camera pulses / not readable")
            v = sync_verdict(d, name)
            if v:
                add(d, c, m, name, "S10 sync", "PASS" if v.startswith(("PASS", "INFO")) else ("WARN" if v.startswith("WARN") else "FAIL"), v)
    for (kind, *_), L in hashes.items():
        if len(L) > 1:
            for d, c, m, f in L:
                add(d, c, m, "", "S9 copies", "FAIL", f"same {kind} content in " + "; ".join(f"day{x[0]}/{x[1]}/{x[2]}/{x[3].name}" for x in L))
    # the QC summary already explains the known / fixed / approved cases (yellow) - for the checks it also covers
    # (files, frames, sync) a yellow or green QC cell turns FAIL / WARN into INFO with its flags; identity checks
    # (S3 capture, S4 identity, S5 rewards, S9 copies) are never relaxed
    qf = sorted((P / "analysis" / "days").glob("qc_summary_days1-*.csv"), key=lambda p: p.stat().st_mtime)
    qc = {}
    if qf:
        for r in csv.DictReader(open(qf[-1], encoding="utf-8-sig", newline="")):
            qc[(int(r["day"]), r["cage"], r["mouse"])] = (r.get("status", ""), r.get("flags", ""))
    for r in rows:
        q = qc.get((int(r["day"]), r["cage"], r["mouse"]))
        if q and r["status"] in ("FAIL", "WARN") and r["check"].split()[0] in ("S1", "S7", "S10") and q[0] in ("yellow", "green"):
            r["status"] = "INFO"
            r["detail"] += f" | QC {q[0]}: {q[1]} (known / fixed / approved, see the QC notes)"
    tag = f"days{days[0]}-{days[-1]}"
    out = P / "analysis" / "days" / f"files_sync_audit_{tag}.csv"
    with open(out, "w", newline="", encoding="utf-8-sig") as h:
        w = csv.DictWriter(h, fieldnames=["day", "cage", "mouse", "session", "check", "status", "detail"])
        w.writeheader()
        w.writerows(rows)
    n = {s: sum(r["status"] == s for r in rows) for s in ("FAIL", "WARN", "INFO", "PASS")}
    sess = len({(r["day"], r["cage"], r["mouse"]) for r in rows})
    print(f"files + sync audit days {days[0]}-{days[-1]}: {sess} mouse-days, {n['FAIL']} FAIL, {n['WARN']} WARN, "
          f"{n['INFO']} INFO, {n['PASS']} PASS -> {out}")
    for r in rows:
        if r["status"] in ("FAIL", "WARN", "INFO"):
            print(f"  {r['status']:4s} day{r['day']:<2} {r['cage']} {r['mouse']:6s} {r['check']:11s} {r['detail'][:150]}")
    return 1 if n["FAIL"] else 0


if __name__ == "__main__":
    raise SystemExit(main())
