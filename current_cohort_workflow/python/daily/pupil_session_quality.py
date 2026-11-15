"""Per-session pupil quality for the pupil/TST package output (pupil_tst_final): automatic numbers + manual grade.

Automatic numbers (every session folder with a qc.json and hybrid_detections.csv):
  far     fraction of accepted detections whose centre is more than 0.5 x the session's median diameter away from
          the 1-s running median centre -> the detector jumps to another dark spot (lid, slit, shadow)
  dev30   fraction of accepted detections whose diameter differs > 30 % from the 1-s running median diameter
  held    fraction of frames in which the tracker held its last value (rejected or no pupil)
  no_blob qc.json fracNoBlob
  auto_grade = no_pupil if no_blob > 0.5;  ok if far <= 0.04 and dev30 <= 0.04;  else check (look at qc_sheet.png)
  The thresholds come from the visual check of all 16 day-1/2 sessions (2026-11-15): the 7 sessions that were good
  by eye had far <= 0.026 and dev30 <= 0.028; each of the other 9 had far >= 0.032 or dev30 >= 0.049.

Manual grade (wins over the automatic one): pupil_grades.csv next to the day folders
(columns day, cage, mouse, grade, reason, checked). Grades:
  ok        the pupil is tracked                                                    -> used
  noisy     the detections are the pupil, but with many dips / held frames          -> used, marked
            (lick spout across the eye, blinks, a few wrong blobs)
  guess     the detection is often NOT the pupil (inflamed or narrowed eye, an eye too dark to see the pupil edge)
            -> not used as a measurement; pupil_guess.csv (pupil_guess_trace.py) is the labelled best guess
  no_pupil  no pupil in most frames                                                -> not used

Output: pupil_session_quality.csv in the pupil root (one row per session; use_for_analysis = grade in ok/noisy),
read by quick_analysis.py. Sessions with auto_grade 'check' and no manual grade are printed: look at their
qc_sheet.png and add a row to pupil_grades.csv.

Alignment check (added 2026-11-15): when the video has more or fewer frames than camera pulses, the package's
alignment (pupil_align_ttl, the original Step 01a logic) matches them at the END, while sync_sessions.py decides from
the video file's creation time (normally video frame 0 = first pulse). align_shift_s = camera time the package gives a
video frame minus the time sync_sessions.py gives it (median over frames; 0 when they agree). A non-zero shift means
that in concat_out_*/combined_pupil_digital.csv (read by Step 02) licks and rewards are at the wrong video frames:
a positive shift puts them that much too EARLY on the pupil trace, a negative one too LATE. The daily analysis uses the sync times and is not affected. Nothing is changed here; the session is listed.
The check is made on the NEWEST concat_out (the one prepare_step02_base.py copies for Step 02). When that is a
sync-aligned copy (concat_out_<stamp>_syncaligned, made by fix_concat_alignment.py for a session the user approved in
alignment_fixes.csv), align_check says so and align_shift_package_s keeps the shift of the package's own folder.

usage: python pupil_session_quality.py --pupil-root <pupil_tst_final> [--sync-root <sync>] [--days 1 2]
"""
import argparse, json
from pathlib import Path

import numpy as np
import pandas as pd

USE = {"ok": True, "noisy": True, "guess": False, "no_pupil": False}
FAR_MAX, DEV_MAX, NOBLOB_MAX = 0.04, 0.04, 0.5


def session_numbers(d: Path) -> dict:
    Q = json.loads((d / "qc.json").read_text(encoding="utf-8"))
    out = dict(no_blob=float(Q.get("fracNoBlob", np.nan)), rejected=float(Q.get("fracRejected", np.nan)),
               median_diam_px=float(Q.get("medianDiameter", np.nan)), detector=Q.get("detector", "unet"),
               n_frames=int(Q.get("nFrames", 0)), far=np.nan, dev30=np.nan, held=np.nan)
    r2 = d / "pupil_diameter_FINAL_U-Net_results2.csv"
    if r2.exists():
        R = pd.read_csv(r2)
        col = next((c for c in R.columns if "iam" in c), None)
        if col is not None:
            dia = R[col].to_numpy(float)
            out["held"] = float(np.mean(np.r_[False, np.diff(dia) == 0]))
            out["median_diam_px"] = float(np.nanmedian(dia))
    det = d / "hybrid_detections.csv"
    if det.exists() and np.isfinite(out["median_diam_px"]) and out["median_diam_px"] > 0:
        H = pd.read_csv(det)
        acc = H["accepted_src_1unet_2dark_0none"].to_numpy() > 0
        if acc.sum() >= 30:
            med = out["median_diam_px"]
            cx = pd.Series(np.where(acc, H["best_cx_px"], np.nan))
            cy = pd.Series(np.where(acc, H["best_cy_px"], np.nan))
            bd = pd.Series(np.where(acc, H["best_diam_px"], np.nan))
            k = 61  # ~1 s at the camera rate (about 30 frames/s)
            mx = cx.rolling(k, center=True, min_periods=10).median()
            my = cy.rolling(k, center=True, min_periods=10).median()
            md = bd.rolling(k, center=True, min_periods=10).median()
            dist = np.hypot(cx - mx, cy - my) / med
            ok = dist.notna()
            out["far"] = float(np.mean(dist[ok] > 0.5)) if ok.any() else np.nan
            rel = (bd / md - 1).abs()
            ok = rel.notna()
            out["dev30"] = float(np.mean(rel[ok] > 0.3)) if ok.any() else np.nan
    return out


SYNC_SUFFIX = "_syncaligned"   # folders written by fix_concat_alignment.py


def run_shift(run: Path, frames_csv: Path):
    """(shift_s, n_frames_compared) between one concat_out and sync_sessions.py's frames.csv."""
    if not (run / "combined_pupil_digital.csv").exists():
        return np.nan, 0
    C = pd.read_csv(run / "combined_pupil_digital.csv", usecols=["Frame", "CamTime_s"])
    F = pd.read_csv(frames_csv, usecols=["frame", "t_saleae_s"]).dropna()
    m = C.assign(frame=C["Frame"] - 1).merge(F, on="frame")
    d = (m["CamTime_s"] - m["t_saleae_s"]).dropna()
    return (float(d.median()), int(len(d))) if len(d) else (np.nan, 0)


def usable_runs(session_dir: Path):
    """concat_out folders a reader may use, oldest first: the package's, and sync-aligned copies that passed
    fix_concat_alignment.py's checks (fix_source.json present). prepare_step02_base.py copies the last one."""
    runs = sorted(p for p in session_dir.glob("concat_out_*") if p.is_dir())
    return [p for p in runs if not p.name.endswith(SYNC_SUFFIX) or
            ((p / "fix_source.json").exists() and (p / "combined_pupil_digital.csv").exists())]


def align_shift(session_dir: Path, sync_dir: Path):
    """(shift_s of the newest usable concat_out, shift_s of the package's newest own concat_out, name of the newest if
    it is a sync-aligned copy else '')."""
    runs = usable_runs(session_dir)
    fcs = sorted(sync_dir.glob("*_frames.csv")) if sync_dir.is_dir() else []
    if not runs or not fcs:
        return np.nan, np.nan, ""
    sh, _ = run_shift(runs[-1], fcs[0])
    if not runs[-1].name.endswith(SYNC_SUFFIX):
        return sh, sh, ""
    pkg = [p for p in runs if not p.name.endswith(SYNC_SUFFIX)]
    psh = run_shift(pkg[-1], fcs[0])[0] if pkg else np.nan
    return sh, psh, runs[-1].name


def shift_words(sh: float) -> str:
    """align_shift_s > 0: the table gives each video frame a LATER camera time than the truth, so a lick or reward
    lands on an EARLIER video frame: it sits sh s too early on the pupil trace (and too late when sh < 0)."""
    return f"licks/rewards sit {abs(sh):.2f} s too {'early' if sh > 0 else 'late'} on the pupil trace"


def auto_grade(n: dict):
    if np.isfinite(n["no_blob"]) and n["no_blob"] > NOBLOB_MAX:
        return "no_pupil", f"no pupil in {100 * n['no_blob']:.0f} % of frames"
    if not (np.isfinite(n["far"]) and np.isfinite(n["dev30"])):
        return "check", "no per-frame detections (hybrid_detections.csv): look at qc_sheet.png"
    if n["far"] <= FAR_MAX and n["dev30"] <= DEV_MAX:
        return "ok", f"stable (far {100 * n['far']:.1f} %, dev30 {100 * n['dev30']:.1f} %)"
    return "check", (f"unstable detections (far {100 * n['far']:.1f} %, dev30 {100 * n['dev30']:.1f} %): "
                     "look at qc_sheet.png")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pupil-root", required=True)
    ap.add_argument("--sync-root", help="sync_sessions.py output (default: <pupil-root>/../sync)")
    ap.add_argument("--days", type=int, nargs="*", help="only these days (default: all); other rows are kept")
    ap.add_argument("--report-days", type=int, nargs="*",
                    help="print the grade lines only for these days (alignment problems are printed for every day)")
    a = ap.parse_args(argv)
    root = Path(a.pupil_root)
    sync_root = Path(a.sync_root) if a.sync_root else root.parent / "sync"
    gp = root / "pupil_grades.csv"
    G = pd.read_csv(gp, dtype={"cage": str, "mouse": str}, encoding="utf-8-sig") if gp.exists() else pd.DataFrame()
    outp = root / "pupil_session_quality.csv"
    old = pd.read_csv(outp, dtype={"cage": str, "mouse": str}) if outp.exists() else pd.DataFrame()
    rows = []
    for q in sorted(root.glob("day*/*/*/pupil_outputs_general/qc.json")):
        day_s, cage, mouse = q.parts[-5], q.parts[-4], q.parts[-3]
        if not day_s[3:].isdigit():
            continue
        day = int(day_s[3:])
        if a.days and day not in a.days:
            continue
        n = session_numbers(q.parent)
        ag, why = auto_grade(n)
        sh, psh, fixed = align_shift(q.parent.parent, sync_root / day_s / cage / mouse)
        if not np.isfinite(sh):
            chk = ""
        elif abs(sh) >= 0.02 and fixed:
            chk = (f"MISMATCH: the sync-aligned copy {fixed} no longer matches sync_sessions.py ({shift_words(sh)}; "
                   "frames.csv changed after the fix?): run fix_concat_alignment.py again")
        elif abs(sh) >= 0.02:
            chk = f"MISMATCH: in the package's concat_out {shift_words(sh)}"
        elif fixed:
            pk = f"in the package's own concat_out {shift_words(psh)}" if np.isfinite(psh) else \
                "the package's own concat_out was not found"
            chk = f"ok (sync-aligned copy {fixed}; {pk})"
        else:
            chk = "ok"
        r = dict(day=day, cage=cage, mouse=mouse, **n, auto_grade=ag, auto_reason=why,
                 align_shift_s=sh, align_shift_package_s=psh, align_check=chk)
        m = G[(G["day"].astype(int) == day) & (G["cage"].astype(str) == cage) &
              (G["mouse"].astype(str).str.lower() == mouse.lower())] if len(G) else pd.DataFrame()
        if len(m):
            g = str(m.iloc[-1]["grade"]).strip().lower()
            r.update(grade=g, grade_source="manual (pupil_grades.csv)", reason=str(m.iloc[-1].get("reason", "")))
        else:
            r.update(grade=ag, grade_source="automatic", reason=why)
        r["use_for_analysis"] = USE.get(r["grade"], False)
        rows.append(r)
    T = pd.DataFrame(rows)
    if len(old) and a.days:
        old = old[~old["day"].astype(int).isin(a.days)]
        T = pd.concat([old, T], ignore_index=True).sort_values(["day", "cage", "mouse"])
    cols = ["day", "cage", "mouse", "grade", "use_for_analysis", "grade_source", "reason", "auto_grade", "auto_reason",
            "no_blob", "rejected", "held", "far", "dev30", "median_diam_px", "n_frames", "detector",
            "align_shift_s", "align_shift_package_s", "align_check"]
    T = T[[c for c in cols if c in T.columns]]
    T.to_csv(outp, index=False, float_format="%.4f")
    print(f"pupil quality -> {outp}")
    for _, r in T.iterrows():
        show = a.report_days if a.report_days else a.days
        if show and int(r["day"]) not in show:
            continue
        flag = "  <- look at qc_sheet.png, then add a row to pupil_grades.csv" \
            if (r["grade_source"] == "automatic" and r["grade"] == "check") else ""
        print(f"  day{int(r['day'])} {r['cage']} {r['mouse']:7s} {r['grade']:8s} ({r['grade_source']}): {r['reason']}{flag}")
    mis = T[T["align_check"].astype(str).str.startswith("MISMATCH")] if "align_check" in T else T.iloc[:0]
    if a.days:   # (with --report-days only, every day's alignment problems are still listed)
        mis = mis[mis["day"].astype(int).isin(a.days)]
    approved = set()
    apf = root / "alignment_fixes.csv"
    if apf.exists():
        A = pd.read_csv(apf, dtype=str, encoding="utf-8-sig").fillna("")
        approved = {(int(x.day), x.cage.strip(), x.mouse.strip()) for x in A.itertuples() if str(x.day).strip()}
    for _, r in mis.iterrows():
        k = (int(r["day"]), str(r["cage"]), str(r["mouse"]))
        if k in approved:
            advice = "it is approved in alignment_fixes.csv: run fix_concat_alignment.py, then prepare_step02_base.py"
        else:
            advice = ("the package matched video and camera pulses at the end, sync_sessions.py at the start; Step 02 "
                      "reads this table. Ask the user; if approved, add the session to alignment_fixes.csv and run "
                      "fix_concat_alignment.py")
        print(f"  !! day{int(r['day'])} {r['cage']} {r['mouse']}: {r['align_check']} ({advice})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
