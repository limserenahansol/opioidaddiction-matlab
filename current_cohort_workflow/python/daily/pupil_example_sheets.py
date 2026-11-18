#!/usr/bin/env python3
"""
pupil_example_sheets.py - example frames showing the pupil that the analysis tables use (added 2026-11-16).

Why: the user wants to see the tracked pupil of every mouse and session on example frames and correct it ("show me
your new tracked pupil per each mouse per each session as example so I can correct"; then "instead of whole ppt show
me just example of each mouse each session").

Two layouts (the same N frames per session, spread evenly: 1/12, 3/12, ... 11/12 of the session for N = 6):
  * default: one sheet per day, rows = mice  -> analysis/days/dayN/figs_extra/dayN_pupil_examples.png
  * --per-session: one larger image per session (2 rows x 3 frames for N = 6)
                   -> analysis/days/dayN/figs_extra/pupil_examples/dayN_<cage>_<mouse>.png
                   (a subfolder, so the PPT appendix does not pick them up)
Each tile = the eye box (ROI) of that frame, contrast-stretched, with the tracked pupil as a green circle (centre +
equivalent diameter). A best-effort trace that exists but is NOT approved is drawn as an orange dotted circle for
comparison (it is not in any table). Frames where the table has no pupil say "no pupil (gap)".
The table pupil is the one quick_analysis.py uses:
  * grade best_effort (approved in besteffort_sessions.csv) -> pupil_outputs_besteffort/besteffort_CLEAN.csv
  * otherwise the package's pupil_outputs_general/pupil_diameter_FINAL_U-Net_results_CLEAN.csv
Drawing convention (as the green check image of 2026-11-16 that the user called almost correct): crop =
imcrop(frame, roi) = cols x..x+w, rows y..y+h (1-based); CenterX_px is 1-based with pixel centres at integers, so
the circle centre is at (CenterX_px - 0.5, CenterY_px - 0.5) crop-pixel widths from the crop's top-left edge (a
half-pixel more exact than the first sheets of 2026-11-16), radius = diameter / 2; video frame k (1-based) is at
container time (k - 1) / 60 s. Nothing is changed.

usage: python pupil_example_sheets.py --project <project> --days 1 2 3 [--n 6] [--per-session]
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io as sio
from PIL import Image, ImageDraw, ImageFont

sys.path.insert(0, str(Path(__file__).resolve().parent))
import concat_choice as CC  # noqa: E402

PACKAGE_CSV = "pupil_diameter_FINAL_U-Net_results_CLEAN.csv"
TW, TH = 330, 264          # tile box of the day sheet (px)
TWS, THS = 600, 480        # tile box of the per-session image (px)
LAB_W = 250                # row label width (day sheet)
GREEN, ORANGE, GREY = (40, 220, 60), (255, 150, 0), (110, 110, 110)


def font(size, bold=False):
    for f in (("arialbd.ttf" if bold else "arial.ttf"), "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(f, size)
        except OSError:
            continue
    return ImageFont.load_default()


def ffmpeg():
    import imageio_ffmpeg
    return imageio_ffmpeg.get_ffmpeg_exe()


def roi_crop(roi):
    x, y, w, h = [int(round(v)) for v in roi]
    return x - 1, y - 1, w + 1, h + 1


def grab(video: Path, crop, k: int):
    """the eye box of video frame k (1-based). Seeks HALF A FRAME EARLY: the mp4 stamps frames at a 1/60000 time
    base and every 250th step is 999 ticks, so frame k sits up to 224 ticks before (k-1)/60 and a seek to (k-1)/60
    returned frame k+1 (found by the review of 2026-11-16; checked by pixels against an in-order decode).
    crop ...:exact=1 keeps odd crop sizes exact."""
    x0, y0, cw, ch = crop
    cmd = [ffmpeg(), "-v", "error", "-ss", f"{max(0.0, (k - 1.5) / 60.0):.4f}", "-i", str(video), "-frames:v", "1",
           "-vf", f"crop={cw}:{ch}:{x0}:{y0}:exact=1,format=gray", "-f", "rawvideo", "-pix_fmt", "gray", "-"]
    raw = subprocess.run(cmd, capture_output=True).stdout
    return np.frombuffer(raw[: cw * ch], np.uint8).reshape(ch, cw) if len(raw) >= cw * ch else None


def grades(pupil_root: Path):
    q = pupil_root / "pupil_session_quality.csv"
    if not q.exists():
        return {}
    Q = pd.read_csv(q, dtype={"cage": str, "mouse": str}, encoding="utf-8-sig")
    return {(int(r.day), r.cage, r.mouse): str(r.grade) for r in Q.itertuples()}


def pupil_csv(sd: Path, root: Path, grade: str):
    if grade == "best_effort":
        p = CC.besteffort_pupil_csv(sd, root)
        if p is not None:
            return p, "best effort"
    hits = sorted(sd.rglob(PACKAGE_CSV), key=lambda p: p.stat().st_mtime)
    return (hits[-1], "package") if hits else (None, "")


def pending_be(sd: Path, root: Path):
    """best-effort CLEAN of a session that is NOT approved (drawn for comparison only), else None"""
    p = sd / "pupil_outputs_besteffort" / "besteffort_CLEAN.csv"
    return p if p.exists() and CC.besteffort_pupil_csv(sd, root) is None else None


def load_trace(p: Path) -> pd.DataFrame:
    T = pd.read_csv(p)
    dcol = "SmoothedDiameter_px" if "SmoothedDiameter_px" in T.columns else "Diameter_pixels"
    return pd.DataFrame(dict(d=pd.to_numeric(T[dcol], errors="coerce"),
                             cx=pd.to_numeric(T["CenterX_px"], errors="coerce"),
                             cy=pd.to_numeric(T["CenterY_px"], errors="coerce")))


def stretch(img):
    lo, hi = np.percentile(img, [1, 99.5])
    return np.clip((img.astype(float) - lo) * 255.0 / max(1.0, hi - lo), 0, 255).astype(np.uint8)


def _ok(r):
    return r is not None and np.isfinite(r["d"]) and np.isfinite(r["cx"]) and np.isfinite(r["cy"])


def tile(img, row, size=(TW, TH), row2=None, width=2, mark_missing=True):
    """contrast-stretched eye box; green circle = table pupil (row), orange dotted = pending best effort (row2)"""
    tw, th = size
    if img is None:
        t = Image.new("RGB", (tw, th), (235, 235, 235))
        ImageDraw.Draw(t).text((10, th // 2), "frame not read", fill=GREY, font=font(16))
        return t
    s = min(tw / img.shape[1], th / img.shape[0])
    w, h = int(img.shape[1] * s), int(img.shape[0] * s)
    ox, oy = (tw - w) // 2, (th - h) // 2
    t = Image.new("RGB", (tw, th), (0, 0, 0))
    t.paste(Image.fromarray(stretch(img)).convert("RGB").resize((w, h), Image.LANCZOS), (ox, oy))
    dr = ImageDraw.Draw(t)
    # CenterX_px is 1-based with pixel centres at integers; the crop's pixel j covers [j - 1, j] in tile units / s.
    # PIL takes coordinates as pixel centres (-0.5) and draws the outline inward from the box (+width/2 on the radius)
    if _ok(row2):
        cx, cy = ox + (row2["cx"] - 0.5) * s - 0.5, oy + (row2["cy"] - 0.5) * s - 0.5
        r = row2["d"] / 2 * s + (width + 1) / 2
        for a in range(0, 360, 24):
            dr.arc([cx - r, cy - r, cx + r, cy + r], a, a + 12, fill=ORANGE, width=width + 1)
    if _ok(row):
        cx, cy = ox + (row["cx"] - 0.5) * s - 0.5, oy + (row["cy"] - 0.5) * s - 0.5
        r = row["d"] / 2 * s + width / 2
        dr.ellipse([cx - r, cy - r, cx + r, cy + r], outline=GREEN, width=width)
        dr.line([cx - 5, cy, cx + 5, cy], fill=GREEN, width=1)
        dr.line([cx, cy - 5, cx, cy + 5], fill=GREEN, width=1)
    elif mark_missing:
        dr.rectangle([0, 0, 200, 30], fill=(0, 0, 0))
        dr.text((6, 5), "no pupil found", fill=ORANGE, font=font(18, True))
    return t


def collect(project: Path, day: int, n: int, pupil_root: Path = None, compare_root: Path = None):
    """per mouse of the day: what to draw (table pupil, pending best effort, frames) + the frame grabs to do.
    pupil_root: another package output tree (e.g. pupil_tst_teach, the new package of 2026-11-17) - its package
    CLEAN is drawn green; compare_root: a tree whose package CLEAN is drawn orange (e.g. the original package)."""
    root = pupil_root or CC.active_pupil_root(project)       # 2026-11-17: the tree the tables use
    G = grades(root)
    man = pd.read_csv(project / "cohort_manifest.csv", dtype=str) if (project / "cohort_manifest.csv").exists() else None
    order = ([(r.cage_folder, r.mouse) for r in man.sort_values(["cage_folder", "mouse"]).itertuples()] if man is not None
             else [(c, m) for c in ("9669m", "9670f") for m in ("black", "green1", "green2", "white")])
    jobs, rows = [], []
    for cage, mouse in order:
        vdir = project / "video" / f"day{day}" / cage / mouse
        vids = sorted(list(vdir.glob("*.mp4")) + list(vdir.glob("*.avi")), key=lambda p: p.stat().st_size) \
            if vdir.is_dir() else []
        sd = CC.session_dir(root, day, cage, mouse)
        grade = G.get((day, cage, mouse), "-")
        base = dict(day=day, cage=cage, mouse=mouse, grade=grade, src="", tiles=[], frames=[])
        if not vids or sd is None or not (vdir / "roi.mat").exists():
            rows.append(dict(base, note="no video - pupil absent" if not vids else "no pupil output / eye box"))
            continue
        p, src = pupil_csv(sd, root, grade)
        if p is None:
            rows.append(dict(base, note="no pupil trace"))
            continue
        P = load_trace(p)
        if compare_root is not None:
            sd2 = CC.session_dir(compare_root, day, cage, mouse)
            hits = sorted(sd2.rglob(PACKAGE_CSV), key=lambda q: q.stat().st_mtime) if sd2 is not None else []
            P2 = load_trace(hits[-1]) if hits else None
            if P2 is not None and len(P2) != len(P):
                P2 = None
        else:
            pb = pending_be(sd, root)
            P2 = load_trace(pb) if pb is not None and len(pd.read_csv(pb, usecols=[0])) == len(P) else None
        crop = roi_crop(sio.loadmat(vdir / "roi.mat")["roi"].ravel().astype(float))
        ks = [min(max(1, int(round((j + 0.5) / n * len(P))) + 1), len(P)) for j in range(n)]
        rec = dict(base, src=src, note="", ks=ks, P=P, P2=P2, med=float(np.nanmedian(P["d"])),
                   med2=float(np.nanmedian(P2["d"])) if P2 is not None else np.nan,
                   vis=100.0 * float(np.isfinite(P["d"]).mean()), crop=crop, frames=[None] * n)
        rows.append(rec)
        for j, k in enumerate(ks):
            jobs.append((rec, j, vids[-1], crop, k))
    with ThreadPoolExecutor(max_workers=6) as ex:
        imgs = list(ex.map(lambda a: grab(a[2], a[3], a[4]), jobs))
    for (rec, j, _, _, _), img in zip(jobs, imgs):
        rec["frames"][j] = img
    return rows


def day_sheet(project: Path, day: int, n: int, out: Path):
    rows = collect(project, day, n)
    HEAD, CAP, GAP = 92, 46, 10
    W = LAB_W + n * (TW + GAP) + 10
    H = HEAD + len(rows) * (TH + CAP + GAP) + 40
    S = Image.new("RGB", (W, H), "white")
    dr = ImageDraw.Draw(S)
    dr.text((12, 10), f"Day {day} - the pupil the analysis tables use (green circle = tracked pupil; diameter below)",
            fill=(0, 0, 0), font=font(26, True))
    dr.text((12, 48), "Eye box contrast-stretched. If the green circle is not on the black pupil, draw where the pupil is "
                      "(or mark 'no pupil'). Column = time in the session.", fill=(60, 60, 60), font=font(18))
    for i, rec in enumerate(rows):
        y = HEAD + i * (TH + CAP + GAP)
        dr.text((12, y + 10), f"{rec['cage']} {rec['mouse']}", fill=(0, 0, 0), font=font(24, True))
        dr.text((12, y + 44), f"grade: {rec['grade'].replace('_', ' ')}", fill=(40, 40, 40), font=font(18))
        if rec["src"]:
            dr.text((12, y + 70), f"source: {rec['src']}", fill=(40, 40, 40), font=font(18))
            dr.text((12, y + 96), f"median {rec['med']:.0f} px", fill=(40, 40, 40), font=font(18))
            dr.text((12, y + 122), f"pupil in {rec['vis']:.0f} % of frames", fill=(40, 40, 40), font=font(18))
        if not rec["frames"]:
            dr.text((LAB_W + 10, y + TH // 2 - 10), rec["note"], fill=GREY, font=font(22, True))
            continue
        for j, img in enumerate(rec["frames"]):
            x = LAB_W + j * (TW + GAP)
            k = rec["ks"][j]
            S.paste(tile(img, rec["P"].iloc[k - 1]), (x, y))
            d = rec["P"]["d"].iloc[k - 1]
            mm, ss = divmod((k - 1) / 30.26, 60)
            txt = f"{j + 1}  ·  {int(mm):02d}:{int(ss):02d}  ·  " + (f"{d:.0f} px" if np.isfinite(d) else "no pupil")
            dr.text((x + 4, y + TH + 6), txt, fill=(0, 0, 0), font=font(18, True))
            dr.text((x + 4, y + TH + 27), f"frame {k}", fill=GREY, font=font(14))
    dr.text((12, H - 30), "time = frame / 30.26 Hz (real camera rate) · circle = centre and equivalent diameter of the "
                          "tracked blob · tools\\pupil_example_sheets.py", fill=GREY, font=font(15))
    out.parent.mkdir(parents=True, exist_ok=True)
    S.save(out)
    return out


def session_sheets(project: Path, day: int, n: int, out_dir: Path, pupil_root=None, compare_root=None,
                   green_label="the pupil in the tables", orange_label="best effort, NOT approved yet"):
    """one larger image per session (2 rows x 3 frames for n = 6)"""
    outs = []
    for rec in collect(project, day, n, pupil_root, compare_root):
        if not rec["frames"]:
            continue
        cols = 3 if n >= 3 else n
        nr = math.ceil(n / cols)
        HEAD, CAP, GAP = 132, 58, 14
        W = cols * (TWS + GAP) + GAP
        H = HEAD + nr * (THS + CAP + GAP) + 36
        S = Image.new("RGB", (W, H), "white")
        dr = ImageDraw.Draw(S)
        dr.text((GAP, 10), f"Day {day} · {rec['cage']} {rec['mouse']} · grade {rec['grade'].replace('_', ' ')} · "
                           f"table = {rec['src']} (median {rec['med']:.0f} px, pupil in {rec['vis']:.0f} % of frames)",
                fill=(0, 0, 0), font=font(28, True))
        legend = f"green circle = {green_label}"
        if rec["P2"] is not None:
            legend += f"   ·   orange dotted = {orange_label} (median {rec['med2']:.0f} px)"
        dr.text((GAP, 54), legend, fill=(40, 40, 40), font=font(21))
        dr.text((GAP, 90), "If a circle is not on the black pupil, draw where the pupil is (or write 'no pupil').",
                fill=(90, 90, 90), font=font(19))
        for j, img in enumerate(rec["frames"]):
            r_, c_ = divmod(j, cols)
            x, y = GAP + c_ * (TWS + GAP), HEAD + r_ * (THS + CAP + GAP)
            k = rec["ks"][j]
            row2 = rec["P2"].iloc[k - 1] if rec["P2"] is not None else None
            S.paste(tile(img, rec["P"].iloc[k - 1], (TWS, THS), row2, width=3), (x, y))
            d = rec["P"]["d"].iloc[k - 1]
            mm, ss = divmod((k - 1) / 30.26, 60)
            txt = f"{j + 1}  ·  {int(mm):02d}:{int(ss):02d}  ·  " + (f"{d:.0f} px" if np.isfinite(d) else "no pupil")
            if row2 is not None:
                txt += "  ·  best effort " + (f"{row2['d']:.0f} px" if _ok(row2) else "no pupil")
            dr.text((x + 4, y + THS + 6), txt, fill=(0, 0, 0), font=font(22, True))
            dr.text((x + 4, y + THS + 34), f"frame {k}", fill=GREY, font=font(16))
        out = out_dir / f"day{day}_{rec['cage']}_{rec['mouse']}.png"
        out.parent.mkdir(parents=True, exist_ok=True)
        S.save(out)
        outs.append(out)
    return outs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--days", type=int, nargs="+", required=True)
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--per-session", action="store_true")
    ap.add_argument("--pupil-root", help="another package output tree to draw green (e.g. pupil_tst_teach)")
    ap.add_argument("--compare-root", help="a package output tree to draw orange (e.g. pupil_tst_final)")
    ap.add_argument("--out-name", default="pupil_examples", help="sub-folder of figs_extra for --per-session")
    ap.add_argument("--green-label", default="the pupil in the tables")
    ap.add_argument("--orange-label", default="best effort, NOT approved yet")
    a = ap.parse_args(argv)
    project = Path(a.project).resolve()
    for d in a.days:
        figs = project / "analysis" / "days" / f"day{d}" / "figs_extra"
        if a.per_session:
            pr = Path(a.pupil_root).resolve() if a.pupil_root else None
            cr = Path(a.compare_root).resolve() if a.compare_root else None
            for o in session_sheets(project, d, a.n, figs / a.out_name, pr, cr, a.green_label, a.orange_label):
                print(f"  {o}")
        else:
            print(f"  day {d}: {day_sheet(project, d, a.n, figs / f'day{d}_pupil_examples.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
