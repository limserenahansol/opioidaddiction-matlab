"""Daily eye-box (ROI) check for the pupil tracker: one PNG per day, one row per mouse.

Mice sit a little differently every day, and the package re-uses the previous day's ROI for a new day. This sheet shows
3 frames per session (25 / 50 / 75 % of the video), contrast stretched on the dark range so the eye and pupil are
visible, with the ROI box (yellow). Each row title gives the pupil tracker's result from qc.json (no pupil = % of frames
without an accepted pupil; red when > 50 %).

If an eye is not inside its box: redraw that ROI (RUN_ALL_18DAYS.bat with "redraw_roi": ["dayN cage_mouse"]) and the
session is redone. A narrowed or closed eye (e.g. eye inflammation) cannot be fixed with the box; it shows as "no pupil".

usage: python pupil_roi_check.py --root <project> --day N [--pupil-root <pupil_tst_final>] [--out <png>]
"""
import argparse, json, subprocess
from pathlib import Path

import imageio_ffmpeg
import numpy as np
import scipy.io as sio
from PIL import Image, ImageDraw

W, H = 1920, 1080


def grab(ff, video, t):
    raw = subprocess.run([ff, "-hide_banner", "-loglevel", "error", "-ss", f"{t:.3f}", "-i", str(video), "-frames:v", "1",
                          "-f", "rawvideo", "-pix_fmt", "gray", "-"], capture_output=True, check=True).stdout
    return np.frombuffer(raw[: W * H], np.uint8).reshape(H, W) if len(raw) >= W * H else None


def duration(ff, video):
    r = subprocess.run([ff, "-hide_banner", "-i", str(video)], capture_output=True, text=True).stderr
    try:
        hh, mm, ss = r.split("Duration: ")[1].split(",")[0].split(":")
        return int(hh) * 3600 + int(mm) * 60 + float(ss)
    except (IndexError, ValueError):
        return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True); ap.add_argument("--day", type=int, required=True)
    ap.add_argument("--pupil-root"); ap.add_argument("--out")
    a = ap.parse_args(argv)
    root = Path(a.root); pupil_root = Path(a.pupil_root) if a.pupil_root else root / "pupil_tst_final"
    out = Path(a.out) if a.out else root / "analysis" / "days" / f"day{a.day}" / "figs" / "01_pupil_roi_check.png"
    vids = sorted((root / "video" / f"day{a.day}").glob("*/*/*.mp4")) + sorted((root / "video" / f"day{a.day}").glob("*/*/*.avi"))
    if not vids:
        print(f"(ROI check: no videos in video\\day{a.day})"); return 0
    ff = imageio_ffmpeg.get_ffmpeg_exe()
    grades = {}
    pq = pupil_root / "pupil_session_quality.csv"
    if pq.exists():
        import csv
        with open(pq, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                if row.get("day", "").strip() == str(a.day):
                    grades[(row["cage"], row["mouse"].lower())] = (row.get("grade", ""), row.get("reason", ""))
    TW, TH, LAB = 400, 300, 26
    sheet = Image.new("RGB", (3 * TW + 20, len(vids) * (TH + LAB + 8) + 50), "white")
    dr = ImageDraw.Draw(sheet)
    dr.text((8, 6), f"day {a.day}: eye box (ROI) check - every eye should be inside its yellow box", fill=(0, 0, 0))
    dr.text((8, 22), "(dark range stretched; a narrowed / closed eye cannot be fixed with the box)", fill=(90, 90, 90))
    for r, v in enumerate(vids):
        cage, mouse = v.parts[-3], v.parts[-2]
        rf = v.parent / "roi.mat"
        roi = sio.loadmat(rf)["roi"].ravel().astype(float) if rf.exists() else None
        q = pupil_root / f"day{a.day}" / cage / mouse / "pupil_outputs_general" / "qc.json"
        note, bad = "pupil not run yet", False
        if q.exists():
            try:
                j = json.loads(q.read_text(encoding="utf-8")); nb = 100 * float(j.get("fracNoBlob", np.nan))
                note = f"no pupil {nb:.0f} % of frames ({j.get('detector', 'unet')})"; bad = nb > 50
            except (ValueError, OSError):
                note = "qc.json unreadable"
        # session grade from pupil_session_quality.csv (2026-11-15): ok / noisy / guess / check / no_pupil
        g = grades.get((cage, mouse.lower()))
        if g:
            note += f"   |   grade {g[0]}" + (" (look at qc_sheet.png)" if g[0] == "check" else "")
            bad = bad or g[0] in ("guess", "no_pupil", "check")
        y_row = 50 + r * (TH + LAB + 8)
        dr.text((8, y_row + 6), f"{cage} {mouse}   ROI {'none yet' if roi is None else ' '.join(str(int(v_)) for v_ in roi)}   |   {note}",
                fill=(200, 0, 0) if bad else (0, 0, 0))
        d = duration(ff, v)
        for c, frac in enumerate((0.25, 0.5, 0.75)):
            F = grab(ff, v, frac * d) if d else None
            if F is None:
                continue
            if roi is not None:
                x, y, w, h = roi; cx, cy = x - 1 + w / 2, y - 1 + h / 2; ww, hh = 2.2 * w, 2.2 * h
            else:
                cx, cy, ww, hh = W / 2, H / 2, W, H
            x0 = int(np.clip(cx - ww / 2, 0, W - 1)); y0 = int(np.clip(cy - hh / 2, 0, H - 1))
            x1 = int(np.clip(cx + ww / 2, 1, W)); y1 = int(np.clip(cy + hh / 2, 1, H))
            S = F[y0:y1, x0:x1].astype(float); lo, hi = np.percentile(S, [1, 75])
            im = Image.fromarray((np.clip((S - lo) / max(1, hi - lo), 0, 1) * 255).astype(np.uint8)).convert("RGB")
            sc = min(TW / im.size[0], TH / im.size[1]); im = im.resize((int(im.size[0] * sc), int(im.size[1] * sc)), Image.BICUBIC)
            di = ImageDraw.Draw(im)
            if roi is not None:
                di.rectangle([(x - 1 - x0) * sc, (y - 1 - y0) * sc, (x - 1 - x0 + w) * sc, (y - 1 - y0 + h) * sc], outline=(255, 220, 0), width=2)
            di.text((4, 4), f"t = {frac * d:.0f} s", fill=(255, 0, 0))
            sheet.paste(im, (c * (TW + 10), y_row + LAB))
    out.parent.mkdir(parents=True, exist_ok=True); sheet.save(out)
    print("ROI check:", out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
