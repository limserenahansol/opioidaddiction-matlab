"""wholeeye_teach_export.py - frames for drawing the WHOLE EYE outline of one mouse (9670f black, eye damage since
day 1; user 2026-11-17: "we can just use whole eye as alternate"; at most 3 frames per session).
Per session: 72 frames spread over the session -> the 3 with the darkest / middle / brightest eye box (so the dark
slit, the in-between and the cloudy eye are all drawn), at least 3 min apart. The crop = the package's eye box
x1.5 around its centre. Runs at Idle priority.
Writes <out>/frames/*.png, <out>/frames.json and, with --page, <out>/index.html for tools/pupil_teach_server.py.
usage: python wholeeye_teach_export.py --project P --cage 9670f --mouse black --days 1 2 4 --out DIR --page PAGE"""
import argparse
import ctypes
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from PIL import Image

try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)   # Idle
except Exception:
    pass
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pupil_example_sheets as PES   # noqa: E402  grab(), stretch(), roi_crop()

FPS_REAL = 30.26


def find_video(project, d, cage, mouse):
    vd = project / "video" / f"day{d}" / cage / mouse
    c = sorted(vd.glob("*.mp4")) if vd.is_dir() else []
    if not c:   # not filed yet (today): IC Capture saves to the Documents folder
        c = sorted(Path.home().joinpath("Documents").glob(f"{cage}{mouse}_*day{d}*.mp4"))
    return c[0] if c else None


def find_roi(project, d, cage, mouse):
    """the package's eye box of that day, else of the nearest earlier day that has one"""
    for dd in range(d, 0, -1):
        for root in ("pupil_tst_teach", "pupil_tst_final"):
            q = project / root / f"day{dd}" / cage / mouse / "pupil_outputs_general" / "qc.json"
            if q.exists():
                try:
                    return json.loads(q.read_text(encoding="utf-8"))["roi"], dd
                except Exception:
                    continue
        r = project / "video" / f"day{dd}" / cage / mouse / "roi.mat"
        if dd == d and r.exists():   # a new day: its own eye box, set before the package has run
            from scipy.io import loadmat
            return [float(v) for v in loadmat(r)["roi"].ravel()], dd
    return None, None


def n_frames(video):
    r = subprocess.run([PES.ffmpeg(), "-hide_banner", "-i", str(video)], capture_output=True, text=True, errors="replace")
    import re
    m = re.search(r"Duration: (\d+):(\d+):([\d.]+)", r.stderr)
    return int(round((int(m.group(1)) * 3600 + int(m.group(2)) * 60 + float(m.group(3))) * 60)) if m else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True)
    ap.add_argument("--cage", required=True)
    ap.add_argument("--mouse", required=True)
    ap.add_argument("--days", type=int, nargs="+", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--page")
    ap.add_argument("--frames", type=int, nargs="+",
                    help="exact frame numbers (one day only) instead of the darkest / middle / brightest pick "
                         "(added 2026-11-18: late day-1 frames of 9670f black, where the eye turns white/cloudy)")
    ap.add_argument("--reason", default="", help="text shown with the --frames frames")
    a = ap.parse_args(argv)
    if a.frames and len(a.days) != 1:
        ap.error("--frames needs exactly one day in --days")
    project, out = Path(a.project), Path(a.out)
    (out / "frames").mkdir(parents=True, exist_ok=True)
    items, sessions = [], []
    for d in a.days:
        video = find_video(project, d, a.cage, a.mouse)
        roi, roi_day = find_roi(project, d, a.cage, a.mouse)
        if video is None or roi is None:
            print(f"  day{d}: no video or no eye box - skipped"); continue
        n = n_frames(video)
        x0, y0, cw, ch = PES.roi_crop(roi)
        W, H = int(round(cw * 1.5)), int(round(ch * 1.5))
        X0 = max(0, min(1920 - W, int(round(x0 + cw / 2 - W / 2))))
        Y0 = max(0, min(1080 - H, int(round(y0 + ch / 2 - H / 2))))
        crop = (X0, Y0, W, H)
        ks = [] if a.frames else np.linspace(0.03 * n, 0.97 * n, 72).astype(int)
        cand = []
        for k in ks:
            img = PES.grab(video, crop, int(k))
            if img is None:
                continue
            inner = img[(y0 - Y0):(y0 - Y0 + ch), (x0 - X0):(x0 - X0 + cw)]
            cand.append((float(inner.mean()), int(k), img))
        cand.sort(key=lambda c: c[0])
        picked = []
        for k in (a.frames or []):
            img = PES.grab(video, crop, int(k))
            if img is not None:
                picked.append((float("nan"), int(k), img, None))
        for q, why in (() if a.frames else ((0.10, "어두운 편"), (0.50, "중간"), (0.90, "밝은 편"))):
            order = sorted(range(len(cand)), key=lambda i: abs(i - q * (len(cand) - 1)))
            for i in order:
                if all(abs(cand[i][1] - p[1]) >= 3 * 60 * FPS_REAL for p in picked):
                    picked.append((cand[i][0], cand[i][1], cand[i][2], why)); break
        session = video.stem
        for mean, k, img, why in sorted(picked, key=lambda p: p[1]):
            fid = f"d{d}_{a.cage}_{a.mouse}_f{k}"
            Image.fromarray(PES.stretch(img)).save(out / "frames" / f"{fid}.png")
            items.append(dict(id=fid, session=session, day=d, cage=a.cage, mouse=a.mouse, frame=k, image_frame=k,
                              t_s=round((k - 1) / FPS_REAL, 2),
                              reason=(f"눈 전체 그리기 · 눈 상자가 {why} 프레임" if why else
                                      (a.reason or "눈 전체 그리기 · 지정한 프레임")),
                              run_s=0, w=W, h=H, img=f"frames/{fid}.png", table=None, table_src=None, be=None,
                              be_src=None, ask=True, crop_full=[X0, Y0, W, H], eye_box_full=[x0, y0, cw, ch],
                              box_day=roi_day, video=str(video)))
        sessions.append(dict(session=session, day=d, cage=a.cage, mouse=a.mouse, frames=n, picked=len(picked),
                             not_found_pct="-", wrong_pct="-"))
        print(f"  day{d} {session}: {n} frames, eye box from day {roi_day}, picked {[p[1] for p in picked]}")
    data = dict(made_by="tools/wholeeye_teach_export.py", fps_real=FPS_REAL,
                coords="crop pixels, 0-based: x right, y down; centre of the top-left pixel = (0, 0); "
                       "crop_full = [x0, y0, w, h] of the crop in the 1920x1080 frame (0-based)",
                sessions=sessions, frames=items)
    (out / "frames.json").write_text(json.dumps(data, indent=1, ensure_ascii=False), encoding="utf-8")
    if a.page:
        page = Path(a.page).read_text(encoding="utf-8").replace("/*FRAMES_JSON*/", json.dumps(data).replace("</", "<\\/"))
        (out / "index.html").write_text(
            '<!doctype html><html lang="ko"><head><meta charset="utf-8"><meta name="viewport" '
            'content="width=device-width,initial-scale=1,viewport-fit=cover"><style>:root{color-scheme:light;'
            'padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}'
            'body{margin:0;font:14px system-ui,sans-serif}img{max-width:100%}[hidden]{display:none!important}'
            '</style></head><body>' + page + '</body></html>', encoding="utf-8")
    print(f"  {len(items)} frames -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
