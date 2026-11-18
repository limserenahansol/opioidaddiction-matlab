"""wholeeye_frame_sheet.py - a frame sheet for drawing the whole-eye outline of one mouse by hand (user 2026-11-17:
9670f black and 9669m black have eye damage; the user draws on this picture in the browser pane).
One row per day, 4 frames (10/35/60/85 % of the session); a tile = the package's eye box x2 around its centre,
shown x2; 6-px gaps. The exact layout goes to <name>_layout.json so the drawn strokes can be mapped back to frame
pixels. Idle priority.
usage: python wholeeye_frame_sheet.py --cage 9669m --mouse black --days 1 2 3"""
import argparse
import ctypes
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

try:
    ctypes.windll.kernel32.SetPriorityClass(ctypes.windll.kernel32.GetCurrentProcess(), 0x40)
except Exception:
    pass
T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import pupil_example_sheets as PES          # noqa: E402
import wholeeye_teach_export as WE          # noqa: E402

P = T.parent
QS = (0.10, 0.35, 0.60, 0.85)
S = 2


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cage", required=True)
    ap.add_argument("--mouse", required=True)
    ap.add_argument("--days", type=int, nargs="+", required=True)
    a = ap.parse_args(argv)
    rows, layout = [], []
    for d in a.days:
        video = WE.find_video(P, d, a.cage, a.mouse)
        roi, roi_day = WE.find_roi(P, d, a.cage, a.mouse)
        if video is None or roi is None:
            print(f"  day{d}: no video or eye box - skipped"); continue
        n = WE.n_frames(video)
        x, y, w, h = [int(round(v)) for v in roi]
        cx, cy = x - 1 + w / 2, y - 1 + h / 2
        cw, ch = 2 * w, 2 * h
        x0, y0 = int(round(cx - cw / 2)), int(round(cy - ch / 2))
        tiles, frames = [], []
        for q in QS:
            k = int(q * n)
            a_ = PES.grab(video, (x0, y0, cw, ch), k).astype(float)
            lo, hi = np.percentile(a_, 1), np.percentile(a_, 99.5)
            a_ = np.clip((a_ - lo) / max(1.0, hi - lo), 0, 1) ** 0.8
            im = Image.fromarray((a_ * 255).astype(np.uint8)).resize((cw * S, ch * S), Image.BILINEAR).convert("RGB")
            dr = ImageDraw.Draw(im)
            bx, by = (x - 1 - x0) * S, (y - 1 - y0) * S
            dr.rectangle([bx, by, bx + w * S, by + h * S], outline=(255, 200, 0), width=2)
            dr.text((6, 4), f"day{d} frame {k} ({int(q * 100)}%)", fill=(255, 255, 0), font=PES.font(18))
            tiles.append(im); frames.append(k)
        rows.append(tiles)
        layout.append(dict(day=d, video=str(video), eye_box=[x, y, w, h], box_day=roi_day, crop0=[x0, y0, cw, ch],
                           tile_w=cw * S, tile_h=ch * S, frames=frames))
    W = max(sum(t.width for t in r) + 6 * (len(r) - 1) for r in rows)
    H = sum(r[0].height for r in rows) + 6 * (len(rows) - 1)
    sheet = Image.new("RGB", (W, H), (30, 30, 30))
    yy = 0
    for r, L in zip(rows, layout):
        xx = 0
        L["row_top"] = yy
        L["tile_x"] = []
        for t in r:
            sheet.paste(t, (xx, yy)); L["tile_x"].append(xx); xx += t.width + 6
        yy += r[0].height + 6
    name = f"{a.cage}_{a.mouse}_frames"
    out = P / "analysis" / "wholeeye" / f"{name}.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(out)
    (P / "analysis" / "wholeeye" / f"{name}_layout.json").write_text(
        json.dumps(dict(sheet=str(out), size=[W, H], scale_shown=S, rows=layout), indent=1), encoding="utf-8")
    print(out, sheet.size)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
