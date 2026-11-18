"""wholeeye_strokes.py - the user's whole-eye outlines, taken from a browser-pane screenshot of a frame sheet
(tools/wholeeye_frame_sheet.py) on which the user drew with the pane's annotation pen (orange strokes).
The sheet's place in the screenshot is found from the screenshot itself (the non-black block); every tile's orange
pixels -> small gaps closed -> the loop filled -> eroded by about half a stroke (= the stroke's centre line) -> the
outline polygon in video-frame pixels (0-based pixel centres). A tile without orange is not used.
Writes analysis/wholeeye/user_outlines_<cage>_<mouse>.json and ..._check.png (outlines back on the real frames).
usage: python wholeeye_strokes.py --shot SCREENSHOT.png --layout analysis/wholeeye/9669m_black_frames_layout.json
       --cage 9669m --mouse black"""
import argparse
import json
import sys
from pathlib import Path

import cv2
import numpy as np
from PIL import Image, ImageDraw

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import pupil_example_sheets as PES   # noqa: E402

P = T.parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--shot", required=True)
    ap.add_argument("--layout", required=True)
    ap.add_argument("--cage", required=True)
    ap.add_argument("--mouse", required=True)
    a = ap.parse_args(argv)
    L = json.loads(Path(a.layout).read_text(encoding="utf-8"))
    W, H = L["size"]
    rgb = np.array(Image.open(a.shot).convert("RGB")).astype(int)
    hsv = cv2.cvtColor(rgb.astype(np.uint8), cv2.COLOR_RGB2HSV)            # H 0-179
    orange = (hsv[..., 0] <= 15) & (hsv[..., 1] >= 100) & (hsv[..., 2] >= 110)
    # the sheet = the block of rows/columns that are not the near-black page background
    content = (rgb.max(axis=2) > 28) & ~orange
    rows = np.where(content.mean(axis=1) > 0.25)[0]
    y0s, y1s = rows.min(), rows.max() + 1
    cols = np.where(content[y0s:y1s].mean(axis=0) > 0.5)[0]
    x0s, x1s = cols.min(), cols.max() + 1
    sx, sy = W / (x1s - x0s), H / (y1s - y0s)
    print(f"sheet in the screenshot: x {x0s}-{x1s}, y {y0s}-{y1s}; scale {sx:.4f} / {sy:.4f} sheet px per screenshot px")
    scale = (sx + sy) / 2
    out, tiles = [], []
    for R in L["rows"]:
        X0c, Y0c, cw, ch = R["crop0"]
        S = L["scale_shown"]
        for i, k in enumerate(R["frames"]):
            # this tile in screenshot pixels
            tx0 = x0s + R["tile_x"][i] / scale
            ty0 = y0s + R["row_top"] / scale
            tx1, ty1 = tx0 + R["tile_w"] / scale, ty0 + R["tile_h"] / scale
            m = np.zeros(orange.shape, np.uint8)
            ys, xs = np.nonzero(orange)
            sel = (xs >= tx0) & (xs < tx1) & (ys >= ty0) & (ys < ty1)
            if sel.sum() < 15:
                continue
            m[ys[sel], xs[sel]] = 1
            stroke = cv2.morphologyEx(m, cv2.MORPH_CLOSE, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5)))
            cs, _ = cv2.findContours(stroke, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            filled = np.zeros_like(stroke)
            cv2.drawContours(filled, [max(cs, key=cv2.contourArea)], -1, 1, -1)
            closed_loop = filled.sum() > 1.6 * stroke.sum()
            if not closed_loop:   # an open stroke: its convex hull
                hull = cv2.convexHull(np.column_stack([xs[sel], ys[sel]]).astype(np.int32))
                filled = np.zeros_like(stroke); cv2.fillPoly(filled, [hull], 1)
            half = cv2.erode(filled, cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (3, 3)))   # ~ the stroke centre
            cs, _ = cv2.findContours(half, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
            c = max(cs, key=cv2.contourArea)[:, 0, :].astype(float)
            # screenshot pixel centre -> sheet -> tile -> crop -> full frame (0-based pixel centres)
            gx = (c[:, 0] + 0.5 - x0s) * scale - R["tile_x"][i]
            gy = (c[:, 1] + 0.5 - y0s) * scale - R["row_top"]
            fx, fy = X0c + gx / S - 0.5, Y0c + gy / S - 0.5
            poly = [[round(float(x), 1), round(float(y), 1)] for x, y in zip(fx, fy)]
            out.append(dict(day=R["day"], cage=a.cage, mouse=a.mouse, frame=int(k), video=R["video"], outline_full=poly,
                            tile=i, closed_loop=bool(closed_loop),
                            read_from=f"orange pixels of {Path(a.shot).name} (the user's pen strokes)"))
            # check tile
            xs_, ys_ = fx, fy
            mx, my = (xs_.min() + xs_.max()) / 2, (ys_.min() + ys_.max()) / 2
            Wc, Hc = 240, 150
            CX0, CY0 = int(mx - Wc / 2), int(my - Hc / 2)
            img = PES.grab(Path(R["video"]), (CX0, CY0, Wc, Hc), int(k))
            im = Image.fromarray(PES.stretch(img)).resize((Wc * 3, Hc * 3), Image.BILINEAR).convert("RGB")
            dr = ImageDraw.Draw(im)
            dr.line([((x - CX0 + 0.5) * 3, (y - CY0 + 0.5) * 3) for x, y in poly + [poly[0]]], fill=(255, 60, 200), width=3)
            dr.text((8, 6), f"{a.cage} {a.mouse} day{R['day']} frame {k}" + ("" if closed_loop else " (open stroke: hull)"),
                    fill=(255, 255, 0), font=PES.font(22))
            tiles.append(im)
            print(f"  day{R['day']} f{k}: {len(poly)} points, x {fx.min():.0f}-{fx.max():.0f}, y {fy.min():.0f}-{fy.max():.0f}"
                  f"{'' if closed_loop else ' (open stroke)'}")
    base = P / "analysis" / "wholeeye" / f"user_outlines_{a.cage}_{a.mouse}"
    base.with_suffix(".json").write_text(json.dumps(out, indent=1), encoding="utf-8")
    if tiles:
        cols_ = 3
        tw, th = tiles[0].size
        nr = (len(tiles) + cols_ - 1) // cols_
        sheet = Image.new("RGB", (cols_ * tw + (cols_ - 1) * 8, nr * th + (nr - 1) * 8), (25, 25, 25))
        for j, t in enumerate(tiles):
            sheet.paste(t, ((j % cols_) * (tw + 8), (j // cols_) * (th + 8)))
        sheet.save(str(base) + "_check.png")
    print(f"{len(out)} outlines -> {base.with_suffix('.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
