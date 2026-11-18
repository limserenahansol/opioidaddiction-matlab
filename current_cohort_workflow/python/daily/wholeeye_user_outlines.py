"""wholeeye_user_outlines.py - the user's whole-eye outlines of 9670f black (2026-11-17), drawn on the frame sheet
analysis/wholeeye/9670f_black_frames.png shown in the browser pane (orange strokes on a 1051x965 screenshot).
Claude read the outline points from that screenshot by eye (about +-3 screenshot px = about +-4 video px) and maps
them to video-frame pixels here. The sheet itself is exactly known (scratchpad black_eye_frames.py):
  rows day1 / day2 / day4; tiles at 10/35/60/85 % of the session; tile = the eye box x2 around its centre, shown x2;
  6-px gaps; the sheet (2898x1644) was shown 1051 px wide, centred vertically in the 965-px screenshot.
Output: analysis/wholeeye/user_outlines_20261117.json (full-frame pixel coordinates, 0-based pixel centres) and
analysis/wholeeye/user_outlines_check.png (the outlines drawn back on the real frames, for the user to confirm)."""
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import pupil_example_sheets as PES   # noqa: E402

P = T.parent
SCALE = 2898 / 1051                    # sheet px per screenshot px
TOP = (965 - 1644 / SCALE) / 2         # screenshot y of the sheet's top edge
SESS = {  # day: (video, eye box x, y, w, h (package, 1-based), frames in the session, sheet row top)
    1: (P / "video/day1/9670f/black/9670fblack__fr3day1.mp4", (1186, 236, 150, 120), 27537, 0),
    2: (P / "video/day2/9670f/black/9670fblack_fr5day2_2.mp4", (1100, 215, 180, 144), 27536, 486),
    4: (Path.home() / "Documents/9670fblack_prday4.mp4", (1100, 215, 180, 144), 27537, 1068),
}
if not SESS[4][0].exists():
    SESS[4] = (P / "video/day4/9670f/black/9670fblack_prday4.mp4",) + SESS[4][1:]
QS = (0.10, 0.35, 0.60, 0.85)
# (day, tile index 0-3): outline points read from the screenshot (x, y in screenshot px)
READ = {
    (1, 1): [(312, 277), (318, 266), (335, 256), (355, 250), (375, 248), (382, 252), (390, 257), (391, 267), (382, 268),
             (370, 274), (350, 280), (330, 284), (315, 283)],
    (1, 2): [(494, 300), (498, 288), (512, 272), (528, 256), (540, 247), (556, 240), (566, 242), (569, 256), (565, 270),
             (552, 278), (530, 286), (510, 293)],
    (2, 0): [(99, 478), (102, 465), (104, 451), (120, 449), (140, 449), (152, 450), (155, 455), (148, 461), (142, 466),
             (125, 470), (110, 474)],
    (2, 1): [(349, 474), (358, 464), (372, 454), (390, 449), (410, 448), (425, 452), (428, 458), (421, 464), (400, 469),
             (375, 472)],
    (2, 3): [(871, 479), (876, 463), (888, 450), (910, 443), (935, 442), (960, 445), (972, 452), (974, 458), (962, 466),
             (940, 476), (915, 485), (890, 487), (876, 484)],
    (4, 0): [(124, 672), (130, 656), (140, 648), (155, 645), (172, 648), (186, 654), (184, 660), (170, 665), (150, 670),
             (132, 673)],
    (4, 1): [(377, 675), (383, 661), (393, 648), (408, 640), (420, 638), (430, 644), (442, 652), (455, 657), (452, 664),
             (432, 669), (405, 673), (385, 676)],
    (4, 2): [(638, 660), (641, 645), (650, 636), (665, 633), (680, 636), (697, 641), (713, 644), (718, 650), (714, 656),
             (695, 659), (668, 660), (645, 661)],
    (4, 3): [(896, 676), (900, 665), (912, 653), (928, 647), (945, 643), (958, 645), (963, 652), (958, 659), (940, 665),
             (918, 671), (900, 676)],
}


def crop_origin(box):
    x, y, w, h = box
    cx, cy = x - 1 + w / 2, y - 1 + h / 2
    return int(round(cx - w)), int(round(cy - h)), 2 * w, 2 * h      # the sheet's tile = the eye box x2


out, tiles = [], []
for (d, i), pts in READ.items():
    video, box, n, row_top = SESS[d]
    x0, y0, cw, ch = crop_origin(box)
    tile_w = 2 * cw                                                   # tiles are shown x2
    tx0 = i * (tile_w + 6)
    k = int(QS[i] * n)
    full = []
    for sx, sy in pts:
        gx, gy = sx * SCALE, (sy - TOP) * SCALE                       # sheet px
        cx_, cy_ = (gx - tx0) / 2, (gy - row_top) / 2                 # crop px (the tile is the crop x2)
        full.append([round(x0 + cx_, 1), round(y0 + cy_, 1)])        # full-frame px, 0-based
    out.append(dict(day=d, cage="9670f", mouse="black", frame=k, video=str(video), outline_full=full,
                    tile=i, read_from="screenshot of the frame sheet (Claude read the user's orange strokes)"))
    # check image: the eye box x1.5 around the outline, x3, the outline in magenta
    xs, ys = [p[0] for p in full], [p[1] for p in full]
    mx, my = (min(xs) + max(xs)) / 2, (min(ys) + max(ys)) / 2
    W, H = 240, 150
    X0, Y0 = int(mx - W / 2), int(my - H / 2)
    img = PES.grab(video, (X0, Y0, W, H), k)
    im = Image.fromarray(PES.stretch(img)).resize((W * 3, H * 3), Image.BILINEAR).convert("RGB")
    dr = ImageDraw.Draw(im)
    dr.line([((p[0] - X0 + 0.5) * 3, (p[1] - Y0 + 0.5) * 3) for p in full + [full[0]]], fill=(255, 60, 200), width=3)
    dr.text((8, 6), f"day{d} frame {k}", fill=(255, 255, 0), font=PES.font(22))
    tiles.append(im)
(P / "analysis/wholeeye").mkdir(parents=True, exist_ok=True)
(P / "analysis/wholeeye/user_outlines_20261117.json").write_text(json.dumps(out, indent=1), encoding="utf-8")
cols = 3
Wt, Ht = tiles[0].size
sheet = Image.new("RGB", (cols * Wt + (cols - 1) * 8, 3 * Ht + 16), (25, 25, 25))
for j, t in enumerate(tiles):
    sheet.paste(t, ((j % cols) * (Wt + 8), (j // cols) * (Ht + 8)))
sheet.save(P / "analysis/wholeeye/user_outlines_check.png")
print("outlines:", len(out), "->", P / "analysis/wholeeye/user_outlines_20261117.json")
for o in out:
    xs, ys = [p[0] for p in o["outline_full"]], [p[1] for p in o["outline_full"]]
    print(f"  day{o['day']} f{o['frame']}: x {min(xs):.0f}-{max(xs):.0f}, y {min(ys):.0f}-{max(ys):.0f}")
