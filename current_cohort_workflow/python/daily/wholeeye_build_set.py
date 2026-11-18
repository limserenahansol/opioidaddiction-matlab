"""wholeeye_build_set.py - ONE training set (and the tracking crops) for the whole-eye model of the mice whose pupil
cannot be seen (user 2026-11-17: 9670f black and 9669m black have eye damage; one model for both).
The crop of a session = the package's eye box of that day (or the nearest earlier day with one) x1.4 around its
centre, so the whole eye fits even where it reaches past the box. The same crop is used for training and tracking.
Inputs: analysis/wholeeye/user_outlines_*.json (tools/wholeeye_user_outlines.py, tools/wholeeye_strokes.py).
Outputs: analysis/wholeeye/train_set/ (sample_NN.png = raw grey crop, mask_NN.png, meta.json) and
         analysis/wholeeye/crops.json ({"<cage>_<mouse>_day<N>": {cage, mouse, day, video, crop [x0 y0 w h] 0-based,
         box_day}}) for every day 1..--last-day that has a finished video.
usage: python wholeeye_build_set.py [--mice 9670f/black 9669m/black] [--last-day 4]"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw

T = Path(__file__).resolve().parent
sys.path.insert(0, str(T))
import pupil_example_sheets as PES          # noqa: E402
import wholeeye_teach_export as WE          # noqa: E402  find_video(), find_roi()

P = T.parent
GROW = 1.4


def session_crop(d, cage, mouse):
    video = WE.find_video(P, d, cage, mouse)
    roi, box_day = WE.find_roi(P, d, cage, mouse)
    if video is None or roi is None:
        return None
    x0, y0, cw, ch = PES.roi_crop(roi)
    W, H = int(round(cw * GROW)), int(round(ch * GROW))
    X0 = max(0, min(1920 - W, int(round(x0 + cw / 2 - W / 2))))
    Y0 = max(0, min(1080 - H, int(round(y0 + ch / 2 - H / 2))))
    return dict(cage=cage, mouse=mouse, day=d, video=str(video), crop=[X0, Y0, W, H], box_day=box_day,
                eye_box=[int(v) for v in roi])


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--mice", nargs="+", help="cage/mouse ... (default: the mice approved in <pupil tree>/"
                                              "wholeeye_mice.csv, else 9670f/black 9669m/black)")
    ap.add_argument("--last-day", type=int, help="default: every day folder in video/")
    ap.add_argument("--crops-only", action="store_true", help="only crops.json (the daily routine; training set kept)")
    a = ap.parse_args(argv)
    if not a.mice:
        sys.path.insert(0, str(T))
        import concat_choice as CC
        rows = CC.read_rows(CC.active_pupil_root(P) / CC.WHOLE_MICE, warn=None, need=("cage", "mouse"))
        a.mice = [f"{r['cage'].strip()}/{r['mouse'].strip()}" for r in rows if CC.is_approved(r.get("approved"))] \
            or ["9670f/black", "9669m/black"]
    last = a.last_day or max([int(p.name[3:]) for p in (P / "video").glob("day*") if p.name[3:].isdigit()] or [4])
    crops = {}
    for cm in a.mice:
        cage, mouse = cm.split("/")
        for d in range(1, last + 1):
            c = session_crop(d, cage, mouse)
            if c:
                crops[f"{cage}_{mouse}_day{d}"] = c
    (P / "analysis/wholeeye/crops.json").write_text(json.dumps(crops, indent=1), encoding="utf-8")
    if a.crops_only:
        print("tracking crops:", {k: v["crop"] for k, v in crops.items()})
        return 0
    out = P / "analysis/wholeeye/train_set"
    out.mkdir(parents=True, exist_ok=True)
    for f in out.glob("*.png"):
        f.unlink()
    outlines = []
    for f in sorted((P / "analysis/wholeeye").glob("user_outlines_*.json")):
        for o in json.loads(f.read_text(encoding="utf-8")):
            o.setdefault("cage", "9670f"); o.setdefault("mouse", "black")
            o["source_file"] = f.name
            outlines.append(o)
    meta = []
    for j, o in enumerate(outlines, 1):
        c = crops.get(f"{o['cage']}_{o['mouse']}_day{o['day']}")
        if c is None:
            print(f"  no crop for {o['cage']} {o['mouse']} day{o['day']} - skipped"); continue
        X0, Y0, W, H = c["crop"]
        img = PES.grab(Path(c["video"]), (X0, Y0, W, H), int(o["frame"]))
        m = Image.new("L", (W, H), 0)
        ImageDraw.Draw(m).polygon([(x - X0, y - Y0) for x, y in o["outline_full"]], fill=255, outline=255)
        Image.fromarray(img).save(out / f"sample_{j:02d}.png")
        m.save(out / f"mask_{j:02d}.png")
        meta.append(dict(img=f"sample_{j:02d}.png", mask=f"mask_{j:02d}.png", cage=o["cage"], mouse=o["mouse"],
                         day=o["day"], frame=o["frame"], crop=c["crop"], eye_px=int((np.array(m) > 0).sum()),
                         source=o["source_file"]))
        print(f"  sample {j}: {o['cage']} {o['mouse']} day{o['day']} frame {o['frame']}, eye {meta[-1]['eye_px']} px")
    (out / "meta.json").write_text(json.dumps(meta, indent=1), encoding="utf-8")
    print("tracking crops:", {k: v["crop"] for k, v in crops.items()})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
