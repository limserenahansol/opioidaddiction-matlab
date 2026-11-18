"""wholeeye_marks_to_outlines.py - turn the outlines the user drew on the whole-eye drawing page (marks.json of a
folder made by wholeeye_teach_export.py, served by pupil_teach_server.py) into a user_outlines_<name>.json that
wholeeye_build_set.py reads (2026-11-18: 9670f black, late day-1 frames, "the white part is eye").
Page polygons are in crop pixels (0-based); outline_full = polygon + crop_full[x0, y0] (0-based full-frame pixels),
the same convention as wholeeye_user_outlines.py / wholeeye_strokes.py. Answers other than a drawing (e.g. 'eye
fully closed') are listed and skipped. Read-only for the page folder.
usage: python wholeeye_marks_to_outlines.py --folder analysis/wholeeye/teach_d1_white --name 9670f_d1white"""
import argparse
import json
from pathlib import Path

T = Path(__file__).resolve().parent
P = T.parent


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--name", required=True)
    a = ap.parse_args(argv)
    F = Path(a.folder) if Path(a.folder).is_absolute() else P / a.folder
    marks = json.loads((F / "marks.json").read_text(encoding="utf-8"))
    frames = {f["id"]: f for f in json.loads((F / "frames.json").read_text(encoding="utf-8"))["frames"]}
    out = []
    for fid, m in sorted(marks.items()):
        fr = frames.get(fid)
        poly = (m.get("ellipse") or {}).get("poly") if isinstance(m.get("ellipse"), dict) else None
        if fr is None or m.get("choice") != "draw" or not poly or len(poly) < 3:
            print(f"  {fid}: answer '{m.get('choice')}' - no outline, skipped")
            continue
        x0, y0 = fr["crop_full"][0], fr["crop_full"][1]
        out.append(dict(day=int(fr["day"]), cage=fr["cage"], mouse=fr["mouse"], frame=int(fr["frame"]),
                        video=fr["video"], outline_full=[[round(x + x0, 2), round(y + y0, 2)] for x, y in poly],
                        tile=None, closed_loop=True, drawn_at=m.get("at"),
                        read_from=f"drawing page {F.relative_to(P) if P in F.parents else F}\\marks.json "
                                  "(the user's free-hand outline)"))
        print(f"  {fid}: {len(poly)} points, area {int(m['ellipse'].get('area', 0))} px")
    dst = P / "analysis" / "wholeeye" / f"user_outlines_{a.name}.json"
    dst.write_text(json.dumps(out, indent=1), encoding="utf-8")
    print(f"{len(out)} outline(s) -> {dst}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
