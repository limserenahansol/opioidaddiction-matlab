"""Put tail-suspension (TST) and hot-plate (HP) videos of one experiment day into the cohort folders.

    python arrange_assay_videos.py --project <cohort folder> --day 5 --src <folder with the day's assay videos> [--apply]

  <project>\\tst_video\\day<N>\\<cage>\\<mouse>\\<original name>   (TST: the pupil/TST package scores every video with
                                                                 'tst' in its path under tst_video_root)
  <project>\\hp_video\\day<N>\\<cage>\\<mouse>\\<original name>    (HP: scored by hand; the ManualScoringResults*hot*.mat
                                                                 goes into the same folder)
The assay is read from the file name ('tst' or 'hot'/'hp'), the mouse from '<4 digits><sex letter><colour>' as in the
task-log names (9670fgreen1_tst.mp4, 9669m_black_hp_day5.mov). A name that gives no mouse, or a colour that is not a
known mouse folder of that cage, is listed and NOT copied. Files are copied (originals stay). Dry run unless --apply.
Use the experiment day the assay belongs to (the Step 02 day), e.g. TST/HP after the day-5 session -> --day 5.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
from pathlib import Path

VIDEO_EXT = {".mp4", ".avi", ".mov", ".mkv", ".wmv", ".m4v"}
MOUSE = re.compile(r"(?P<cage>\d{4}[mfMF]?)[\s_-]*(?P<color>black|orange|red|white|green\d*|blue|yellow|pink|purple)", re.I)


def known_mice(project: Path):
    """cage -> set of mouse folder names, from the arranged task logs (python\\dayN\\<cage>\\<mouse>)."""
    out = {}
    for p in (project / "python").glob("day*/*/*"):
        if p.is_dir():
            out.setdefault(p.parent.name.lower(), set()).add(p.name.lower())
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project", required=True, type=Path)
    ap.add_argument("--day", required=True, type=int)
    ap.add_argument("--src", required=True, type=Path)
    ap.add_argument("--apply", action="store_true")
    a = ap.parse_args(argv)
    mice = known_mice(a.project)
    n = 0
    for v in sorted(p for p in a.src.rglob("*") if p.is_file() and p.suffix.lower() in VIDEO_EXT):
        low = v.name.lower()
        assay = "tst" if "tst" in low else ("hp" if re.search(r"(hot|(^|[^a-z])hp([^a-z]|$))", low) else None)
        m = MOUSE.search(v.stem)
        if assay is None or m is None:
            print(f"SKIP {v.name}: {'no tst/hot/hp in the name' if assay is None else 'no <cage><sex><colour> in the name'}")
            continue
        cage, color = m.group("cage").lower(), m.group("color").lower()
        cage_dirs = {c for c in mice if c[:4] == cage[:4]}
        cage_dir = cage if cage in mice else (next(iter(cage_dirs)) if len(cage_dirs) == 1 else None)
        if cage_dir is None or color not in mice.get(cage_dir, set()):
            cands = sorted(x for x in mice.get(cage_dir or "", set()) if x.startswith(color))
            if cage_dir and len(cands) == 1:
                print(f"NOTE {v.name}: colour '{color}' read as '{cands[0]}' (the only {cage_dir} mouse starting so)")
                color = cands[0]
            else:
                print(f"SKIP {v.name}: {cage}/{color} is not a known mouse folder ({sorted(mice.get(cage_dir or '', []))})")
                continue
        dst = a.project / f"{assay}_video" / f"day{a.day}" / cage_dir / color / v.name
        if dst.exists():
            print(f"have {dst.relative_to(a.project)}")
            continue
        print(f"{'copy' if a.apply else 'would copy'} {v.name} -> {dst.relative_to(a.project)}")
        if a.apply:
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(v, dst)
        n += 1
    print(f"{n} video(s) {'copied' if a.apply else 'to copy (dry run)'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
