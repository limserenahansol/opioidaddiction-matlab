#!/usr/bin/env python3
r"""
render_deck_lo.py - slide previews with LibreOffice headless and its OWN profile (made 2026-11-19).

The user's rule (mrpe package docs/07, rule 23): "Slides: build with python-pptx only. Never drive PowerPoint through
Office COM: it closes the user's open documents. Render previews with LibreOffice headless, using its own profile."
LibreOffice was installed by the user on 2026-11-19 ("I installed").

  python render_deck_lo.py DECK.pptx [--out FOLDER] [--png]
    -> FOLDER/DECK.pdf (every slide; open it to check the layout)
    --png also writes one PNG per slide: each slide is saved as its own 1-slide copy (python-pptx) and converted, so
          a fix can be sent as one slide + its PNG (rule 24). The deck itself is never changed.
The profile lives in %LOCALAPPDATA%\claude_lo_profile (separate from the user's LibreOffice), and the call never
touches PowerPoint. Exit code 0 when the PDF (and PNGs) were made.
"""
from __future__ import annotations

import argparse
import copy
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

SOFFICE = Path(r"C:\Program Files\LibreOffice\program\soffice.exe")
PROFILE = Path(os.environ.get("LOCALAPPDATA", tempfile.gettempdir())) / "claude_lo_profile"


def convert(src: Path, fmt: str, out: Path) -> Path | None:
    out.mkdir(parents=True, exist_ok=True)
    cmd = [str(SOFFICE), f"-env:UserInstallation={PROFILE.as_uri()}", "--headless", "--norestore", "--nolockcheck",
           "--convert-to", fmt, "--outdir", str(out), str(src)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    dst = out / (src.stem + "." + fmt.split(":")[0])
    if r.returncode != 0 or not dst.exists():
        print(f"  !! LibreOffice could not convert {src.name} to {fmt}: {r.stdout.strip()} {r.stderr.strip()}")
        return None
    return dst


def one_slide_copies(deck: Path, tmp: Path) -> list:
    """one pptx per slide (the other slides removed from a copy); the original deck is not touched"""
    from pptx import Presentation
    n = len(Presentation(str(deck)).slides)
    paths = []
    for k in range(n):
        prs = Presentation(str(deck))
        ids = prs.slides._sldIdLst
        for i, sid in reversed(list(enumerate(list(ids)))):
            if i != k:
                prs.part.drop_rel(sid.rId)
                ids.remove(sid)
        p = tmp / f"slide{k + 1}.pptx"
        prs.save(str(p))
        paths.append(p)
    return paths


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("deck")
    ap.add_argument("--out", default=None)
    ap.add_argument("--png", action="store_true")
    a = ap.parse_args(argv)
    deck = Path(a.deck).resolve()
    if not SOFFICE.exists():
        print(f"  !! LibreOffice not found at {SOFFICE}")
        return 2
    out = Path(a.out).resolve() if a.out else deck.parent / f"{deck.stem}_preview"
    pdf = convert(deck, "pdf", out)
    if pdf is None:
        return 1
    print(pdf)
    if a.png:
        with tempfile.TemporaryDirectory() as td:
            hi = ('png:impress_png_Export:{"PixelWidth":{"type":"long","value":"1920"},'
                  '"PixelHeight":{"type":"long","value":"1080"}}')
            for p in one_slide_copies(deck, Path(td)):
                png = convert(p, hi, out) or convert(p, "png", out)     # 1920 x 1080, else LibreOffice's default
                if png is None:
                    return 1
                print(png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
