#!/usr/bin/env python3
"""
pupil_teach_server.py - serves the pupil drawing page on this computer and saves every answer to a file (2026-11-16).

Why: the drawing page (tools/pupil_teach_page.html, filled by pupil_teach_export.py) could not be published as a
claude.ai artifact from this session (the organization asks a person to approve each new artifact). Run locally, the
page saves each answer straight into the project, where Claude reads it.

  GET  /                 the page folder (index.html, frames/*.png)
  GET  /api/marks        all answers: {frame_id: answer}
  PUT  /api/marks/<id>   save one answer (JSON object; <id> must be a frame of frames.json)
  DELETE /api/marks/<id> remove one answer
Answers: <folder>/marks.json. Every change is first appended to <folder>/marks_history.jsonl, then marks.json is
rewritten atomically (retried if Windows / OneDrive holds the file). marks.json is read STRICTLY: if it cannot be
read or is not a JSON object, nothing is written and the page gets an error (a broken file is never replaced by a
file with one answer). Listens on 127.0.0.1 only, one request at a time; a second copy cannot take the same port.

usage: python pupil_teach_server.py --folder <project>/analysis/pupil_teach [--port 8765]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import threading
import time
from datetime import datetime
from functools import partial
from http.server import HTTPServer, SimpleHTTPRequestHandler
from pathlib import Path

ID_RE = re.compile(r"^[A-Za-z0-9_.:@+~-]{1,120}$")
LOCK = threading.Lock()


class StoreError(Exception):
    pass


class Store:
    def __init__(self, folder: Path):
        self.folder = folder
        self.path = folder / "marks.json"
        self.hist = folder / "marks_history.jsonl"
        self.frames = folder / "frames.json"
        self._ids, self._mt = set(), None

    @property
    def ids(self) -> set:
        mt = self.frames.stat().st_mtime
        if mt != self._mt:       # the export was run again: read the frame ids again
            self._ids = {f["id"] for f in json.loads(self.frames.read_text(encoding="utf-8")).get("frames", [])}
            self._mt = mt
        return self._ids

    def load(self) -> dict:
        """strict: a missing file is empty; an unreadable or non-object file raises StoreError"""
        if not self.path.exists():
            return {}
        try:
            d = json.loads(self.path.read_text(encoding="utf-8-sig"))
        except (OSError, ValueError, UnicodeDecodeError) as e:
            raise StoreError(f"marks.json could not be read ({e}); nothing was changed")
        if not isinstance(d, dict):
            raise StoreError("marks.json is not a JSON object; nothing was changed")
        return d

    def change(self, fid: str, body):
        with LOCK:
            m = self.load()
            with open(self.hist, "a", encoding="utf-8") as fh:     # the history first: nothing is ever lost
                fh.write(json.dumps({"when": datetime.now().isoformat(timespec="seconds"), "id": fid,
                                     "body": body}, ensure_ascii=False) + "\n")
            if body is None:
                m.pop(fid, None)
            else:
                m[fid] = body
            tmp = self.path.with_suffix(".json.tmp")
            tmp.write_text(json.dumps(m, indent=1, ensure_ascii=False), encoding="utf-8")
            for i in range(40):
                try:
                    os.replace(tmp, self.path)
                    return
                except PermissionError:
                    time.sleep(0.05)            # OneDrive or a reader holds the file for a moment
            raise StoreError("marks.json is locked by another program; the answer is in marks_history.jsonl")


class Handler(SimpleHTTPRequestHandler):
    store: Store = None

    def log_message(self, fmt, *args):          # quiet: only changes are printed
        pass

    def _json(self, code: int, obj):
        b = json.dumps(obj, ensure_ascii=False).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Content-Length", str(len(b)))
        self.end_headers()
        self.wfile.write(b)

    def end_headers(self):                       # never cache the page or the list while drawing
        if not self.path.startswith("/frames/"):
            self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def _id(self):
        p = self.path.split("?")[0]
        if not p.startswith("/api/marks/"):
            return None
        fid = p[len("/api/marks/"):]
        return fid if ID_RE.match(fid) and fid in self.store.ids else ""

    def do_GET(self):
        if self.path.split("?")[0] == "/api/marks":
            try:
                return self._json(200, self.store.load())
            except StoreError as e:
                return self._json(500, {"error": str(e)})
        return super().do_GET()

    def _change(self, fid, body):
        try:
            self.store.change(fid, body)
        except (StoreError, OSError) as e:
            print(f"  !! {fid}: {e}", flush=True)
            return self._json(500, {"error": str(e)})
        print(f"  {'saved' if body is not None else 'removed'} {fid}" + (f": {body.get('choice')}" if body else ""),
              flush=True)
        return self._json(200, {"ok": True})

    def do_PUT(self):
        fid = self._id()
        if not fid:
            return self._json(404 if fid is None else 400, {"error": "unknown frame id"})
        try:
            n = int(self.headers.get("Content-Length", "0"))
            body = json.loads(self.rfile.read(n).decode("utf-8")) if 0 < n <= 65536 else None
        except (ValueError, UnicodeDecodeError):
            body = None
        if not isinstance(body, dict):
            return self._json(400, {"error": "the answer must be a JSON object (at most 64 KB)"})
        return self._change(fid, body)

    def do_DELETE(self):
        fid = self._id()
        if not fid:
            return self._json(404 if fid is None else 400, {"error": "unknown frame id"})
        return self._change(fid, None)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--folder", required=True)
    ap.add_argument("--port", type=int, default=8765)
    a = ap.parse_args(argv)
    folder = Path(a.folder).resolve()
    Handler.store = Store(folder)
    Handler.store.load()                         # fail now, not on the first answer, if marks.json is broken
    HTTPServer.allow_reuse_address = False       # Windows: SO_REUSEADDR would let a second copy share the port
    srv = HTTPServer(("127.0.0.1", a.port), partial(Handler, directory=str(folder)))
    print(f"pupil drawing page: http://127.0.0.1:{a.port}/  (answers -> {folder / 'marks.json'})", flush=True)
    srv.serve_forever()


if __name__ == "__main__":
    main()
