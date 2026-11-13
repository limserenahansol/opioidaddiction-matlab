"""Local paths for the one-click pupil scripts.

Paths are NOT committed (repo rule). Copy one_click_paths.example.json to one_click_paths.json (git-ignored) in this folder
and fill in your folders. Every script also accepts the paths as command-line arguments.
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def load():
    p = HERE / "one_click_paths.json"
    if not p.exists():
        sys.exit(f"Missing {p}\nCopy one_click_paths.example.json to one_click_paths.json and set your folders.")
    cfg = json.loads(p.read_text(encoding="utf-8"))
    bad = [k for k, v in cfg.items() if isinstance(v, str) and v.startswith("<")]
    if bad:
        sys.exit(f"Set these entries in {p}: {bad}")
    return cfg


def get(key, argv_index=None):
    """Command-line argument argv_index if given, else the config entry."""
    if argv_index is not None and len(sys.argv) > argv_index:
        return Path(sys.argv[argv_index])
    return Path(load()[key])
