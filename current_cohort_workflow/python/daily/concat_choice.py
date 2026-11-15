"""Which concat_out folder of a session is THE table (the one Step 02 gets), for every tool the same way.

A session folder <pupil_root>/dayN/<cage>/<mouse>/ can hold:
  concat_out_<stamp>                  the pupil/TST package's table (pupil_align_ttl, original code)
  concat_out_<stamp>_syncaligned      fix_concat_alignment.py: frame<->pulse match from sync_sessions.py
  concat_out_<stamp>_fullpupil        fullpupil_sessions.py: pupil columns = the full pupil (whole dark area)
Corrected copies are OPTIONS for rare sessions (the user, 2026-11-15). A copy is used only while
  * it passed its tool's checks (fix_source.json present, combined_pupil_digital.csv present), and
  * its session is approved NOW in its list (alignment_fixes.csv for _syncaligned, fullpupil_sessions.csv for
    _fullpupil), and a _fullpupil copy built on a _syncaligned table also needs the alignment approval.
Withdrawing an approval (row removed, or 'approved' not yes / a date) therefore takes the session back to the
package's table. Approved = the 'approved' cell starts with 'yes' or with a date (YYYY-MM-DD).
Newest = latest stamp; for equal stamps the package's table < _syncaligned < _fullpupil.
"""
import csv
import json
import re
from pathlib import Path

SYNC, FULL = "_syncaligned", "_fullpupil"
DERIVED = (SYNC, FULL)
LISTS = {SYNC: "alignment_fixes.csv", FULL: "fullpupil_sessions.csv"}
_STAMP = re.compile(r"concat_out_(\d{8}_\d{6})")


def is_approved(v) -> bool:
    v = str(v or "").strip().lower()
    return v.startswith("yes") or bool(re.match(r"\d{4}-\d{2}-\d{2}", v))


def norm_key(day, cage, mouse):
    return int(str(day).strip()), str(cage).strip().lower(), str(mouse).strip().lower()


def approved_sessions(pupil_root: Path, list_name: str, warn=print) -> set:
    """{(day, cage, mouse)} approved in <pupil_root>/<list_name> (lower case). A malformed file gives an empty set and
    a warning, never an exception (it must not stop the grading of all other sessions)."""
    p = Path(pupil_root) / list_name
    out = set()
    if not p.exists():
        return out
    try:
        with open(p, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        for r in rows:
            r = {str(k).strip().lower(): v for k, v in r.items() if k is not None}
            d = str(r.get("day", "")).strip()
            if d.isdigit() and r.get("cage") and r.get("mouse") and is_approved(r.get("approved")):
                out.add(norm_key(d, r["cage"], r["mouse"]))
    except Exception as exc:   # noqa: BLE001
        if warn:
            warn(f"  !! {list_name} could not be read ({exc}); treated as empty")
    return out


def session_key(session_dir: Path):
    d = Path(session_dir)
    day = d.parts[-3]
    return norm_key(day[3:], d.parts[-2], d.parts[-1]) if day[3:].isdigit() else None


def _rank(p: Path):
    m = _STAMP.match(p.name)
    suf = 2 if p.name.endswith(FULL) else 1 if p.name.endswith(SYNC) else 0
    return (m.group(1) if m else p.name, suf, p.name)


def _source(p: Path) -> dict:
    try:
        return json.loads((p / "fix_source.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}


def usable_runs(session_dir: Path, pupil_root: Path, approvals=None):
    """usable concat_out folders of one session, oldest first (see the module docstring). approvals = optional
    {suffix: set} from approved_sessions, to avoid re-reading the lists for every session."""
    sd = Path(session_dir)
    if approvals is None:
        approvals = {s: approved_sessions(pupil_root, n) for s, n in LISTS.items()}
    key = session_key(sd)
    runs = sorted((p for p in sd.glob("concat_out_*") if p.is_dir() and (p / "combined_pupil_digital.csv").exists()),
                  key=_rank)
    out = []
    for p in runs:
        suf = next((s for s in DERIVED if p.name.endswith(s)), "")
        if not suf:
            out.append(p)
            continue
        if not (p / "fix_source.json").exists() or key not in approvals.get(suf, set()):
            continue
        if suf == FULL and str(_source(p).get("base_run", "")).endswith(SYNC) and key not in approvals.get(SYNC, set()):
            continue          # built on a sync-aligned table whose approval was withdrawn
        out.append(p)
    return out


def newest_usable(session_dir: Path, pupil_root: Path, approvals=None):
    r = usable_runs(session_dir, pupil_root, approvals)
    return r[-1] if r else None


def sha1(p: Path) -> str:
    import hashlib
    h = hashlib.sha1()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def fullpupil_trace_valid(session_dir: Path) -> bool:
    """the session's full-pupil trace (fullpupil_sessions.py) exists, was made from the package's CURRENT output
    (qc.json and CLEAN file unchanged since) and has not been changed or cut since it was made"""
    sd = Path(session_dir)
    od, g = sd / "pupil_outputs_fullpupil", sd / "pupil_outputs_general"
    src, clean = od / "fullpupil_source.json", od / "pupil_fullpupil_CLEAN.csv"
    if not (src.exists() and clean.exists() and (od / "fullpupil_qc.json").exists()):
        return False
    try:
        s = json.loads(src.read_text(encoding="utf-8"))
        i, o = s.get("inputs", {}), s.get("outputs", {})
        return (i.get("qc_json") == sha1(g / "qc.json") and
                i.get("package_clean") == sha1(g / "pupil_diameter_FINAL_U-Net_results_CLEAN.csv") and
                o.get("clean_sha1") == sha1(clean))
    except (OSError, ValueError):
        return False


def session_dir(pupil_root: Path, day: int, cage: str, mouse: str):
    """<pupil_root>/day<day>/<cage>/<mouse>, matching cage and mouse without regard to case; None if not there"""
    d = Path(pupil_root) / f"day{int(day)}"
    if not d.is_dir():
        return None
    c = next((p for p in d.iterdir() if p.is_dir() and p.name.lower() == str(cage).lower()), None)
    m = next((p for p in c.iterdir() if p.is_dir() and p.name.lower() == str(mouse).lower()), None) if c else None
    return m


def package_runs(session_dir: Path):
    return sorted((p for p in Path(session_dir).glob("concat_out_*")
                   if p.is_dir() and not p.name.endswith(DERIVED) and (p / "combined_pupil_digital.csv").exists()),
                  key=_rank)
