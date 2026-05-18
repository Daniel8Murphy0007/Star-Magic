"""
_emit_closure_json.py — one-line JSON-emit retrofit shim for session scripts.

Purpose
-------
Many `_session*.py` scripts print closures only to stdout and are flagged
PARSE_FAIL by `_uqff_program.py`. Importing this module and calling
`emit(label, predicted, observed)` (or `emit_many(...)`) writes a
canonical `_session{N}_closures.json` file in the repo root, which the
master runner's `_parse_session_json` walker will pick up automatically
and promote the script to status OK_JSON.

Auto-detects the session ID from the calling script's filename
(`_session{NNN}_*.py`). No CLI args, no fitting, no rounding — values
are stored at full float precision.

Usage (single closure)
----------------------
    from _emit_closure_json import emit
    emit("alpha_inv", predicted=137.0359990, observed=137.0359991)

Usage (multiple)
----------------
    from _emit_closure_json import emit_many
    emit_many([
        {"label": "M_pl_kg",       "predicted": 2.176e-8,    "observed": 2.176434e-8},
        {"label": "h_planck_Js",   "predicted": 6.62607e-34, "observed": 6.62607015e-34},
    ])

Output file
-----------
Written to: {repo_root}/_session{NNN}_closures.json
Schema (auto-detected by _uqff_program._parse_session_json):
    {
      "session": 305,
      "closures": [
        {
          "label":     "alpha_inv",
          "predicted": 137.0359990,
          "observed":  137.0359991,
          "residual_pct": 7.30e-7
        }
      ]
    }
"""

import json
import re
import sys
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parent
_SESSION_RE = re.compile(r"_session(\d+)_")


def _detect_session_id() -> str:
    """Find caller's session id from its filename. Walks the call stack."""
    depth = 1
    while True:
        try:
            frame = sys._getframe(depth)
        except ValueError:
            break
        fname = frame.f_globals.get("__file__", "") or ""
        m = _SESSION_RE.search(Path(fname).name)
        if m:
            return m.group(1)
        depth += 1
    # Last resort: scan sys.argv[0]
    m = _SESSION_RE.search(Path(sys.argv[0] if sys.argv else "").name)
    return m.group(1) if m else "NA"


def _residual_pct(predicted: float, observed: float) -> float:
    """Symmetric percent residual; full precision, no rounding."""
    if observed == 0:
        return float("inf") if predicted != 0 else 0.0
    return abs(predicted - observed) / abs(observed) * 100.0


def _emit_path(session_id: str) -> Path:
    return _REPO_ROOT / f"_session{session_id}_closures.json"


def emit(label: str, predicted: float, observed: float, **extras) -> Path:
    """Write a single-closure JSON file. Returns the file path."""
    return emit_many([{"label": label, "predicted": predicted,
                       "observed": observed, **extras}])


def emit_many(closures: list[dict], session_id: str | None = None) -> Path:
    """Write a multi-closure JSON file. Returns the file path.

    Each closure must contain `label`, `predicted`, `observed` (floats).
    All other keys are passed through. Residual % is auto-computed when
    not supplied. Full precision is preserved (no rounding)."""
    sid = session_id or _detect_session_id()
    payload = {"session": int(sid) if sid.isdigit() else sid, "closures": []}
    for c in closures:
        if "label" not in c or "predicted" not in c or "observed" not in c:
            raise ValueError("each closure needs label/predicted/observed")
        entry = dict(c)
        try:
            p = float(entry["predicted"]); o = float(entry["observed"])
        except (TypeError, ValueError):
            raise ValueError(f"predicted/observed not numeric for {entry.get('label')}")
        entry["predicted"] = p
        entry["observed"] = o
        if "residual_pct" not in entry:
            entry["residual_pct"] = _residual_pct(p, o)
        payload["closures"].append(entry)
    out = _emit_path(sid)
    out.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    return out


def append(label: str, predicted: float, observed: float, **extras) -> Path:
    """Append one closure to existing _session{N}_closures.json (creates if missing)."""
    sid = _detect_session_id()
    out = _emit_path(sid)
    if out.exists():
        try:
            existing = json.loads(out.read_text(encoding="utf-8"))
            closures = existing.get("closures") or []
        except Exception:
            closures = []
    else:
        closures = []
    closures.append({"label": label, "predicted": float(predicted),
                     "observed": float(observed),
                     "residual_pct": _residual_pct(float(predicted), float(observed)),
                     **extras})
    return emit_many(closures, session_id=sid)
