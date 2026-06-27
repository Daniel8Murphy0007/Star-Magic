# UQFF Star-Magic v5.28.0 — REST API + Jupyter integration

**Date:** 2026-06-22
**Type:** Minor release (new features, fully backward-compatible)
**Upgrade:** `pip install --upgrade 'uqff[api,jupyter]'`

This is the first **minor** version bump since v5.27.0 (matching semantic
versioning: new features → minor bump). No physics changes; this release
ships two new ways to talk to UQFF.

---

## 🌐 REST API (FastAPI)

```bash
pip install 'uqff[api]'
uqff serve                            # http://localhost:8000
# or:
uvicorn uqff_api:app --reload
```

Auto-generated Swagger UI at `http://localhost:8000/docs`.

### Endpoints

```
GET  /                       landing JSON
GET  /version                version + headline metrics
GET  /status                 production status summary
GET  /surfaces               34 public calculate_* surfaces
GET  /list?filter=X&all=true list closure names
GET  /search?q=X             multi-namespace substring search
GET  /predict/{name}         fetch a closure (searches all 5 namespaces)
GET  /atlas                  CLOSURE_ATLAS summary counts
GET  /healthz                health probe (always 200 if running)
```

Example:

```bash
curl http://localhost:8000/predict/holmlid_D_minus_1
# {"name":"holmlid_D_minus_1","source":"calculate_lenr_full",
#  "value":{"holmlid_D_minus_1":{"KER_eV":630.0,...,"reference":"PAPER_1133"}}}

curl 'http://localhost:8000/search?q=millennium'
# {...all 8 Clay Millennium aliases...}
```

### Why this matters

- **Frontend integrations**: any web or mobile app can now call UQFF over HTTP
- **Multi-language support**: anyone can use UQFF from JavaScript, Go, Rust, etc.
- **Demo deployment**: free-tier hosts (Fly.io, Railway, Cloudflare Workers) can run the API as an always-on public demo
- **Documentation as code**: the Swagger UI at `/docs` is auto-generated from the FastAPI route definitions; lets reviewers explore the API interactively

---

## 📓 Jupyter / IPython integration

```bash
pip install 'uqff[jupyter]'
```

In a Jupyter notebook:

```python
import uqff_jupyter
uqff_jupyter.enable()                  # activate rich HTML rendering

import uqff_pure_calculator as u
u.calculate_status_report({})          # now renders as styled HTML table!
```

Load the `%uqff` line magic for inline lookups:

```python
%load_ext uqff_jupyter
%uqff predict hubble_tension           # styled card with the closure value
%uqff search holmlid                   # multi-namespace search result
%uqff status                           # production summary
```

### Why this matters

- Any UQFF dict result auto-renders as a colour-coded HTML table
- `residual_pct` columns color-code red (>1%) vs. green (<1%)
- Observable lists render as proper tables (first 50 rows)
- `%uqff` magic gives REPL-style quick lookups

---

## 🆕 `uqff serve` CLI subcommand

```
uqff serve [--host 127.0.0.1] [--port 8000] [--reload]
```

Convenience wrapper around `uvicorn`. Same behaviour as
`uvicorn uqff_api:app` but with sensible defaults.

---

## Installation options

```bash
pip install uqff                       # CLI + calculator only
pip install 'uqff[api]'                # + FastAPI REST server
pip install 'uqff[jupyter]'            # + IPython rich display
pip install 'uqff[docs]'               # + Sphinx
pip install 'uqff[test]'               # + coverage / pytest
pip install 'uqff[all]'                # everything
```

---

## What didn't change

- 9 truly-independent primitives, 794 closures, 857/857 fidelity gate, ΔBIC = 94.1
- Physics + CLI surface from v5.27.2 unchanged
- License unchanged (AGPL-3.0 + commercial)
- `uqff_pure_calculator.py` untouched

---

## Files added in v5.28.0

| File | Purpose | Size |
|---|---|---|
| `uqff_api.py` | FastAPI server with 9 endpoints | 7 KB |
| `uqff_jupyter.py` | IPython rich-display + `%uqff` magic | 9 KB |
| `RELEASE_NOTES_v5.28.0.md` | This document | 3 KB |
| `notebooks/04_repl_jupyter_demo.ipynb` | Demo notebook | — |

---

## Quick verification after upgrading

```bash
pip install --upgrade 'uqff[api]'
uqff version                                # → uqff 5.28.0
uqff serve --port 8000 &                    # launch server in background
curl http://localhost:8000/version          # → JSON metrics
curl http://localhost:8000/predict/yang_mills   # → 1.736
kill %1                                     # stop server
```

---

**Author:** Daniel T. Murphy / Star-Magic Research Program
**Copyright:** © 2025-2026, all rights reserved.
