"""UQFF Star-Magic REST API.

A small FastAPI server exposing UQFF closures over HTTP.

Install + run:
    pip install 'uqff[api]'
    uqff serve                          # starts on http://localhost:8000

Or directly:
    pip install fastapi 'uvicorn[standard]'
    python -m uqff_api                  # same as uqff serve
    # OR
    uvicorn uqff_api:app --reload

Endpoints (auto-documented at http://localhost:8000/docs):
    GET  /                              landing JSON
    GET  /version                       version + metrics
    GET  /status                        production status summary
    GET  /surfaces                      list 34 public calculate_* surfaces
    GET  /list?filter=X&all=true        list closure names
    GET  /search?q=X                    multi-namespace search
    GET  /predict/{name}                fetch a closure by name
    GET  /atlas                         return CLOSURE_ATLAS overview
    GET  /healthz                       health check (always 200 if running)
"""
import sys
from typing import Optional

try:
    from fastapi import FastAPI, HTTPException, Query
    from fastapi.responses import JSONResponse, HTMLResponse
except ImportError as e:  # noqa
    raise SystemExit(
        "FastAPI not installed. Install with: pip install 'uqff[api]'"
    ) from e

import uqff_pure_calculator as u

# Reuse the CLI's multi-namespace search helpers
try:
    from uqff_cli import (
        _all_paradox_keys, _all_millennium_keys, _all_lenr_keys,
        _all_nuclear_keys, _all_bucket_observables,
        _try_paradox, _try_lenr_full, _try_nuclear, _try_bucket_observable,
        _VERSION,
    )
except ImportError:
    _VERSION = "5.48.1"
    def _all_paradox_keys(): return sorted(u.PARADOX_TO_CLOSURE.keys())
    def _all_millennium_keys(): return sorted(getattr(u, "PARADOX_TO_MILLENNIUM", {}).keys())
    def _all_lenr_keys():
        try: return sorted(u.calculate_lenr_full({})["value"].keys())
        except Exception: return []
    def _all_nuclear_keys():
        try: return sorted(u.calculate_nuclear_magic({})["value"].keys())
        except Exception: return []
    def _all_bucket_observables():
        return {}
    def _try_paradox(name): return u.calculate_paradox({"paradox": name}).get("value")
    def _try_lenr_full(name): return None
    def _try_nuclear(name): return None
    def _try_bucket_observable(name): return None


app = FastAPI(
    title="UQFF Star-Magic REST API",
    description=(
        "HTTP interface for the UQFF Star-Magic physics framework. "
        "9 truly-independent primitives; 794 paradox closures + 8 Clay Millennium "
        "prize problems + complete SM 12-fermion spectrum + 18 ΛCDM observables. "
        "Auto-generated Swagger UI at /docs. License: AGPL-3.0 + commercial. "
        "See https://github.com/Daniel8Murphy0007/Star-Magic"
    ),
    version=_VERSION,
    license_info={
        "name": "AGPL-3.0 + commercial dual license",
        "url": "https://github.com/Daniel8Murphy0007/Star-Magic/blob/main/LICENSE",
    },
)


@app.get("/", tags=["info"])
def root():
    """Landing page metadata."""
    return {
        "name": "UQFF Star-Magic REST API",
        "version": _VERSION,
        "docs": "/docs",
        "openapi": "/openapi.json",
        "endpoints": [
            "/version", "/status", "/surfaces", "/list", "/search",
            "/predict/{name}", "/atlas", "/healthz",
        ],
        "github": "https://github.com/Daniel8Murphy0007/Star-Magic",
        "license": "AGPL-3.0 + commercial",
    }


@app.get("/healthz", tags=["info"])
def healthz():
    """Always 200 if the server is running."""
    return {"status": "ok"}


@app.get("/version", tags=["info"])
def version():
    """Version + headline metrics."""
    summary = {}
    try:
        summary = u.calculate_status_report({})["value"]["summary"]
    except Exception:
        pass
    return {
        "uqff_version": _VERSION,
        "python_version": sys.version.split()[0],
        "closures_paradox_dispatch": summary.get("total_closures"),
        "truly_independent_primitives": summary.get("truly_independent_primitives"),
        "derivative_primitives": summary.get("derivative_primitives"),
    }


@app.get("/status", tags=["info"])
def status():
    """Full production status summary."""
    try:
        return u.calculate_status_report({})["value"]
    except Exception as e:  # noqa
        raise HTTPException(status_code=500, detail=f"status_report failed: {e}")


@app.get("/surfaces", tags=["discovery"])
def surfaces():
    """List the 34 public calculate_* surfaces."""
    publics = sorted(n for n in dir(u) if n.startswith("calculate_"))
    return {"surfaces": publics, "count": len(publics)}


@app.get("/list", tags=["discovery"])
def list_keys(
    filter: Optional[str] = Query(None, description="substring filter (case-insensitive)"),
    all: bool = Query(False, description="include all 5 namespaces (default: paradox only)"),
):
    """List closure names. Use ?all=true to include all 5 namespaces."""
    groups: dict[str, list[str]] = {}
    if all:
        groups["PARADOX_TO_CLOSURE"] = _all_paradox_keys()
        groups["PARADOX_TO_MILLENNIUM"] = _all_millennium_keys()
        groups["calculate_lenr_full"] = _all_lenr_keys()
        groups["calculate_nuclear_magic"] = _all_nuclear_keys()
        groups.update(_all_bucket_observables())
    else:
        groups["PARADOX_TO_CLOSURE"] = _all_paradox_keys()

    if filter:
        flt = filter.lower()
        for k in list(groups.keys()):
            groups[k] = [n for n in groups[k] if flt in n.lower()]
            if not groups[k]:
                del groups[k]

    total = sum(len(v) for v in groups.values())
    return {"total": total, "namespaces": len(groups), "by_source": groups}


@app.get("/search", tags=["discovery"])
def search(q: str = Query(..., description="search substring (case-insensitive)")):
    """Search all 5 closure namespaces for matches."""
    needle = q.lower().strip()
    if not needle:
        raise HTTPException(status_code=400, detail="query string ?q= is required")

    results: dict[str, list[str]] = {}
    for source, keys in [
        ("PARADOX_TO_CLOSURE", _all_paradox_keys()),
        ("PARADOX_TO_MILLENNIUM", _all_millennium_keys()),
        ("calculate_lenr_full", _all_lenr_keys()),
        ("calculate_nuclear_magic", _all_nuclear_keys()),
    ]:
        hits = [k for k in keys if needle in k.lower()]
        if hits:
            results[source] = hits

    bucket = _all_bucket_observables()
    for surf, names in bucket.items():
        hits = [n for n in names if needle in n.lower()]
        if hits:
            results[surf] = hits

    total = sum(len(v) for v in results.values())
    return {"substring": needle, "total_matches": total, "by_source": results}


@app.get("/predict/{name}", tags=["lookup"])
def predict(name: str):
    """Fetch a closure value by name. Searches all 5 namespaces in order."""
    lookup_name = name.lower().strip()
    for source_name, fn in [
        ("PARADOX_TO_CLOSURE", _try_paradox),
        ("calculate_lenr_full", _try_lenr_full),
        ("calculate_nuclear_magic", _try_nuclear),
        ("bucket_observables", _try_bucket_observable),
    ]:
        try:
            value = fn(lookup_name)
        except Exception:
            value = None
        if value is not None:
            return {"name": name, "source": source_name, "value": value}

    raise HTTPException(
        status_code=404,
        detail=f"Closure '{name}' not found in any of the 5 UQFF namespaces. "
               f"Use /search?q={name} to find candidates."
    )


@app.get("/atlas", tags=["info"])
def atlas():
    """Return the CLOSURE_ATLAS summary numbers (counts only)."""
    paradox_n = len(u.PARADOX_TO_CLOSURE)
    mill_n = len(getattr(u, "PARADOX_TO_MILLENNIUM", {}))
    mill_derive = len(getattr(u, "_MILLENNIUM_DERIVE", {}))
    publics_n = sum(1 for n in dir(u) if n.startswith("calculate_"))
    try:
        lenr_n = len(u.calculate_lenr_full({})["value"])
    except Exception:
        lenr_n = 0
    try:
        nuclear_n = len(u.calculate_nuclear_magic({})["value"])
    except Exception:
        nuclear_n = 0
    bucket = _all_bucket_observables()
    total_bucket_obs = sum(len(v) for v in bucket.values())
    return {
        "paradox_dispatch_keys": paradox_n,
        "millennium_aliases": mill_n,
        "millennium_derivation_functions": mill_derive,
        "calculate_lenr_full_subkeys": lenr_n,
        "calculate_nuclear_magic_subkeys": nuclear_n,
        "bucket_observables_total": total_bucket_obs,
        "bucket_observables_per_surface": {k: len(v) for k, v in bucket.items()},
        "public_calculate_surfaces": publics_n,
    }


def run(host: str = "127.0.0.1", port: int = 8000, reload: bool = False):
    """Launch the API via uvicorn."""
    try:
        import uvicorn
    except ImportError as e:
        raise SystemExit(
            "uvicorn not installed. Install with: pip install 'uqff[api]'"
        ) from e
    print(f"\nUQFF Star-Magic REST API v{_VERSION}")
    print(f"  Listening on http://{host}:{port}")
    print(f"  Swagger UI:  http://{host}:{port}/docs")
    print(f"  OpenAPI:     http://{host}:{port}/openapi.json")
    print("  Stop:        Ctrl+C")
    print()
    uvicorn.run("uqff_api:app", host=host, port=port, reload=reload)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Launch the UQFF REST API")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true", help="auto-reload on code change")
    args = parser.parse_args()
    run(host=args.host, port=args.port, reload=args.reload)
