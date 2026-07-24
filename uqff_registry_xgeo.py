"""uqff_registry_xgeo — Cross-Geometry Derivation Campaign queue builder (XGEO).

Drains the R1 bulk-triage debt per Daniel's 2026-07-24 verdict:
  - Full derivation campaign over all 348 cross-geometry tasks
    (116 assimilation_dispatch observables x 3 non-owner geometries;
     the x3 numeric-mode axis is synthetic in _solve_via_dispatch and
     is collapsed here — one task per (observable, geometry) cell).
  - Terminal status for unrouted cells: CROSS_GEOMETRY_PENDING.

Inputs (read-only): assimilation_dispatch.DISPATCH, OVERDETERMINATION_MAP.csv
                    (protected baseline — never written).
Route rulings ledger: UNIFIED_REGISTRY_XGEO_ROUTES.csv (append-only; Daniel
  supplies route_formula + route_paper per cell across campaign sessions;
  this builder merges them on regeneration — Rule 10: Daniel provides the
  information, the session assembles it).

Emits (idempotent): UNIFIED_REGISTRY_XGEO_QUEUE.csv (348 rows)
"""
import csv
import os
import re
import sys

_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _DIR)
import assimilation_dispatch as _ad

GEOMS = ["qcalcgeom", "bsfg", "dpm", "d26"]
_PRIM_TOKENS = ["D_PHYS", "D_BSFG", "D_CRIT", "N_CH", "SO_5", "A_5",
                "F_TRZ", "PHI_RES", "SSQ", "BETA_I", "K_MEX", "S_26",
                "RHO_SCM", "RHO_UA", "OMEGA_SCM", "KAPPA", "PI"]

ROUTES_FILE = os.path.join(_DIR, "UNIFIED_REGISTRY_XGEO_ROUTES.csv")
QUEUE_FILE = os.path.join(_DIR, "UNIFIED_REGISTRY_XGEO_QUEUE.csv")


def _primitives_used(formula):
    hay = re.sub(r"[^A-Za-z0-9_]", "_", (formula or "")).upper()
    return ";".join(t for t in _PRIM_TOKENS if t in hay)


def _load_routes():
    routes = {}
    if os.path.exists(ROUTES_FILE):
        with open(ROUTES_FILE, encoding="utf-8", newline="") as f:
            for r in csv.DictReader(f):
                key = (r["observable"].strip(), r["target_geometry"].strip())
                routes[key] = r
    return routes


def build_queue():
    routes = _load_routes()
    rows = []
    for name in sorted(_ad.DISPATCH):
        rec = _ad.DISPATCH[name]
        owner = rec["owner_geometry"]
        for g in GEOMS:
            if g == owner:
                continue
            ruling = routes.get((name, g))
            if ruling and ruling.get("route_formula", "").strip():
                status = ruling.get("status", "").strip() or "XGEO_ROUTED"
                route_formula = ruling["route_formula"].strip()
                route_paper = ruling.get("route_paper", "").strip()
            else:
                status = "CROSS_GEOMETRY_PENDING"
                route_formula = ""
                route_paper = ""
            rows.append({
                "observable": name,
                "domain": rec["domain"],
                "owner_geometry": owner,
                "target_geometry": g,
                "owner_formula": rec["uqff_formula"],
                "owner_value": repr(rec["uqff_value"]),
                "target": "" if rec.get("target") is None else repr(rec["target"]),
                "primary_source": rec.get("primary_source", ""),
                "primitives_used": _primitives_used(rec["uqff_formula"]),
                "route_status": status,
                "route_formula": route_formula,
                "route_paper": route_paper,
            })
    return rows


def cell_status_map():
    """(observable, geometry) -> registry status, for registry_generator merge."""
    return {(r["observable"], r["target_geometry"]): r["route_status"]
            for r in build_queue()}


def main():
    rows = build_queue()
    cols = ["observable", "domain", "owner_geometry", "target_geometry",
            "owner_formula", "owner_value", "target", "primary_source",
            "primitives_used", "route_status", "route_formula", "route_paper"]
    with open(QUEUE_FILE, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    if not os.path.exists(ROUTES_FILE):
        with open(ROUTES_FILE, "w", encoding="utf-8", newline="") as f:
            csv.writer(f).writerow(["observable", "target_geometry",
                                    "route_formula", "route_paper", "status"])
    from collections import Counter
    st = Counter(r["route_status"] for r in rows)
    print(f"tasks={len(rows)} | " + ", ".join(f"{k}={v}" for k, v in sorted(st.items())))
    import hashlib
    print("queue_sha256=" + hashlib.sha256(open(QUEUE_FILE, "rb").read()).hexdigest()[:16])


if __name__ == "__main__":
    main()
