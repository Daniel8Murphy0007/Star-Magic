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
CONFIRMATIONS_FILE = os.path.join(_DIR, "UNIFIED_REGISTRY_XGEO_CONFIRMATIONS.csv")

# STANDING RULE (Daniel, 2026-07-24): LEADING-PRIMITIVE GEOMETRY RULE.
# A second (independent) route's geometry home is determined by the leading
# structural primitive of its dominant term:
#   A_5-led / rotational        -> dpm   (PAPER_2116 A_5 rotational lineage)
#   D_crit-led / 26-structure   -> d26
#   D_BSFG-led                  -> bsfg
#   F_U-machinery compositions on the {F_TRZ, SSq, Phi_res, beta_i} constant
#   set (F_TRZ ladders, kernel-K forms)  -> qcalcgeom (PAPER_1203/2124)
# Owners keep their corpus assignments; the rule classifies SECOND routes only.
# Promotions are appended to UNIFIED_REGISTRY_XGEO_ROUTES.csv with status
# XGEO_INDEPENDENT and win over identity routes on regeneration.

# XGEO-U: documented independent formula pairs for the same observable
# (corpus-sourced: dispatch S-variant families + notes 'Alternative' markers +
#  exact-name dispatch closure matches). All pairs share the OWNER geometry per
# the corpus's own assignments -> class SAME_GEOMETRY_CONFIRMATION (R1 category:
# alternates = numbered independent confirmations). Cell-level XGEO_INDEPENDENT
# upgrades require a corpus/Daniel geometry ruling per pair; none exists yet.
CONFIRMATION_PAIRS = [
    ("mp_me_ratio", "m_p/m_e = D_BSFG*pi^5", 6 * 3.141592653589793 ** 5,
     "A_5^2/2 + D_BSFG^2 = 1836 EXACT-integer", 60 ** 2 / 2 + 6 ** 2,
     1836.15267343, "S266;atlas-note"),
    ("alpha_s_M_Z", "1/(K_MEX*D_phys + F_TRZ)", 1.0 / ((25.0 / 12) * 4 + 0.1),
     "F_TRZ*K_MEX*SSQ - F_TRZ^3*Phi_res", 0.1 * (25.0 / 12) * 0.57 - 0.1 ** 3 * (5.0 / 6),
     0.1179, "S348;S378"),
    ("SM_higgs_lambda_S377", "F_TRZ*K_MEX*SSQ + F_TRZ^3*K_MEX*N_CH*SSQ",
     0.1 * (25.0 / 12) * 0.57 + 0.1 ** 3 * (25.0 / 12) * 9 * 0.57,
     "F_TRZ + F_TRZ^2*K_MEX + F_TRZ^2 - F_TRZ^3 (S441 variant)",
     0.1 + 0.01 * (25.0 / 12) + 0.01 - 0.001, 0.1293, "S377;S441"),
    ("SM_top_yukawa_S376", "1 - F_TRZ^2 (at m_t scale)", 1 - 0.01,
     "Phi_res + F_TRZ + F_TRZ^2 - F_TRZ^3*K_MEX (S440 variant; vs its OWN target 0.94 the residual is 0.13 pct — the 5.27 pct column value is vs S376's 0.9936 target, different scale)",
     5.0 / 6 + 0.1 + 0.01 - 0.001 * (25.0 / 12), 0.9936, "S376;S440"),
    ("hubble_tension", "H0 grid (2*SO_5+2)*F_TRZ^19 vs local anchor (PAPER_2093/2125)",
     22 * 0.1 ** 19, "1/12 EXACT tilt closure (PAPER_1156; dispatch p9_h_tension)",
     1.0 / 12, None, "PAPER_2093;PAPER_1156;PAPER_2125"),
]


def write_confirmations():
    rows = []
    for obs, f1, v1, f2, v2, target, src in CONFIRMATION_PAIRS:
        r1 = "" if target in (None, 0) else f"{abs(v1 - target) / abs(target) * 100:.4f}"
        r2 = "" if target in (None, 0) else f"{abs(v2 - target) / abs(target) * 100:.4f}"
        rows.append({"observable": obs, "route_1_formula": f1, "route_1_value": repr(v1),
                     "route_2_formula": f2, "route_2_value": repr(v2),
                     "target": "" if target is None else repr(target),
                     "route_1_residual_pct": r1, "route_2_residual_pct": r2,
                     "classification": "SAME_GEOMETRY_CONFIRMATION",
                     "sources": src})
    cols = ["observable", "route_1_formula", "route_1_value", "route_2_formula",
            "route_2_value", "target", "route_1_residual_pct", "route_2_residual_pct",
            "classification", "sources"]
    with open(CONFIRMATIONS_FILE, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    return rows


def confirmations_map():
    """observable -> confirmations text, for registry_generator merge."""
    return {obs: f"2 independent routes ({src}) — see UNIFIED_REGISTRY_XGEO_CONFIRMATIONS.csv"
            for obs, _f1, _v1, _f2, _v2, _t, src in CONFIRMATION_PAIRS}


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
    conf = write_confirmations()
    from collections import Counter
    st = Counter(r["route_status"] for r in rows)
    print(f"tasks={len(rows)} | " + ", ".join(f"{k}={v}" for k, v in sorted(st.items())))
    print(f"confirmations={len(conf)}")
    import hashlib
    for fn in (QUEUE_FILE, CONFIRMATIONS_FILE):
        print(os.path.basename(fn) + "_sha256=" +
              hashlib.sha256(open(fn, "rb").read()).hexdigest()[:16])


if __name__ == "__main__":
    main()
