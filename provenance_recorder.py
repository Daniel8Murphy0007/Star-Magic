"""
provenance_recorder.py
Part of UQFF Assimilation Geometry — Round 12 EXPANSION_PLAN deliverable.

Status: bundled in `uqff` PyPI package (Star-Magic repo) for academic peer-review
        access — supports the 3-university first peer review + the 8 NASA-Roses
        grant evaluation panels.

Future: relocatable to https://github.com/Daniel8Murphy0007/Aetheric-Propulsion
        for commercial-tier extraction. See ASSIMILATION_GEOMETRY_ATLAS.md for
        the extraction procedure.

Dependencies (internal):  none
Dependencies (external):  none

License:  AGPL-3.0-or-later OR Commercial (dual; identical to parent uqff)

----------------------------------------------------------------------------
PROVENANCE RECORDER
----------------------------------------------------------------------------
Builds the provenance_chain field that ships in every QCalcGeom.solve()
result dict. The chain is an ordered list of human-readable strings that
trace each derivation from the locked canonical primitives to the final
value, so any peer reviewer can audit the line of reasoning.

Format of each chain line:
    "<primitive_or_step> : <value_or_expression> [(<paper_citation>)]"

Examples:
    "rho_SCm : 7.09e-37 J/m^3 (PAPER_1271)"
    "K_MEX   : 25/12 (PAPER_1522 EXACT)"
    "Phi_res : 5/6 (PAPER_1159 G6 closure)"
    "poincare: 1/2 + F_TRZ * Phi_res_5_6 = 7/12 (PAPER_646)"
    "geometry: dpm (Di-Pseudo-Monopole 26-state mediator)"
    "numeric : symbolic (sympy exact Rational)"
"""

# Canonical citation strings for the 11 locked primitives + key derivative facts.
PRIMITIVE_CITATIONS = {
    "D_phys":      "PAPER_1167 integer primitive D_phys = 4",
    "D_BSFG":      "PAPER_1521 derivative: D_crit - 2*SO_5 = 6 EXACT",
    "D_crit":      "PAPER_1167 integer primitive D_crit = 26 (bosonic-string critical dim)",
    "N_CH":        "PAPER_1167 integer primitive N_CH = 9 (channel dimension)",
    "SO_5":        "PAPER_1167 integer primitive SO_5 = 10",
    "A_5":         "PAPER_1167 integer primitive A_5 = 60 (icosahedral group order)",
    "F_TRZ":       "PAPER_1160 G7 closure F_TRZ = 1/|SO(5)| = 1/10 EXACT",
    "Phi_res_5_6": "PAPER_1159 G6 closure Phi_res = [SSq]/Omega_Lambda = 5/6 EXACT",
    "K_MEX":       "PAPER_1522 derivative K_MEX = Phi_5/6 * SO_5 / D_phys = 25/12 EXACT",
    "SSQ":         "PAPER_1154 first-principles SSq = 0.57",
    "beta_i":      "PAPER_1203 canonical beta_i = 0.6029",
    "rho_SCm":     "PAPER_1271 canonical anchor rho_SCm = 7.09e-37 J/m^3",
    "rho_UA":      "PAPER_1167 structural rho_UA = 10 * rho_SCm",
    "T_10000":     "PAPER_1182 Riemann zero index 10000 along the critical line",
    "A_26":        "PAPER_1155 amplification A_26 = sum_{i=1..26} i^6 = 1,307,797,101 EXACT",
}

# Per-Millennium closure derivation chain elements (the structural reasoning).
CLOSURE_DERIVATION = {
    "yang_mills": {
        "geometry":  "bsfg",
        "formula":   "m_gap = 2 * D_phys * Lambda_QCD",
        "papers":    ["PAPER_1318"],
        "primitives": ["D_phys"],
    },
    "riemann": {
        "geometry":  "d26",
        "formula":   "t_10000 along the critical line; 26-fold radial derivative selects index 10000",
        "papers":    ["PAPER_1182", "PAPER_1080"],
        "primitives": ["D_crit", "T_10000"],
    },
    "navier_stokes": {
        "geometry":  "bsfg",
        "formula":   "enstrophy_cap = 1 - F_TRZ * D_BSFG / D_phys = 1 - (1/10)(6)/(4) = 17/20",
        "papers":    ["PAPER_1148"],
        "primitives": ["F_TRZ", "D_BSFG", "D_phys"],
    },
    "hodge": {
        "geometry":  "dpm",
        "formula":   "H_closure = (D_phys + D_BSFG) / SO_5 = (4+6)/10 = 1",
        "papers":    ["PAPER_1203_Nuclear", "PAPER_1230"],
        "primitives": ["D_phys", "D_BSFG", "SO_5"],
    },
    "poincare": {
        "geometry":  "dpm",
        "formula":   "P_closure = 1/2 + F_TRZ * Phi_res_5_6 = 1/2 + (1/10)(5/6) = 7/12",
        "papers":    ["PAPER_646"],
        "primitives": ["F_TRZ", "Phi_res_5_6"],
    },
    "p_vs_np": {
        "geometry":  "d26",
        "formula":   "P_closure = 1 - F_TRZ ** N_CH = 1 - (1/10)^9 = 0.999999999",
        "papers":    ["PAPER_1162"],
        "primitives": ["F_TRZ", "N_CH"],
    },
    "bsd": {
        "geometry":  "qcalcgeom",
        "formula":   "L'(E,1) via Cremona 37a1 anchors; UniversalGravity simultaneous solver",
        "papers":    ["PAPER_1149"],
        "primitives": [],
        "external_anchors": ["Cremona_37a1: Omega, R, Sha, c_p, |E_torsion|"],
    },
    "black_hole_info": {
        "geometry":  "qcalcgeom",
        "formula":   "Page recovery fraction from F_U=0 chain at M_BH = 1.989e30 kg",
        "papers":    ["PAPER_594", "PAPER_1213"],
        "primitives": [],
        "external_anchors": ["M_BH = 1.989e30 kg (M_Sun)"],
    },
}


def build_chain(observable, primary, alternate_paths=None):
    """Build a provenance_chain list for an observable.

    Args:
        observable: closure name
        primary: dict with at least {value, geometry, numeric, primary_source}
                 from the canonical owning geometry x numeric backend
        alternate_paths: optional dict from QCalcGeom.solve() for cross-reference

    Returns:
        list of strings; each is one step in the derivation chain
    """
    chain = []

    if observable not in CLOSURE_DERIVATION:
        chain.append(f"observable: {observable} (no canonical derivation registered)")
        chain.append(f"geometry  : {primary.get('geometry', '?')}")
        chain.append(f"numeric   : {primary.get('numeric', '?')}")
        if primary.get("value") is not None:
            chain.append(f"value     : {primary['value']}")
        return chain

    spec = CLOSURE_DERIVATION[observable]

    # Header: which geometry x numeric produced the canonical value
    chain.append(
        f"closure   : {observable} (Clay Millennium prize problem)"
    )
    chain.append(
        f"geometry  : {spec['geometry']} (canonical owner; "
        f"primary source: {','.join(spec['papers'])})"
    )
    chain.append(
        f"numeric   : {primary.get('numeric', '?')} backend"
    )

    # Show each primitive citation that feeds the formula
    for p in spec["primitives"]:
        cit = PRIMITIVE_CITATIONS.get(p, "(no citation registered)")
        chain.append(f"primitive : {p} -> {cit}")

    # External anchors (only for closures that depend on something outside UQFF)
    for anchor in spec.get("external_anchors", []):
        chain.append(f"anchor    : {anchor}")

    # The closure formula
    chain.append(f"formula   : {spec['formula']}")

    # Final value
    if primary.get("value") is not None:
        chain.append(f"value     : {primary['value']}")

    # Overdetermination summary, if alternate paths were collected
    if alternate_paths:
        N = sum(
            1
            for g in alternate_paths
            for n in alternate_paths[g]
            if alternate_paths[g][n].get("value") is not None
        )
        chain.append(f"overdet_N : {N} chains computed across the 4-geometry x 3-numeric matrix")

    return chain


__all__ = ["build_chain", "PRIMITIVE_CITATIONS", "CLOSURE_DERIVATION"]
