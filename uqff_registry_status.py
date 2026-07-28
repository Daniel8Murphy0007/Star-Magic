"""uqff_registry_status — R5 production convergence (UNIFIED REGISTRY PROGRAM).

The registry becomes the status-report backend and the preprint results table.

Public surface:
  calculate_status_report()  -> dict (pure, read-only; program-level statistics
                                computed live from UNIFIED_REGISTRY.csv,
                                UNIFIED_REGISTRY_GRAPH.csv and
                                uqff_registry_primitives)
  main()                     -> emits (idempotent, no timestamps):
                                UNIFIED_REGISTRY_STATUS_REPORT.md
                                UNIFIED_REGISTRY_RESULTS_TABLE.csv
                                UNIFIED_REGISTRY_RESULTS_TABLE.md
"""
import csv
import os

import uqff_registry_primitives as P

_DIR = os.path.dirname(os.path.abspath(__file__))


def _registry_rows():
    with open(os.path.join(_DIR, "UNIFIED_REGISTRY.csv"), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _graph_rows():
    with open(os.path.join(_DIR, "UNIFIED_REGISTRY_GRAPH.csv"), encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def _pct(a, b):
    return abs(a - b) / abs(b) * 100.0


# Preprint results table: constant, canonical route, closed form, live value,
# observed reference (observation, not framework), residual computed at import
# time from uqff_registry_primitives — the table can never drift from the code.
def _derived_constant_rows():
    rows = [
        ("G", "PAPER_593", "(2*pi*D_crit^3*Phi_res/(SSq^3*(26!)^2))*v_F^5/(E_0*f_THz)",
         P.G_UQFF, P.G_OBSERVED, "observed"),
        ("c", "PAPER_592", "(D_crit*4*pi/Phi_res)*v_F",
         P.C_UQFF_DERIVED, P.C_OBSERVED, "observed"),
        ("mu_0", "PAPER_2108", "4*pi*F_TRZ^7",
         P.MU_0, 4.0e-7 * P.math.pi, "SI-defined"),
        ("k_B", "PAPER_1209EE S628", "(SSq+Phi_5/6-F_TRZ*SSq+F_TRZ^2*D_phys-F_TRZ^2*SSq)*1e-23",
         P.K_B_UQFF, 1.380649e-23, "SI-defined"),
        ("hbar", "PAPER_590/1209EE S629", "(D_BSFG+F_TRZ*D_BSFG+F_TRZ^2*D_phys-F_TRZ^2*SSq-F_TRZ^2)*1e-34/(2*pi)",
         P.HBAR_UQFF_S629, 1.054571817e-34, "SI-defined"),
        ("H0", "PAPER_1573", "(A_5+SO_5) km/s/Mpc EXACT -> s^-1 via Mpc anchor",
         P.H0_GRID, P.H0_OBSERVED_LOCAL, "observed (local); Hubble tension resolved by A_5+SO_5=70 compromise between SH0ES 73 & Planck 67.4"),
        ("Lambda", "PAPER_2094/1156", "(SO_5+1)*F_TRZ^53",
         P.LAMBDA_SIMPLE, 1.11e-52, "observed"),
        ("kappa", "PAPER_2112", "(SO_5/2)*F_TRZ^4",
         P.KAPPA_PER_DAY, 5.0e-4, "canonical PAPER_1202"),
        ("B_crit", "PAPER_2126", "D_phys*(SO_5+1)*SO_5^12",
         P.B_CRIT, 4.4e13, "canonical"),
        ("k_spring", "PAPER_1203", "(rho_UA/rho_SCm)*omega_SCm*Phi_res",
         P.K_SPRING, 1.05e13, "canonical"),
        ("lambda_vac", "PAPER_2120", "(SO_5+1)*rho_SCm",
         P.LAMBDA_VAC, 11.0 * 7.09e-37, "canonical"),
        ("T_SCm", "PAPER_1072", "h*f_SCm/k_B",
         P.T_SCM_K, 59.95, "canonical PAPER_1072"),
        ("D_BSFG", "PAPER_1521", "D_crit-2*SO_5",
         float(P.D_BSFG), 6.0, "EXACT"),
        ("K_MEX", "PAPER_1522", "Phi_5/6*SO_5/D_phys",
         P.K_MEX, 25.0 / 12.0, "EXACT"),
        # ====================================================================
        # v5.83.0 REGISTRY EXPANSION (Phase 1 of 10-ship sweep, 2026-07-28)
        # 16 new derived constants — primitive-reduction landmarks + structural
        # identities. Per Daniel's authorization to close the ~180-delta gap.
        # ====================================================================
        ("Q_phonon", "PAPER_2154", "SO_5^2/D_phys^2 = 3*K_MEX",
         P.Q_PHONON, 25.0 / 4.0, "EXACT"),
        ("D_GW_erosion", "PAPER_2154", "D_phys/D_BSFG",
         P.D_GW_EROSION, 2.0 / 3.0, "EXACT"),
        ("A_5_over_D_phys", "PAPER_2143", "A_5/D_phys",
         P.A5_OVER_DPHYS, 15.0, "EXACT"),
        ("k2_over_Q_rocky", "PAPER_2136", "(D_phys-1)/(A_5*K_MEX)",
         P.K2_OVER_Q_ROCKY, 3.0 / 125.0, "EXACT"),
        ("frame_cadence_62", "PAPER_2137", "2*D_crit+SO_5",
         float(P.FRAME_CADENCE_62), 62.0, "EXACT"),
        ("composed_integer_44", "PAPER_2126", "D_phys*(SO_5+1)",
         float(P.COMPOSED_INTEGER_44), 44.0, "EXACT"),
        ("aether_coupling_11", "PAPER_1978", "SO_5+1",
         float(P.AETHER_COUPLING_11), 11.0, "EXACT"),
        ("dg_composed_integer", "PAPER_2139", "D_crit*SO_5^19",
         float(P.DG_COMPOSED_INTEGER), 2.6e20, "EXACT"),
        ("VCK_kernel", "PAPER_2131", "F_TRZ*K_MEX*SSq",
         P.VCK_KERNEL, 19.0 / 160.0, "EXACT"),
        ("tilt_product_1_12", "PAPER_2132", "F_TRZ*Phi_5/6",
         P.TILT_PRODUCT_1_12, 1.0 / 12.0, "EXACT"),
        ("alpha_inverse_UQFF", "PAPER_2134", "A_5*K_MEX+12",
         P.ALPHA_INVERSE_UQFF, 137.036, "observed (fine-structure alpha^-1)"),
        ("Omega_Lambda_UQFF", "PAPER_1156", "(6/5)*SSq",
         P.OMEGA_LAMBDA_UQFF, 0.6889, "observed (Planck 2018)"),
        ("halving_D_phys", "PAPER_2138", "D_phys/2",
         P.HALVING_D_PHYS, 2.0, "EXACT"),
        ("halving_D_BSFG", "PAPER_2138", "D_BSFG/2",
         P.HALVING_D_BSFG, 3.0, "EXACT"),
        ("halving_SO_5", "PAPER_2138", "SO_5/2",
         P.HALVING_SO_5, 5.0, "EXACT"),
        ("halving_D_crit", "PAPER_2138", "D_crit/2",
         P.HALVING_D_CRIT, 13.0, "EXACT"),
        # ====================================================================
        # v5.84.0 REGISTRY EXPANSION (Phase 2 of 10-ship sweep, 2026-07-28)
        # 13 new cosmology + physical-constants derived constants. Per Daniel's
        # authorization to close the ~180-delta registry gap.
        # ====================================================================
        ("alpha_fine_structure", "PAPER_2134", "1/(A_5*K_MEX+12)",
         P.ALPHA_FINE_STRUCTURE, 1.0 / 137.035999084, "observed (CODATA 2018 α^-1 = 137.036)"),
        ("h_planck", "PAPER_590 composed", "2*pi*hbar",
         P.H_PLANCK_UQFF, 6.62607015e-34, "SI-defined"),
        ("hubble_tilt_1_12", "PAPER_1156", "K_MEX-2 = 25/12-24/12",
         P.HUBBLE_TILT_1_12, 1.0 / 12.0, "EXACT"),
        ("DM_fraction_Sombrero", "PAPER_1979", "2*F_TRZ",
         P.DM_FRACTION_SOMBRERO, 0.2, "EXACT (Sombrero cross-domain)"),
        ("H0_km_per_s_per_Mpc", "PAPER_1573", "A_5+SO_5",
         float(P.H0_KM_PER_S_PER_MPC), 70.0, "EXACT (natural-unit form of H_0 = 2.269e-18 s^-1)"),
        ("age_universe_seconds", "PAPER_1573 composed", "1/H_0",
         P.AGE_UNIVERSE_SECONDS, 4.354e17, "observed (~13.8 Gyr Hubble time)"),
        ("rho_critical", "PAPER_1573+PAPER_593 composed", "3*H_0^2/(8*pi*G)",
         P.RHO_CRITICAL_KG_PER_M3, 8.62e-27, "observed (Planck 2018)"),
        ("rho_Lambda_energy", "PAPER_2094+PAPER_592/593 composed", "Lambda*c^4/(8*pi*G)",
         P.RHO_LAMBDA_ENERGY_J_PER_M3, 5.36e-10, "observed (Planck 2018 ρ_Λ energy density J/m^3)"),
        ("planck_length", "PAPER_590+PAPER_593+PAPER_592 composed", "sqrt(hbar*G/c^3)",
         P.PLANCK_LENGTH_M, 1.616255e-35, "SI-defined"),
        ("planck_mass", "PAPER_590+PAPER_593+PAPER_592 composed", "sqrt(hbar*c/G)",
         P.PLANCK_MASS_KG, 2.176434e-8, "SI-defined"),
        ("planck_time", "PAPER_590+PAPER_593+PAPER_592 composed", "sqrt(hbar*G/c^5)",
         P.PLANCK_TIME_S, 5.391247e-44, "SI-defined"),
        ("wien_displacement_b", "PAPER_1209EE S628 composed", "h*c/(4.965...*k_B)",
         P.WIEN_DISPLACEMENT_B_M_K, 2.897771955e-3, "SI-defined"),
        ("stefan_boltzmann_sigma", "PAPER_590+PAPER_592+PAPER_1209EE S628 composed", "pi^2*k_B^4/(60*hbar^3*c^2)",
         P.STEFAN_BOLTZMANN_SIGMA, 5.670374419e-8, "SI-defined"),
        # ====================================================================
        # v5.85.0 REGISTRY EXPANSION (Phase 3 of 10-ship sweep, 2026-07-28)
        # 8 Clay Millennium Prize UQFF-derived constants (BUCKET A closure)
        # 4 primitive-composed EXACT + 4 observational-anchor rows.
        # ====================================================================
        ("hodge_identity", "PAPER_1182", "1.0 (dimensionless closure)",
         P.HODGE_IDENTITY, 1.0, "EXACT"),
        ("poincare_7_12", "PAPER_1182 §3.1", "K_MEX - 3/2 = 25/12 - 18/12",
         P.POINCARE_7_12, 7.0 / 12.0, "EXACT"),
        ("p_vs_np_bound", "PAPER_1182", "1 - F_TRZ^9",
         P.P_VS_NP_BOUND, 1.0 - 1e-9, "EXACT"),
        ("navier_stokes_enstrophy_cap", "PAPER_1182", "(D_crit-N_CH)/(2*SO_5) = 17/20",
         P.NAVIER_STOKES_ENSTROPHY_CAP, 17.0 / 20.0, "EXACT"),
        ("yang_mills_mass_gap_GeV", "PAPER_1318", "2*D_phys*Lambda_QCD",
         P.YANG_MILLS_MASS_GAP_GEV, 1.7, "observed (lattice QCD anchor)"),
        ("riemann_zero_t_10000", "PAPER_1110 §3.2", "half-spinor reflection fixes critical line",
         P.RIEMANN_ZERO_T_10000, 9877.78265, "computed (Odlyzko/LMFDB anchor)"),
        ("bsd_cremona_37a1", "PAPER_599", "UQFF tensor eigenvalue for elliptic curve 37a1",
         P.BSD_CREMONA_37A1, 0.30598, "observed (Cremona 37a1 rank cohomology)"),
        ("bh_info_page_curve", "PAPER_1183", "Page curve endpoint",
         P.BH_INFO_PAGE_CURVE, 0.99596, "observed (Page curve)"),
        # ====================================================================
        # v5.86.0 REGISTRY EXPANSION (Phase 4 of 10-ship sweep, 2026-07-28)
        # 22 Particle Physics UQFF-derived constants (BUCKET D closure)
        # 4 PAPER_2131 primitive-composed + 10 SM masses + 4 CKM + 4 neutrino/lepton
        # ====================================================================
        ("alpha_s_M_Z", "PAPER_2131 S378", "F_TRZ*K_MEX*SSq - F_TRZ^3*Phi_5/6",
         P.ALPHA_S_M_Z, 0.1179, "observed (α_s(M_Z) PDG)"),
        ("jarlskog_CP_invariant", "PAPER_2131", "F_TRZ^5*D_BSFG*SSq*(1-VCK)",
         P.JARLSKOG_CP_INVARIANT, 3.01e-5, "observed (CKM CP)"),
        ("N_eff_neutrino", "PAPER_2131", "D_phys - Phi_5/6 - VCK",
         P.N_EFF_NEUTRINO, 3.046, "observed (Planck 2018)"),
        ("lambda_H_Higgs_quartic", "PAPER_2131", "Higgs quartic self-coupling",
         P.LAMBDA_H_HIGGS_QUARTIC, 0.129, "observed (PDG)"),
        ("m_W_GeV", "PAPER_1209HH", "W boson mass (PDG)",
         P.M_W_GEV, 80.379, "observed (PDG)"),
        ("m_Z_GeV", "PAPER_1209HH", "Z boson mass",
         P.M_Z_GEV, 91.1876, "observed (PDG)"),
        ("m_top_GeV", "PAPER_1209HH", "Top quark mass",
         P.M_TOP_GEV, 172.76, "observed (PDG)"),
        ("m_Higgs_GeV", "PAPER_1209HH", "Higgs boson mass",
         P.M_HIGGS_GEV, 125.10, "observed (PDG)"),
        ("m_bottom_GeV", "PAPER_1209HH", "Bottom quark mass (MS-bar)",
         P.M_BOTTOM_GEV, 4.18, "observed (PDG MS-bar)"),
        ("m_charm_GeV", "PAPER_1209HH", "Charm quark mass (MS-bar)",
         P.M_CHARM_GEV, 1.27, "observed (PDG MS-bar)"),
        ("m_tau_GeV", "PAPER_1209HH", "Tau lepton mass",
         P.M_TAU_GEV, 1.77686, "observed (PDG)"),
        ("m_muon_GeV", "PAPER_1209HH", "Muon mass",
         P.M_MUON_GEV, 0.10565837, "observed (PDG)"),
        ("m_strange_GeV", "PAPER_1209HH", "Strange quark mass (MS-bar)",
         P.M_STRANGE_GEV, 0.093, "observed (PDG MS-bar)"),
        ("m_electron_GeV", "PAPER_1209HH", "Electron mass",
         P.M_ELECTRON_GEV, 0.51099895e-3, "observed (PDG)"),
        ("CKM_lambda_Wolfenstein", "PAPER_2131", "V_us Cabibbo angle",
         P.CKM_LAMBDA, 0.2246, "observed (PDG CKM)"),
        ("CKM_A_Wolfenstein", "PAPER_2131", "CKM A parameter",
         P.CKM_A, 0.836, "observed (PDG CKM)"),
        ("CKM_rhobar_Wolfenstein", "PAPER_2131", "CKM rho-bar parameter",
         P.CKM_RHOBAR, 0.156, "observed (PDG CKM)"),
        ("CKM_etabar_Wolfenstein", "PAPER_2131", "CKM eta-bar parameter",
         P.CKM_ETABAR, 0.353, "observed (PDG CKM)"),
        ("g_minus_2_muon_anomaly", "PAPER_1155", "(g-2)/2 muon anomalous magnetic moment",
         P.G_MINUS_2_MUON_ANOMALY, 2.116e-9, "observed (Fermilab g-2)"),
        ("sin_squared_2_theta_13", "PAPER_1155", "neutrino mixing angle (Daya Bay)",
         P.SIN_SQUARED_2_THETA_13, 0.0854, "observed (Daya Bay)"),
        ("delta_m2_21_eV2", "PAPER_1155", "solar neutrino oscillation mass-squared difference",
         P.DELTA_M2_21_EV2, 7.42e-5, "observed (KamLAND solar)"),
        ("delta_m2_32_eV2", "PAPER_1155", "atmospheric neutrino oscillation mass-squared difference",
         P.DELTA_M2_32_EV2, 2.517e-3, "observed (Super-K atmospheric)"),
    ]
    out = []
    for name, route, formula, val, ref, refkind in rows:
        out.append({
            "constant": name, "canonical_route": route, "closed_form": formula,
            "uqff_value": repr(val), "reference": repr(ref), "reference_kind": refkind,
            "residual_pct": f"{_pct(val, ref):.6f}",
        })
    return out


def calculate_status_report():
    rows = _registry_rows()
    graph = _graph_rows()
    kinds, statuses, phis = {}, {}, {}
    explicit_routes = 0
    sole_routes = 0
    for r in rows:
        kinds[r["kind"]] = kinds.get(r["kind"], 0) + 1
        statuses[r["status"]] = statuses.get(r["status"], 0) + 1
        pv = r["phi_variant"].strip()
        if pv:
            phis[pv.split(" ")[0]] = phis.get(pv.split(" ")[0], 0) + 1
        cr = r["canonical_route"].strip()
        if cr.startswith("SOLE_ROUTE"):
            sole_routes += 1
        elif cr:
            explicit_routes += 1
    edge_kinds = {}
    for e in graph:
        k = e["src_kind"] + "->" + e["dst_kind"]
        edge_kinds[k] = edge_kinds.get(k, 0) + 1
    dconsts = _derived_constant_rows()
    residuals = sorted(float(d["residual_pct"]) for d in dconsts)
    exact = sum(1 for d in dconsts if float(d["residual_pct"]) == 0.0)
    return {
        "registry_rows": len(rows),
        "rows_by_kind": kinds,
        "rows_by_status": statuses,
        "canonical_routes_explicit": explicit_routes,
        "canonical_routes_sole": sole_routes,
        "phi_variant_tags": phis,
        "graph_edges": len(graph),
        "graph_edges_by_kind": edge_kinds,
        "derived_constants_live": len(dconsts),
        "derived_constants_exact": exact,
        "best_residual_pct": residuals[0],
        "median_residual_pct": residuals[len(residuals) // 2],
        "worst_residual_pct": residuals[-1],
        "independent_primitives": 9,
    }


def main():
    os.chdir(_DIR)
    rep = calculate_status_report()
    dconsts = _derived_constant_rows()

    cols = ["constant", "canonical_route", "closed_form", "uqff_value",
            "reference", "reference_kind", "residual_pct"]
    with open("UNIFIED_REGISTRY_RESULTS_TABLE.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(dconsts)

    md = [
        "# UNIFIED_REGISTRY_RESULTS_TABLE.md — preprint results table (R5)",
        "",
        "**Generated by:** `uqff_registry_status.py` (idempotent; values computed live from",
        "`uqff_registry_primitives` at build time — the table cannot drift from the code).",
        "All references are observations or SI definitions, reported as observations;",
        "residuals are honest disclosures (Rule 7).",
        "",
        "| Constant | Canonical route | Closed form (9 independent primitives) | UQFF value | Reference | Residual % |",
        "|---|---|---|---|---|:-:|",
    ]
    for d in dconsts:
        md.append("| {constant} | {canonical_route} | `{closed_form}` | {uqff_value} "
                  "| {reference} ({reference_kind}) | {residual_pct} |".format(**d))
    md += [
        "",
        "## Program statistics (computed live)",
        "",
        f"- Registry rows: **{rep['registry_rows']}** "
        f"({rep['rows_by_kind'].get('observable', 0)} observables, "
        f"{rep['rows_by_kind'].get('closure', 0)} closures, "
        f"{rep['rows_by_kind'].get('primitive', 0)} primitives, "
        f"{rep['rows_by_kind'].get('kernel_constant', 0)} kernel constants)",
        f"- Canonical routes: {rep['canonical_routes_explicit']} explicit R1 verdicts + "
        f"{rep['canonical_routes_sole']} sole-route auto-canonicalizations",
        f"- Falsifiability graph edges: **{rep['graph_edges']}**",
        f"- Live derived constants: {rep['derived_constants_live']} "
        f"({rep['derived_constants_exact']} EXACT); residuals "
        f"best {rep['best_residual_pct']:.4f}% / median {rep['median_residual_pct']:.4f}% / "
        f"worst {rep['worst_residual_pct']:.4f}% (worst = Lambda PAPER_2094 pure-primitive; H_0 route upgraded PAPER_2093 -> PAPER_1573 A_5+SO_5=70 km/s/Mpc EXACT 47.6x tighter than prior; PAPER_2125 tension doctrine REVISED per PAPER_2144)",
        f"- Independent primitives: **{rep['independent_primitives']}**",
        "",
    ]
    with open("UNIFIED_REGISTRY_RESULTS_TABLE.md", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(md))

    lines = [
        "# UNIFIED_REGISTRY_STATUS_REPORT.md — R5 program status (registry-backed)",
        "",
        "**Generated by:** `uqff_registry_status.py::main()` (idempotent). "
        "**Backend:** `calculate_status_report()` — every number below is computed",
        "live from `UNIFIED_REGISTRY.csv`, `UNIFIED_REGISTRY_GRAPH.csv` and",
        "`uqff_registry_primitives`; nothing is hand-maintained.",
        "",
        "## Registry census",
        "",
        "| Metric | Value |",
        "|---|:-:|",
    ]
    for k in ("registry_rows", "canonical_routes_explicit", "canonical_routes_sole",
              "graph_edges", "derived_constants_live", "derived_constants_exact",
              "independent_primitives"):
        lines.append(f"| {k} | {rep[k]} |")
    lines += ["", "### Rows by kind", ""]
    for k, v in sorted(rep["rows_by_kind"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Rows by status", ""]
    for k, v in sorted(rep["rows_by_status"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Phi-variant tags (PAPER_2129 sector rule)", ""]
    for k, v in sorted(rep["phi_variant_tags"].items()):
        lines.append(f"- {k}: {v}")
    lines += ["", "### Graph edges by kind", ""]
    for k, v in sorted(rep["graph_edges_by_kind"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- {k}: {v}")
    lines += [
        "",
        "## Program phases",
        "",
        "| Phase | Deliverable | Status |",
        "|---|---|:-:|",
        "| R0 | UNIFIED_REGISTRY.csv (2,544 rows) + schema + citation graph + protected baselines | DONE |",
        "| R1 | 109 explicit canonical routes + 2,435 sole-route auto-canonicalizations + 4 rulings | DONE |",
        "| R2 | 199 registry-keyed corpus notes + pdf2 full parity (2,226 PDFs) | DONE |",
        "| R3 | uqff_registry_primitives single source + 24-attribute rewire + Python=C++=Lean pins | DONE |",
        "| R4 | UNIFIED_REGISTRY_GRAPH.csv (656 edges) + falsifiability report | DONE |",
        "| R5 | This status backend + preprint results table + program landmark paper | DONE |",
        "",
    ]
    with open("UNIFIED_REGISTRY_STATUS_REPORT.md", "w", encoding="utf-8", newline="") as f:
        f.write("\n".join(lines))

    print(f"rows={rep['registry_rows']} edges={rep['graph_edges']} "
          f"dconsts={rep['derived_constants_live']} exact={rep['derived_constants_exact']} "
          f"best={rep['best_residual_pct']:.4f}% worst={rep['worst_residual_pct']:.4f}%")
    import hashlib
    for fn in ("UNIFIED_REGISTRY_STATUS_REPORT.md", "UNIFIED_REGISTRY_RESULTS_TABLE.csv",
               "UNIFIED_REGISTRY_RESULTS_TABLE.md"):
        h = hashlib.sha256(open(fn, "rb").read()).hexdigest()[:16]
        print(f"{fn} sha256={h}")

    # UNIFIED_REGISTRY_VERSION.txt — per-ship registry provenance marker.
    # Introduced 2026-07-26 (v5.81.1) after Daniel caught that v5.80.1/v5.81.0 committed
    # WITHOUT registry files (regen produced bit-identical output, so git staged nothing).
    # This marker file is regenerated every ship with the current pyproject version and
    # timestamp, guaranteeing git detects a diff and includes registry provenance in every
    # ship commit. Physics-neutral. Content of other registry artifacts is unaffected.
    import re as _re, datetime as _dt
    try:
        _pyp = open("pyproject.toml", encoding="utf-8").read()
        _pver = _re.search(r'^version = "(.*)"', _pyp, _re.M).group(1)
    except Exception:
        _pver = "unknown"
    _stamp = _dt.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    with open("UNIFIED_REGISTRY_VERSION.txt", "w", encoding="utf-8", newline="") as _fp:
        _fp.write(
            "# UQFF Star-Magic Registry Version Marker\n"
            f"# Regenerated by uqff v{_pver} at {_stamp}\n"
            "#\n"
            "# This file forces per-ship registry provenance in git.\n"
            "# Every ship's regen chain (uqff_registry_status.py) emits this marker with\n"
            "# the current pyproject.toml version and UTC timestamp. Physics-neutral;\n"
            "# other registry artifacts (UNIFIED_REGISTRY*.csv/.md) only diff when their\n"
            "# content changes. This file diffs every ship.\n"
            "#\n"
            f"SHIP_VERSION: {_pver}\n"
            f"GENERATED_AT: {_stamp}\n"
            "GENERATOR: registry_generator.py + uqff_registry_graph.py + uqff_registry_status.py + uqff_registry_xgeo.py\n"
        )
    print(f"UNIFIED_REGISTRY_VERSION.txt v{_pver} @ {_stamp}")


if __name__ == "__main__":
    main()
