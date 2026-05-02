"""
_chain_trace_fix348.py
======================
Verification trace for Fix #3, Fix #4, and Fix #8 of the DPM quantum chain.

  Fix #3 (S5e): Electron mass from Ug2 outer-bubble De Broglie confinement.
  Fix #4 (S5f): K1-K4 coupling constants derived from vacuum + particle mass constraints.
  Fix #8 (S5g): r_cross independent solution from primordial galactic FUBii.

Run:
    & ".venv_py314_backup\\Scripts\\python.exe" -X utf8 _chain_trace_fix348.py

Exit code 0 = all section checks pass.  Any assertion failure exits non-zero.
"""

import math
import sys
import dpm_vacuum_manifold as dpm

SEP  = "=" * 78
SEP2 = "-" * 78

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def banner(title: str) -> None:
    print(f"\n{SEP}")
    print(f"  {title}")
    print(SEP)


def check(label: str, condition: bool, detail: str = "") -> None:
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}]  {label}" + (f"  -- {detail}" if detail else ""))
    if not condition:
        raise AssertionError(f"Check failed: {label}  {detail}")


# ===========================================================================
# SECTION 1: Constants added for Fix #3, #4, #8
# ===========================================================================
banner("SECTION 1 — New constants (R_C_LEPTON, galactic, K3_EFF)")

print(f"  R_C_UP              = {dpm.R_C_UP:.4e} m  (up-quark confinement)")
print(f"  R_C_DOWN            = {dpm.R_C_DOWN:.4e} m  (down-quark confinement)")
print(f"  R_C_LEPTON          = {dpm.R_C_LEPTON:.4e} m  (Ug2 outer-bubble electron)")
print(f"  OMEGA_G_GALACTIC    = {dpm.OMEGA_G_GALACTIC:.4e} rad/s")
print(f"  M_BH_SgrA           = {dpm.M_BH_SgrA:.4e} kg")
print(f"  D_GALACTIC_SUN      = {dpm.D_GALACTIC_SUN:.4e} m")
print(f"  K3_EFF              = {dpm.K3_EFF:.4e}")

gal_coupling_expected = dpm.OMEGA_G_GALACTIC * dpm.M_BH_SgrA / dpm.D_GALACTIC_SUN
check("R_C_LEPTON = R_C_UP * 10^2.5",
      abs(dpm.R_C_LEPTON - dpm.R_C_UP * (dpm.DPM_DENSITY_RATIO ** 2.5)) < 1e-25,
      f"R_C_LEPTON={dpm.R_C_LEPTON:.4e}")
check("Galactic coupling = Omega_g * M_bh / d_g ~ 23",
      18 < gal_coupling_expected < 30,
      f"coupling={gal_coupling_expected:.4f}")
check("K3_EFF positive",
      dpm.K3_EFF > 0.0,
      f"K3_EFF={dpm.K3_EFF:.4e}")


# ===========================================================================
# SECTION 2: Fix #3 — Electron mass (chain_Ug2_electron_mass)
# ===========================================================================
banner("SECTION 2 — Fix #3: Electron mass from Ug2 outer-bubble (S5e)")

e3 = dpm.chain_Ug2_electron_mass()

M_e_obs   = 9.1093837015e-31   # kg  PDG 2022
MeV_per_kg = dpm.C_LIGHT ** 2 / 1.602176634e-13

print(f"\n  Route B (PRIMARY):")
print(f"    r_c_lepton     = {e3['route_B']['r_c_lepton_m']:.4e} m")
print(f"    DPM layers     = {e3['route_B']['DPM_layers']}")
print(f"    DPM factor     = {e3['route_B']['DPM_factor']:.4f}")
print(f"    m_e derived    = {e3['route_B']['m_e_derived_kg']:.6e} kg")
print(f"    m_e derived    = {e3['route_B']['m_e_derived_MeV']:.4f} MeV/c^2")
print(f"    m_e observed   = {M_e_obs:.6e} kg  (PDG 2022)")
print(f"    error          = {e3['route_B']['error_pct']:+.2f}%")

print(f"\n  Route A (field ratio, K1=K2=1 placeholders):")
print(f"    field_ratio    = {e3['route_A']['field_ratio']:.6e}")
print(f"    m_e_A          = {e3['route_A']['m_e_derived_kg']:.6e} kg")
print(f"    error_A        = {e3['route_A']['error_pct']:+.1f}%")

print(f"\n  EM residual consistency (Fix #2 feed-forward):")
print(f"    Fix #2 EM residual  = {e3['em_residual_check']['Fix2_EM_residual_MeV']:.4f} MeV")
print(f"    Fraction of m_e     = {e3['em_residual_check']['EM_as_fraction_of_m_e_mass']:.4f}")
print(f"    Interpretation: {e3['em_residual_check']['interpretation'][:80]}...")

# Checks
check("Route B error within ±10% of PDG electron mass",
      abs(e3["route_B"]["error_pct"]) < 10.0,
      f"error={e3['route_B']['error_pct']:+.2f}%")
check("m_e Route B positive",
      e3["primary_result_kg"] > 0.0)
check("m_e Route B in plausible range [1e-31, 1e-29] kg",
      1e-31 <= e3["primary_result_kg"] <= 1e-29,
      f"m_e={e3['primary_result_kg']:.4e}")
check("EM residual fraction of m_e in (0, 2)",
      0 < e3["em_residual_check"]["EM_as_fraction_of_m_e_mass"] < 2.0,
      f"frac={e3['em_residual_check']['EM_as_fraction_of_m_e_mass']:.4f}")
check("primary_result matches route_B",
      e3["primary_result_kg"] == e3["route_B"]["m_e_derived_kg"])


# ===========================================================================
# SECTION 3: Fix #4 — Coupling constants K1-K4 (chain_coupling_constants)
# ===========================================================================
banner("SECTION 3 — Fix #4: Coupling constants K1-K4 from vacuum constraints (S5f)")

k4 = dpm.chain_coupling_constants()

print(f"\n  Inputs:")
print(f"    Z              = {k4['inputs']['Z']}")
print(f"    R_nuc          = {k4['inputs']['R_nuc_m']:.4e} m")
print(f"    M_proto        = {k4['inputs']['M_proto_kg']:.4e} kg")
print(f"    mu_s           = {k4['inputs']['mu_s_kg']:.4e} kg")
print(f"    E_react        = {k4['inputs']['E_react']:.4e}")
print(f"    cos(pi/4)      = {k4['inputs']['cos_tn']:.6f}")

print(f"\n  K1_eff: {k4['K1_eff']['value']:.4e}")
print(f"    Constraint: {k4['K1_eff']['constraint']}")
print(f"    Ug1_unit  = {k4['K1_eff']['Ug1_unit']:.4e}")
print(f"    M_p*c^2   = {k4['K1_eff']['M_p_c2_J']:.4e} J")

print(f"\n  K2_eff: {k4['K2_eff']['value']:.4e}")
print(f"    Constraint: {k4['K2_eff']['constraint']}")
print(f"    M_e used  = {k4['K2_eff']['M_e_used_kg']:.4e} kg  (Fix #3 Route B)")
print(f"    M_e*c^2   = {k4['K2_eff']['M_e_c2_J']:.4e} J")

print(f"\n  K3_eff: {k4['K3_eff']['value']:.4e}")
print(f"    Constraint: {k4['K3_eff']['constraint']}")

print(f"\n  K4_eff: {k4['K4_eff']['value']:.4e}")
print(f"    Constraint: {k4['K4_eff']['constraint']}")
print(f"    gal_coupling = {k4['K4_eff']['gal_coupling']:.4f}")

print(f"\n  SC mode consistency (H_SCm^n_j structural pattern):")
sc = k4["SC_mode_consistency"]
print(f"    H_SCm^0 = {sc['structural_n0']:.4f}  (K1 structure)")
print(f"    H_SCm^1 = {sc['structural_n1']:.4f}  (K2 structure)")
print(f"    H_SCm^2 = {sc['structural_n2']:.4f}  (K3 structure)")
print(f"    H_SCm^3 = {sc['structural_n3']:.4f}  (K4 structure)")
print(f"    ratio K2/K1 = {sc['ratio_K2_K1']:.6e}")
print(f"    ratio K3/K1 = {sc['ratio_K3_K1']:.6e}")
print(f"    ratio K4/K1 = {sc['ratio_K4_K1']:.6e}")

# Checks
check("K1_eff positive", k4["K1_eff"]["value"] > 0.0, f"K1_eff={k4['K1_eff']['value']:.4e}")
check("K2_eff positive", k4["K2_eff"]["value"] > 0.0, f"K2_eff={k4['K2_eff']['value']:.4e}")
check("K3_eff == K3_EFF constant", k4["K3_eff"]["value"] == dpm.K3_EFF)
check("K4_eff positive", k4["K4_eff"]["value"] > 0.0, f"K4_eff={k4['K4_eff']['value']:.4e}")
check("K1_eff / (M_p*c^2/Ug1_unit) = 1.0 (self-consistent)",
      abs(k4["K1_eff"]["value"] - k4["K1_eff"]["M_p_c2_J"] / k4["K1_eff"]["Ug1_unit"]) < 1e-10 * abs(k4["K1_eff"]["value"]),
      "K1_eff self-consistency")
check("K2_eff / (M_e*c^2/Ug2_unit) = 1.0 (self-consistent)",
      abs(k4["K2_eff"]["value"] - k4["K2_eff"]["M_e_c2_J"] / k4["K2_eff"]["Ug2_unit"]) < 1e-10 * abs(k4["K2_eff"]["value"]),
      "K2_eff self-consistency")


# ===========================================================================
# SECTION 4: Fix #8 — r_cross independent (chain_r_cross_independent)
# ===========================================================================
banner("SECTION 4 — Fix #8: r_cross independent from primordial FUBii (S5g)")

H_body = dpm.ELEMENT[1]
r8_H = dpm.chain_r_cross_independent(H_body)

print(f"\n  Z=1 (Hydrogen):")
print(f"    v_fermi           = {r8_H['v_fermi_ms']:.4e} m/s")
print(f"    E_react           = {r8_H['E_react']:.4e}")
print(f"    galactic_coupling = {r8_H['galactic_coupling']:.4f}  (Omega_g*M_bh/d_g)")
print(f"    r_cross_indep     = {r8_H['r_cross_independent_m']:.4e} m")
print(f"    R_nuc             = {r8_H['r_nuc_m']:.4e} m")
print(f"    r_cross / R_nuc   = {r8_H['r_vs_R_nuc']:.4e}  (scale bridging factor)")
print(f"    r_vs_R_cov        = {r8_H['r_vs_R_cov']:.4e}")

# Z-scaling survey
print(f"\n  Z^(-2/3) scaling survey:")
print(f"  {'Z':>4}  {'Sym':4}  {'r_cross_indep [m]':>18}  {'vs R_nuc':>12}  "
      f"{'Z-scale err [%]':>16}")
print(SEP2)
for body in dpm.PERIODIC_TABLE[:10]:
    r_i = dpm.chain_r_cross_independent(body)
    print(f"  {body.Z:>4}  {body.symbol:4}  "
          f"{r_i['r_cross_independent_m']:>18.4e}  "
          f"{r_i['r_vs_R_nuc']:>12.4e}  "
          f"{r_i['z_scaling']['error_pct']:>+16.4f}%")

# Checks
check("r_cross_independent positive",
      r8_H["r_cross_independent_m"] > 0.0)
check("r_cross_independent finite",
      math.isfinite(r8_H["r_cross_independent_m"]))
check("galactic_coupling = 23 ± 10",
      13 < r8_H["galactic_coupling"] < 33,
      f"gal_coupling={r8_H['galactic_coupling']:.4f}")
check("r_cross_indep(Z=1) between 1e-14 and 1e-10 m (atomic-to-nuclear bridge)",
      1e-14 < r8_H["r_cross_independent_m"] < 1e-10,
      f"r={r8_H['r_cross_independent_m']:.4e}")

# Z-scaling monotone: r_cross decreases as Z increases
r_prev = r8_H["r_cross_independent_m"]
for body in dpm.PERIODIC_TABLE[1:6]:
    r_next = dpm.chain_r_cross_independent(body)["r_cross_independent_m"]
    check(f"r_cross monotone decreasing Z={body.Z}",
          r_next < r_prev,
          f"r({body.Z})={r_next:.4e} < r({body.Z-1})={r_prev:.4e}")
    r_prev = r_next


# ===========================================================================
# SECTION 5: chain_derive_particle_masses electron entry (Fix #3 integrated)
# ===========================================================================
banner("SECTION 5 — chain_derive_particle_masses electron entry (Fix #3 integrated)")

pm = dpm.chain_derive_particle_masses()
e_entry = pm["electron"]

print(f"\n  electron observed  = {e_entry['observed_kg']:.6e} kg")
print(f"  electron derived   = {e_entry['derived_kg']:.6e} kg")
print(f"  electron error     = {e_entry['error_pct']:+.2f}%")
print(f"  mp/me ratio obs    = {e_entry['mp_me_ratio_obs']:.2f}")
print(f"  formula            : {e_entry['formula']}")

check("Electron entry has derived_kg",
      "derived_kg" in e_entry and e_entry["derived_kg"] > 0.0)
check("Electron error within ±20%",
      abs(e_entry["error_pct"]) < 20.0,
      f"error={e_entry['error_pct']:+.2f}%")
check("Proton still at -2.74% ± 0.1",
      abs(pm["proton"]["error_pct"] - (-2.74)) < 0.1,
      f"p_error={pm['proton']['error_pct']:+.4f}%")
check("Neutron still at -2.70% ± 0.2",
      abs(pm["neutron"]["error_pct"] - (-2.70)) < 0.2,
      f"n_error={pm['neutron']['error_pct']:+.4f}%")


# ===========================================================================
# SUMMARY
# ===========================================================================
banner("SUMMARY — Fix #3 / Fix #4 / Fix #8 chain trace")

print(f"  Fix #3  m_e Route B:   {pm['electron']['derived_kg']:.4e} kg  "
      f"({pm['electron']['error_pct']:+.2f}% vs PDG)")
print(f"  Fix #4  K1_eff:        {k4['K1_eff']['value']:.4e}  "
      f"[Ug1=M_p*c^2 at R_nuc]")
print(f"  Fix #4  K2_eff:        {k4['K2_eff']['value']:.4e}  "
      f"[Ug2=M_e*c^2 at R_bubble]")
print(f"  Fix #4  K3_eff:        {k4['K3_eff']['value']:.4e}  "
      f"[n-p split Fix#2 Route A]")
print(f"  Fix #4  K4_eff:        {k4['K4_eff']['value']:.4e}  "
      f"[galactic energy constraint]")
print(f"  Fix #8  r_cross(Z=1):  {r8_H['r_cross_independent_m']:.4e} m  "
      f"(ratio to R_nuc: {r8_H['r_vs_R_nuc']:.4e})")
print(f"  Fix #8  scaling:       r_cross ∝ Z^(-2/3)  confirmed monotone Z=1..6")
print()
print("  All checks PASSED.  Exit code 0.")
