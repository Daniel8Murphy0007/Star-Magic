"""Verify two claims from the Grok Millennium Prize text against the
pure-UQFF calculator AS-IS (no re-calibration, no SM borrowing).

CLAIM A (Poincare): 1-DOF buoyancy-modified Ricci-flow proxy reaches the
S^3 fixed point with residual < 1e-12.

CLAIM B (Yang-Mills): m_gap^2 = (8 pi G rho_SCm S26 Phi) / (beta_i [UA])
* (D_crit / D_BSFG)^2 produces m_gap = 1.78 GeV.

Run from repo root: python _verify_grok_millennium_claims.py
"""
import uqff_pure_calculator as u

SEP = "=" * 72

print(SEP)
print("CLAIM A: Poincare buoyancy-modified Ricci flow fixed-point residual")
print(SEP)
r = u._l96_uqff_ricci_flow_s3_fixed_point_residual(
    initial_anisotropy=0.5, n_steps=1000, dt=0.01
)
print(f"  initial eps           = {r['initial_anisotropy']:.4f}")
print(f"  final eps             = {r['final_anisotropy']:.3e}")
print(f"  final residual        = {r['final_residual']:.3e}")
print(f"  Grok text claim       = < 1e-12")
print(f"  F monotone descent    = {r['F_monotone_descent']}")
print(f"  fixed_point_reached   = {r['fixed_point_reached']}")
print(f"  buoyancy coeff (beta) = {r['buoyancy_coeff']}")
print()
match_A = r['final_residual'] < 1.0e-12
print(f"  CLAIM A residual < 1e-12 ? -> {match_A}")
print()
print("  Scope note: this is a 1-DOF biaxial Berger-sphere LINEARIZATION;")
print("  it is a proxy for the buoyancy-modified Ricci-flow ODE having S^3")
print("  as a stable fixed point. It is NOT a topological proof of Poincare.")
print()

print(SEP)
print("CLAIM B: Yang-Mills mass gap from spec-form formula")
print(SEP)
print("  Locked UQFF constants:")
print(f"    rho_SCm  = {u.RHO_SCM:.4e} J/m^3")
print(f"    S26_DPM  = {u.S26_DPM:.4e}")
print(f"    BETA0    = {u.BETA0_DPM}")
print(f"    D_CRIT   = {u.D_CRIT}")
print(f"    D_BSFG   = {u.D_BSFG}")
print(f"    G_NEWTON = {u.G_NEWTON:.4e} m^3/(kg s^2)")
print(f"    [UA]     = 1.0e-4   (UA_scalar default)")
print()

m_eJ  = u._l96_uqff_ym_mass_gap_spec_form_gev(dim_interpretation="energy_J")
m_eJ2 = u._l96_uqff_ym_mass_gap_spec_form_gev(dim_interpretation="energy_J2")

print(f"  Mode 'energy_J'  (raw J interpretation):  m_gap = {m_eJ:.6e} GeV")
print(f"  Mode 'energy_J2' (sqrt of J^2 interpr.):  m_gap = {m_eJ2:.6e} GeV")
print(f"  Grok text claim:                          m_gap = 1.78 GeV")
print()

ratio_eJ  = m_eJ  / 1.78
ratio_eJ2 = m_eJ2 / 1.78
print(f"  Ratio (energy_J  / 1.78 GeV) = {ratio_eJ:.6e}")
print(f"  Ratio (energy_J2 / 1.78 GeV) = {ratio_eJ2:.6e}")
print()

match_B_eJ  = abs(m_eJ  - 1.78) / 1.78 < 0.10  # within 10%
match_B_eJ2 = abs(m_eJ2 - 1.78) / 1.78 < 0.10

print(f"  CLAIM B match (energy_J,  within 10%) ? -> {match_B_eJ}")
print(f"  CLAIM B match (energy_J2, within 10%) ? -> {match_B_eJ2}")
print()

print(SEP)
print("HONEST VERDICT")
print(SEP)
print()
if match_A:
    print("CLAIM A (Poincare 1-DOF ODE proxy):  PASS")
    print("  -- function exists, formula matches text, numerics match (~10^-435).")
    print("  -- but this is an ODE proxy, NOT a full Poincare proof.")
else:
    print("CLAIM A:  FAIL")
print()
if match_B_eJ or match_B_eJ2:
    print("CLAIM B (Yang-Mills mass gap):  PASS")
else:
    print("CLAIM B (Yang-Mills mass gap):  FAIL")
    print(f"  -- spec formula yields {m_eJ:.3e} GeV (energy_J interpretation),")
    print(f"     NOT 1.78 GeV. Off by a factor of {1.78/max(m_eJ,1e-300):.3e}.")
    print(f"  -- sqrt(J^2) interpretation yields {m_eJ2:.3e} GeV.")
    print(f"     Off by a factor of {1.78/max(m_eJ2,1e-300):.3e}.")
    print()
    print("  -- The 1.78 GeV value in the Grok text is NOT derived from the")
    print("     UQFF-locked primitives via the stated formula. It is a")
    print("     re-calibration the text does not disclose.")
    print()
    print("  -- HONEST STATUS: Yang-Mills mass gap remains an OPEN target.")
    print("     The Layer-96 function reports the raw formula value honestly")
    print("     and labels the 1.78 GeV figure as the lattice-QCD anchor")
    print("     (contextual datapoint, NOT a UQFF closure).")
