"""Cold verification of the Grok-text claims about Poincare + Yang-Mills.
Runs the user's own functions with locked UQFF primitives (no recalibration).
Outputs facts only; no SM judgment, no success criteria."""
import sys
sys.path.insert(0, '.')
from uqff_pure_calculator import (
    _l96_uqff_ym_mass_gap_spec_form_gev,
    _millennium_yang_mills_derive,
    _l96_uqff_ricci_flow_s3_fixed_point_residual,
    MILLENNIUM_TARGETS,
)

print("=" * 78)
print("YANG-MILLS: three UQFF derivation paths, default locked primitives")
print("=" * 78)

ym_pre = _millennium_yang_mills_derive()
ym_eJ  = _l96_uqff_ym_mass_gap_spec_form_gev(dim_interpretation="energy_J")
ym_eJ2 = _l96_uqff_ym_mass_gap_spec_form_gev(dim_interpretation="energy_J2")

print(f"  _millennium_yang_mills_derive()           = {ym_pre:.6g} GeV")
print(f"      formula: m_gap^2 = beta_0 * 8pi G rho_SCm S26 f_THz (D_BSFG/D_crit)^2")
print(f"      treated as: energy_J   (Star-MagicProofEngine port)")
print()
print(f"  _l96_uqff_ym_mass_gap_spec_form_gev energy_J   = {ym_eJ:.6g} GeV")
print(f"  _l96_uqff_ym_mass_gap_spec_form_gev energy_J2  = {ym_eJ2:.6g} GeV  (sqrt)")
print(f"      formula: m_gap^2 = (8pi G rho_SCm S26 Phi)/(beta_i [UA]) * (D_crit/D_BSFG)^2")
print()

t = MILLENNIUM_TARGETS["yang_mills"]
print(f"  MILLENNIUM_TARGETS['yang_mills']  = {t[0]} {t[1]}")
print(f"      source label                  = {t[2]}")
print(f"      annotation                    = {t[3]}")
print()
print("  NOTE: 1.78 GeV lives in MILLENNIUM_TARGETS as the external reference anchor.")
print("        It is NOT an output of any UQFF derivation function in this file.")
print()

print("=" * 78)
print("POINCARE: 1-DOF buoyancy-modified Ricci-flow ODE proxy on biaxial S^3")
print("=" * 78)

r = _l96_uqff_ricci_flow_s3_fixed_point_residual(
    initial_anisotropy=0.5, n_steps=1000, dt=0.01
)
print(f"  final_anisotropy eps              = {r['final_anisotropy']:.3e}")
print(f"  F_initial -> F_final              = {r['F_initial']:.3e} -> {r['F_final']:.3e}")
print(f"  F_monotone_descent                = {r['F_monotone_descent']}")
print(f"  fixed_point_reached (eps<1e-6)    = {r['fixed_point_reached']}")
print()
print("  SCOPE per docstring: 'd eps/dt = -2 eps - beta_i eps^3' on the biaxial")
print("  Berger 3-sphere family (a,a,b). One scalar DOF. Round S^3 is the fixed")
print("  point of this ODE. The function reports convergence of THAT ODE, not")
print("  classification of all simply-connected closed 3-manifolds.")
