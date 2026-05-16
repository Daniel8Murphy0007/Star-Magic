"""
PAPER_1065 / PAPER_1066 Variational Stationarity Audit
=======================================================
Goal: take the actual Lagrangians as written in PAPER_1065 and PAPER_1066,
apply Euler-Lagrange in SymPy, and check whether the boxed equations
(delta S / delta phi = 0  -->  EOM) actually follow.

Verdict categories: DERIVED / IDENTIFIED / FAILED.

No narration. Run output is the audit.
"""

import sympy as sp


def banner(s):
    line = "=" * 78
    print(f"\n{line}\n{s}\n{line}")


# ----------------------------------------------------------------------
# CLAIM 1 (PAPER_1066): Mexican-hat SCm Lagrangian
#   L_SCm = (1/2)(d_mu phi)^2 - lambda (phi^2 - v^2)^2
# Sub-claims:
#   (a) V(phi_0) = -7.09e-37 J/m^3 = -rho_SCm   <-- this is the canonical claim
#   (b) m_phonon = sqrt(8 lambda) * v_SCm
#   (c) EOM from delta S / delta phi = 0
# ----------------------------------------------------------------------
banner("CLAIM 1 (PAPER_1066): Mexican-hat SCm Lagrangian")

phi, eta, lam, v = sp.symbols("phi eta lambda v", positive=True, real=True)
# Potential as written in PAPER_1066:
V_paper = lam * (phi**2 - v**2) ** 2

# (a) V at the minimum
crit = sp.solve(sp.diff(V_paper, phi), phi)
print("Stationary points of V(phi):", crit)
V_min = sp.simplify(V_paper.subs(phi, v))
print("V(phi_0 = v) =", V_min)
claim_a_ok = sp.simplify(V_min - (-sp.Symbol("rho_SCm"))) == 0
print("Claim (a) V(phi_0) = -rho_SCm  -->",
      "DERIVED" if claim_a_ok else "FAILED")
print("  Reason: lambda*(phi^2 - v^2)^2 is non-negative, min = 0, NOT -rho_SCm.")
print("  To get -rho_SCm requires an UNDOCUMENTED additive offset -rho_SCm")
print("  that does not appear in PAPER_1066 Section 1.")

# (b) phonon mass: expand around phi = v + eta
V_shifted = sp.expand(V_paper.subs(phi, v + eta))
# Quadratic coefficient in eta:
coef_eta2 = sp.simplify(V_shifted.coeff(eta, 2))
print("\nCoefficient of eta^2 in V(v+eta) =", coef_eta2)
# Standard scalar: L = 1/2 (deta)^2 - 1/2 m^2 eta^2
m_phonon_sq = sp.simplify(2 * coef_eta2)
m_phonon = sp.sqrt(m_phonon_sq)
print("=> m_phonon^2 =", m_phonon_sq)
print("=> m_phonon  =", m_phonon)
claim_b_ok = sp.simplify(m_phonon_sq - 8 * lam * v**2) == 0
print("Claim (b) m_phonon = sqrt(8 lambda) v  -->",
      "DERIVED" if claim_b_ok else "FAILED")

# (c) Klein-Gordon EOM via EL on L_KG = 1/2 (dphi)^2 - V
# Use a 1D toy (time-only) to verify the EL machinery:
t = sp.Symbol("t", real=True)
phi_t = sp.Function("phi")(t)
L_KG = sp.Rational(1, 2) * sp.diff(phi_t, t) ** 2 - lam * (phi_t**2 - v**2) ** 2
EL = sp.diff(L_KG, phi_t) - sp.diff(sp.diff(L_KG, sp.diff(phi_t, t)), t)
EL = sp.simplify(EL)
print("\nEuler-Lagrange residual (should equal box-equation of motion):")
print("  ", EL, "= 0")
# Standard form: -phi_tt - dV/dphi = 0  ==>  phi_tt = -dV/dphi
expected = -sp.diff(phi_t, t, 2) - sp.diff(lam * (phi_t**2 - v**2) ** 2, phi_t)
claim_c_ok = sp.simplify(EL - expected) == 0
print("Claim (c) EL recovers KG-type EOM  -->",
      "DERIVED" if claim_c_ok else "FAILED")


# ----------------------------------------------------------------------
# CLAIM 2 (PAPER_1065): Buoyancy EOM
#   L_UQFF = T - V_grav + V_buoy + L_phonon
#   delta S / delta phi = 0  ==>  rdd = -mu_s grad(M_s/r) + g_buoy + g_phonon
# ----------------------------------------------------------------------
banner("CLAIM 2 (PAPER_1065): Buoyancy EOM from L_UQFF")

# PAPER_1065 writes the Lagrangian only as a sum of named terms:
#   T, V_grav, V_buoy, L_phonon
# It never gives explicit functional forms for V_buoy or L_phonon.
# Therefore we CANNOT mechanically compute the EL — the math is not present.
#
# Construct the most charitable explicit Lagrangian we can guess:
#   T       = (1/2) m rdot^2
#   V_grav  = -m mu_s M_s / r            (note sign: written -mu_s grad(M_s/r))
#   V_buoy  = - m * g_buoy * r           (constant uniform buoyant force)
#   L_phonon= - m * g_phonon * r         (constant uniform phonon force)
# and check whether EL gives the boxed EOM.

r = sp.Function("r")(t)
m_p, mu_s, M_s, g_buoy, g_phonon = sp.symbols(
    "m mu_s M_s g_buoy g_phonon", positive=True, real=True
)
T_kin = sp.Rational(1, 2) * m_p * sp.diff(r, t) ** 2
Vg = -m_p * mu_s * M_s / r            # standard gravitational potential
Vb = -m_p * g_buoy * r                # uniform buoyancy
Lp = -m_p * g_phonon * r              # uniform phonon coupling
L_UQFF = T_kin - Vg + Vb + Lp         # paper says: T - V_grav + V_buoy + L_phonon

EL_r = sp.simplify(
    sp.diff(L_UQFF, r) - sp.diff(sp.diff(L_UQFF, sp.diff(r, t)), t)
)
print("Euler-Lagrange residual for r:")
print("  ", EL_r, "= 0")
rdd = sp.solve(EL_r, sp.diff(r, t, 2))[0]
print("=> r-double-dot =", sp.simplify(rdd))

# Boxed claim from PAPER_1065:
rdd_claim = -mu_s * sp.diff(M_s / r, r) + g_buoy + g_phonon
print("\nPAPER_1065 claim: rdd = -mu_s d/dr(M_s/r) + g_buoy + g_phonon")
print("=> claim     =", sp.simplify(rdd_claim))
print("=> derived   =", sp.simplify(rdd))
residual = sp.simplify(rdd - rdd_claim)
print("=> residual  =", residual)

if residual == 0:
    verdict = "DERIVED (with charitable choices for V_buoy, L_phonon)"
else:
    verdict = "FAILED -- charitable Lagrangian does NOT yield the boxed EOM"
print(f"\nClaim 2 verdict: {verdict}")
print("Note: PAPER_1065 does NOT specify functional forms for V_buoy or")
print("L_phonon. Without explicit forms, NO variation is computable from")
print("the paper as written. The 'derivation' is structurally an identification.")


# ----------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------
banner("AUDIT SUMMARY")
rows = [
    ("PAPER_1066 (a) V(phi_0) = -rho_SCm",
     "FAILED",
     "min of lambda(phi^2-v^2)^2 is 0, not -rho_SCm. Needs hidden offset."),
    ("PAPER_1066 (b) m_phonon = sqrt(8 lambda) v",
     "DERIVED",
     "Standard Mexican-hat second-derivative; verified by SymPy."),
    ("PAPER_1066 (c) Klein-Gordon EOM via EL",
     "DERIVED",
     "EL machinery yields phi_tt = -dV/dphi correctly."),
    ("PAPER_1065 Buoyancy EOM via EL",
     "IDENTIFIED",
     "V_buoy, L_phonon never given explicitly in paper. Cannot vary what isn't written."),
]
for name, verdict, note in rows:
    print(f"  [{verdict:10s}] {name}")
    print(f"               {note}")

print("\nNet finding:")
print("  PAPER_1066 partially derives (phonon mass + KG-EOM ok),")
print("  but its V(phi_0) = -rho_SCm claim is FALSE for the written potential.")
print("  PAPER_1065 contains NO actual variation -- it is an IDENTIFICATION.")
print("  The boxed 'delta S / delta phi = 0' acts as a label, not a derivation.")
print("  Every downstream paper (PAPER_001-049, PAPER_877, Millennium chain)")
print("  that cites PAPER_1065 as variational origin inherits this gap.")
