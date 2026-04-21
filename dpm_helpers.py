"""
DPM-emergent helpers for Star-Magic codebase.

CANONICAL ONTOLOGY LOCK (v1) alignment:
1) Starting state is zero-mass vacuum; no moving-mass Newtonian seed.
2) Mass emergence precedes motion.
3) Promotion order is fixed: Ug1 -> Ug2 + Ug3 + Ug4 (+ Ug4_i).
4) Gravity family is assembled simultaneously before downstream projections.
5) GM/r^2 is observationally useful but only as a reduced projection term,
   never a foundational seed equation.

This helper computes DPM-emergent terms. For Ug1, the quantity G*M/r^2 is
used as the mass-gradient operator within the dipole formulation:
Ug1 = mu_s * grad(M_s/r), where mu_s = B*r^3.
"""

_G = 6.674e-11  # gravitational constant


def dpm_emergent_ug1(M, r, B=1e-4):
    """Compute Ug1 in DPM form; GM/r^2 is treated as gradient, not seed force."""
    mu_s = B * r ** 3
    return mu_s * _G * M / (r ** 2)  # = B * r * G * M


def dpm_emergent_ug2(M, r, B=1e-4, k2=1.2):
    """Compute Ug2 charge-reactivity shell from the same emergent basis."""
    mu_s = B * r ** 3
    return k2 * mu_s * _G * M / (r ** 2)
