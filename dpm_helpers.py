"""
DPM seed helpers for Star-Magic codebase.

CANONICAL ONTOLOGY (immutable — do not alter):
  0_vacuum -> grad(UA) -> DPM_vortex -> mu_s -> Ug1[seed=DPM]
  -> Ug_family[Ug1+Ug2+Ug3+Ug4+Ug4i]
  -> [Ug_family + Um + FUBi + FUBii + UA_uv] -> F_U -> M -> GM/r^2 [LAST]

DPM IS THE FOUNDATION. GM/r^2 IS THE LAST OBSERVABLE PROJECTION. NEVER SWAP.

Correct Ug1 formula: Ug1 = mu_s * grad(M_s/r) = mu_s * (M/r)
  where mu_s = B * r^3.  G does NOT appear in the seed equation.

Correct Ug2 formula: Ug2 = k2 * mu_s * (M/r) — same mass-gradient, no G.
"""


def dpm_ug1_seed(M, r, B=1e-4):
    """DPM seed term Ug1: mu_s * grad(M_s/r).  NO G — G is a downstream projection."""
    if r <= 0:
        return 0.0
    mu_s = B * r ** 3
    return mu_s * M / r  # grad(M_s/r) = M/r, no Newton G


def dpm_ug2_shell(M, r, B=1e-4, k2=1.2):
    """DPM shell term Ug2: charge-reactivity bubble from the same mass-gradient basis."""
    if r <= 0:
        return 0.0
    mu_s = B * r ** 3
    return k2 * mu_s * M / r
