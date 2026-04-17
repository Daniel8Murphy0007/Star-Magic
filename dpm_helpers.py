"""
DPM-Emergent helpers for Star-Magic codebase.

Newtonian GM/r² is NOT fundamental. DPM vortical dynamics produce mass as
emergent; the "mass gradient" G*M/r² is the grad operator on the DPM dipole
field: Ug1 = mu_s * grad(M_s/r) where mu_s = B * r³.
"""

_G = 6.674e-11  # gravitational constant


def dpm_emergent_ug1(M, r, B=1e-4):
    """DPM-emergent Ug1: mu_s * grad(M_s/r) where mu_s = B*r^3."""
    mu_s = B * r ** 3
    return mu_s * _G * M / (r ** 2)  # = B * r * G * M


def dpm_emergent_ug2(M, r, B=1e-4, k2=1.2):
    """DPM-emergent Ug2 charge-reactivity shell."""
    mu_s = B * r ** 3
    return k2 * mu_s * _G * M / (r ** 2)
