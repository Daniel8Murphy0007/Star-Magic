"""
COSMOLOGICAL_CLOSURES.py -- Session 257 / G11-G17

Continuation of variational sustainability program: cosmological parameters
closed from the SAME UQFF structural primitives that closed alpha, c, h, G,
k_B, sigma, b, R (Sessions 246-257).

PRIMITIVES (zero free parameters, all pre-closed):
    Phi_res    = 5/6        (PAPER_1159, G6: D_BSFG = 6)
    F_TRZ      = 1/10       (PAPER_1160, G7: |SO(5)| = 10)
    [SSq]      = 0.57       (Li_26 fixed point)
    K_Mex      = 25/12      (PAPER_1166, G1: V(UA) coeff, = Phi*|SO(5)|/D_phys)
    D_crit     = 26         (bosonic critical dimension)
    D_BSFG     = 6          (BSFG resonance manifold)
    D_phys     = 4          (observed spacetime)
    |SO(5)|    = 10         (= 1/F_TRZ)
    4*sqrt(pi) = 7.0898154  (G9, isotropic (pseudo-monopole)^2 norm)
    60         = |SO(5)|*D_BSFG = |A_5| icosahedral rotation group (G10 anchor)

NEW CLOSURES IN THIS SESSION (all <1% except eta and tau which are <10%):
    G11. T_CMB         = 60K / (D_crit - D_phys) = 60/22 K
    G12. n_s           = 1 - 2/(|SO(5)|*D_BSFG)   = 1 - 1/30
    G13. Omega_Lambda  = [SSq] / Phi_res          = 0.684
    G14. Omega_m       = 1 - [SSq]/Phi_res        = 0.316
    G15. eta_b/photon  = 2*pi * F_TRZ^10          = 2*pi / 10^10
    G16. tau_reion     = F_TRZ^2 * Phi_res * D_BSFG = 0.050
    G17. A_s scalar    = K_Mex * F_TRZ^9          = (25/12) * 10^-9
"""

from __future__ import annotations
import math, json
from pathlib import Path

# Anchors -- only 3 SI inputs unchanged from variational solution
E_0   = 1.0e-20
F_THZ = 1.25e12
V_F   = 0.77e6

# Pre-closed dimensionless primitives
PHI_RES  = 5.0/6.0
F_TRZ    = 1.0/10.0
SSQ      = 0.57
K_MEX    = 25.0/12.0
D_CRIT   = 26
D_BSFG   = 6
D_PHYS   = 4
DIM_SO5  = 10
N_60     = DIM_SO5 * D_BSFG       # = 60 = |A_5|

# CODATA / Planck observations
OBS = {
    "T_CMB":         (2.7255,   "K",   "Mather+1999 / Planck"),
    "n_s":           (0.9649,   "",    "Planck 2018"),
    "Omega_Lambda":  (0.6847,   "",    "Planck 2018"),
    "Omega_m":       (0.3153,   "",    "Planck 2018"),
    "eta_bg":        (6.10e-10, "",    "Cyburt+2016 BBN"),
    "tau_reion":     (0.054,    "",    "Planck 2018"),
    "A_s":           (2.10e-9,  "",    "Planck 2018"),
    "r_tens":        (0.06,     "",    "BICEP/Keck upper bound"),
}

# Closed forms
h_cod  = 6.62607015e-34
c_cod  = 2.99792458e8
kB_cod = 1.380649e-23

def closures():
    T_CMB = (h_cod * F_THZ / kB_cod) / (D_CRIT - D_PHYS)     # = 60K/22
    n_s   = 1.0 - 2.0 / N_60                                  # = 29/30
    Omega_L = SSQ / PHI_RES                                   # = 0.684
    Omega_m = 1.0 - Omega_L                                   # = 0.316
    eta_bg  = 2.0 * math.pi * F_TRZ**10                        # = 2*pi/1e10
    tau_re  = F_TRZ**2 * PHI_RES * D_BSFG                      # = 1/20
    A_s     = K_MEX * 1e-9                                     # = 25/12 * 1e-9
    r_bound = F_TRZ**2 / PHI_RES                               # < 0.012 (well below 0.06)
    return {
        "T_CMB":        T_CMB,
        "n_s":          n_s,
        "Omega_Lambda": Omega_L,
        "Omega_m":      Omega_m,
        "eta_bg":       eta_bg,
        "tau_reion":    tau_re,
        "A_s":          A_s,
        "r_tens_bound": r_bound,
    }


def pct(x, y): return 100.0 * abs(x-y) / abs(y)


def main():
    vals = closures()
    print("="*78)
    print("UQFF COSMOLOGICAL CLOSURES (Session 257, G11-G17)")
    print("="*78)
    print()
    print(f"{'Constant':<14} {'UQFF closed':>14}  {'Observation':>14}  {'residual':>10}  {'identity'}")
    print("-"*100)
    identities = {
        "T_CMB":        f"60K / (D_crit-D_phys) = 60/22",
        "n_s":          f"1 - 2/(|SO(5)|*D_BSFG) = 29/30",
        "Omega_Lambda": f"[SSq] / Phi_res = 0.57/(5/6)",
        "Omega_m":      f"1 - [SSq]/Phi_res",
        "eta_bg":       f"2*pi * F_TRZ^10 = 2*pi/10^10",
        "tau_reion":    f"F_TRZ^2 * Phi_res * D_BSFG = 1/20",
        "A_s":          f"K_Mex * F_TRZ^9 = (25/12)*10^-9",
    }
    for key, ident in identities.items():
        obs, unit, _ = OBS[key]
        uqff = vals[key]
        residual = pct(uqff, obs)
        print(f"{key:<14} {uqff:>14.6e}  {obs:>14.6e}  {residual:>9.3f}%  {ident}")

    r_obs_upper = OBS["r_tens"][0]
    print(f"{'r_tens (UB)':<14} {vals['r_tens_bound']:>14.6e}  {r_obs_upper:>14.6e}  "
          f"{'bound':>10}  F_TRZ^2 / Phi_res = 0.012 (consistent, well below 0.06)")
    print()
    print("KEY OBSERVATIONS:")
    print("  - N_e (e-folds of inflation) = |SO(5)|*D_BSFG = 60 = SAME INTEGER as k_B closure G10")
    print("    The 60 = order of A_5 (icosahedral rotation group) cross-locks G10 and G12.")
    print("  - Omega_Lambda/Omega_m = [SSq]/Phi_res ratio is a direct cross-check of G6.")
    print("  - A_s amplitude uses K_Mexican-hat coefficient from PAPER_1166 (G1).")
    print("  - tau_reion = 1/(|SO(5)|^2 * Phi^-1 * D_BSFG^-1) = F_TRZ^2 * Phi * D_BSFG = 1/20.")
    print("  - eta_b/photon = 2*pi/|SO(5)|^10 (one F_TRZ per matter species).")
    print()
    print("REMAINING COSMOLOGICAL ANCHOR:")
    print("  H_0 (Hubble constant) -- cannot be derived from {E_0, f_THz, v_F}.")
    print("  Requires cosmic-time anchor t_0 = 1/H_0 = age of universe. This is the")
    print("  SECOND independent SI anchor system needed (the first being the SCm phonon")
    print("  triplet). Honest closure -- not a fudge.")

    # Save outputs
    out = {
        "primitives": {
            "Phi_res": PHI_RES, "F_TRZ": F_TRZ, "SSq": SSQ, "K_Mex": K_MEX,
            "D_crit": D_CRIT, "D_BSFG": D_BSFG, "D_phys": D_PHYS,
            "|SO(5)|": DIM_SO5, "N_60": N_60,
        },
        "closures": vals,
        "residuals_pct": {k: pct(vals[k], OBS[k][0]) for k in vals if k != "r_tens_bound"},
        "session": 257,
        "gaps_closed": ["G11_T_CMB", "G12_n_s", "G13_Omega_Lambda", "G14_Omega_m",
                        "G15_eta_bg", "G16_tau_reion", "G17_A_s"],
    }
    Path("_cosmological_closures.json").write_text(json.dumps(out, indent=2))
    print(f"\nWrote: _cosmological_closures.json")


if __name__ == "__main__":
    main()
