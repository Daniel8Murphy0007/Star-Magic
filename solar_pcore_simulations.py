#!/usr/bin/env python3
"""
solar_pcore_simulations.py — Solar / Planetary Core Composite Density Calculator

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
Computes the Pcore composite density from the SCm vacuum structure.

Gap closed:
  - ρ_Pcore = ρ_SCm · S₂₆^(3)([SSq]) · Φ_{1.25THz}(ω,T) · (F_{U,Bi}/F_U)
  - Per-planet calibration via P_core parameter
  - Temperature-dependent phonon occupation

Physics:
  The core density of a planetary or stellar body ρ_Pcore is set by the
  SCm vacuum condensate density ρ_SCm, modulated by the 26D spectral sum
  S₂₆^(3), the 1.25 THz phonon occupation Φ, and the buoyancy ratio.
  The P_core calibration factor normalises to the observed central density.

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
eV        = 1.602e-19
SOLAR_MASS = 1.989e30

OMEGA_SCM = 2 * PI * 1.25e12
GAMMA_0   = 2 * PI * 0.1e12
SSQ       = 0.57
H_SCM     = 0.99
F_UBI_RATIO = 0.6

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor R_n^{(26,k)}."""
    prefactor = (2 * PI) ** (n / 6.0) / math.factorial(min(n, 170))
    correction = 0.0
    for m in range(1, k + 1):
        inner = 0.0
        for j in range(1, 27):
            sign = (-1) ** (j + 1)
            binom = math.comb(26, j)
            inner += sign * binom * math.factorial(26 - j) / n ** j
        correction += inner / n ** (26 * m)
    return prefactor * (1.0 + correction)


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 28))


# ── §1  PHONON OCCUPATION ─────────────────────────────────────────────────

class PhononOccupation:
    """Phonon spectral weight Φ_{1.25THz}(ω, T).

    Φ(ω, T) = exp(-((ω - ω_SCm)² / (2Γ²))) · n_BE(ω, T)

    where n_BE = 1/(exp(ℏω/k_BT) - 1) is the Bose-Einstein occupation.
    At the SCm resonance (ω = ω_SCm), the Gaussian envelope is unity.
    """

    def __init__(self, omega_scm: float = OMEGA_SCM, gamma: float = GAMMA_0):
        self.omega_scm = omega_scm
        self.gamma = gamma

    def gaussian_envelope(self, omega: float) -> float:
        """Gaussian line shape centred on ω_SCm."""
        dw = omega - self.omega_scm
        return math.exp(-dw ** 2 / (2 * self.gamma ** 2))

    def bose_einstein(self, omega: float, T: float) -> float:
        """Bose-Einstein occupation number n_BE(ω, T)."""
        if T <= 0:
            return 0.0
        x = HBAR * abs(omega) / (K_B * T)
        if x > 500:
            return 0.0
        return 1.0 / (math.exp(x) - 1) if x > 1e-12 else K_B * T / (HBAR * abs(omega))

    def phi(self, omega: float, T: float) -> float:
        """Full phonon occupation Φ(ω, T)."""
        return self.gaussian_envelope(omega) * self.bose_einstein(omega, T)

    def phi_resonance(self, T: float) -> float:
        """Phonon occupation at resonance (ω = ω_SCm)."""
        return self.bose_einstein(self.omega_scm, T)


# ── §2  PCORE COMPOSITE DENSITY ───────────────────────────────────────────

class PcoreCompositeDensity:
    """Compute core density from SCm vacuum condensate.

    ρ_Pcore = ρ_SCm · S₂₆^(3) · Φ(ω_SCm, T_core) · (F_{U,Bi}/F_U) · P_core

    Parameters:
    ───────────
    ρ_SCm :     SCm vacuum condensate density (kg/m³)
    S₂₆^(3) :  Ramanujan-accelerated 26D spectral sum
    Φ :         1.25 THz phonon occupation at core temperature
    F_{UBi}/F_U : buoyancy ratio (canonical 0.6)
    P_core :    per-body calibration factor (Sun=1.0, others << 1)
    """

    def __init__(self, fubi_ratio: float = F_UBI_RATIO):
        self.fubi_ratio = fubi_ratio
        self.phonon = PhononOccupation()

    def rho_scm(self, M: float, R: float) -> float:
        """Estimate ρ_SCm from body mass M (kg) and radius R (m).

        ρ_SCm = M / (4π/3 · R³) is the mean density. The SCm vacuum
        condensate is calibrated to the mean bulk density.
        """
        V = (4.0 / 3.0) * PI * R ** 3
        return M / V if V > 0 else 0.0

    def compute_rho_pcore(self, rho_scm: float, T_core: float,
                          P_core: float = 1.0,
                          fubi_ratio: float = None) -> float:
        """Compute ρ_Pcore."""
        r = fubi_ratio if fubi_ratio is not None else self.fubi_ratio
        phi = self.phonon.phi_resonance(T_core)
        return rho_scm * S26_3RD * phi * r * P_core

    def density_profile(self, rho_scm: float, T_core: float, T_surface: float,
                        R: float, P_core: float = 1.0,
                        n_shells: int = 20) -> list:
        """Radial density profile from core to surface.

        Temperature model: T(r) = T_core · (1 - (r/R)²) + T_surface · (r/R)²
        """
        results = []
        for i in range(n_shells + 1):
            frac = i / n_shells
            r = frac * R
            T = T_core * (1 - frac ** 2) + T_surface * frac ** 2
            phi = self.phonon.phi_resonance(T) if T > 0 else 0.0
            rho = rho_scm * S26_3RD * phi * self.fubi_ratio * P_core
            results.append({
                "r_frac": frac,
                "r_m": r,
                "T_K": T,
                "Phi": phi,
                "rho_Pcore": rho,
            })
        return results

    def compute(self, dataset: dict) -> dict:
        """Full Pcore calculation from dataset.

        Expected dataset keys:
            M       : mass (kg)
            R       : radius (m)
            T_core  : core temperature (K)
            T_surface : surface temperature (K, optional, default 300)
            P_core  : calibration factor (optional, default 1.0)
            FUBi_ratio : buoyancy ratio (optional, default 0.6)
        """
        M = float(dataset.get("M", SOLAR_MASS))
        R = float(dataset.get("R", 6.957e8))
        T_core = float(dataset.get("T_core", 1.57e7))
        T_surface = float(dataset.get("T_surface", 5778))
        P_core = float(dataset.get("P_core", 1.0))
        fubi = float(dataset.get("FUBi_ratio", self.fubi_ratio))
        self.fubi_ratio = fubi

        rho_mean = self.rho_scm(M, R)
        rho_pcore = self.compute_rho_pcore(rho_mean, T_core, P_core, fubi)
        phi_core = self.phonon.phi_resonance(T_core)
        profile = self.density_profile(rho_mean, T_core, T_surface, R, P_core)

        return {
            "rho_mean_kg_m3": rho_mean,
            "rho_Pcore_kg_m3": rho_pcore,
            "T_core_K": T_core,
            "Phi_core": phi_core,
            "S26_3rd": S26_3RD,
            "FUBi_ratio": fubi,
            "P_core": P_core,
            "radial_profile_len": len(profile),
            "primary_equations": [
                "ρ_Pcore = ρ_SCm · S₂₆⁽³⁾ · Φ(ω_SCm, T_core) · (F_{UBi}/F_U) · P_core",
                f"ρ_SCm (mean density) = {rho_mean:.4e} kg/m³",
                f"S₂₆⁽³⁾ = {S26_3RD:.6e}",
                f"Φ(ω_SCm, T_core={T_core:.2e} K) = {phi_core:.6e}",
                f"P_core = {P_core}",
                f"ρ_Pcore = {rho_pcore:.6e} kg/m³",
            ],
        }


# ── §3  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("Solar / Planetary Pcore Composite Density Calculator")
    print("=" * 72)

    bodies = [
        {"name": "Sun",     "M": 1.989e30, "R": 6.957e8, "T_core": 1.57e7,  "T_surface": 5778,  "P_core": 1.0},
        {"name": "Earth",   "M": 5.972e24, "R": 6.371e6, "T_core": 6000,    "T_surface": 288,   "P_core": 1e-3},
        {"name": "Jupiter", "M": 1.898e27, "R": 6.991e7, "T_core": 3.6e4,   "T_surface": 165,   "P_core": 1e-3},
        {"name": "Neptune", "M": 1.024e26, "R": 2.462e7, "T_core": 7100,    "T_surface": 72,    "P_core": 1e-3},
    ]

    calc = PcoreCompositeDensity()
    for body in bodies:
        name = body.pop("name")
        result = calc.compute(body)
        print(f"\n{name}:")
        for eq in result["primary_equations"]:
            print(f"  {eq}")

    print("\n✓ Pcore composite density calculation complete")
