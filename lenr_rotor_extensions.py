#!/usr/bin/env python3
"""
lenr_rotor_extensions.py — LENR Phonon-Rotor Coupling Extensions

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
Standalone LENR rotor module with phonon-enhanced cross sections and COP
(Coefficient of Performance) calculations.

Gap closed:
  - σ_rotor = σ₀·exp(-(ω-ω_SCm)²/(2Γ²))·S₂₆⁽³⁾  (phonon-enhanced cross section)
  - COP phonon coupling: COP = P_out/P_in with SCm resonance boost
  - Widom-Larsen UQFF integration

Physics:
  The LENR rotor cross section is enhanced by the SCm 1.25 THz phonon
  resonance. At the resonance frequency, the Gaussian envelope is unity
  and the cross section is amplified by S₂₆⁽³⁾. The COP includes the
  buoyancy-sector energy contribution (F_{U,Bi}/F_U).

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
eV        = 1.602e-19
M_P       = 1.673e-27       # proton mass
M_E       = 9.109e-31       # electron mass
M_N       = 1.675e-27       # neutron mass

OMEGA_SCM = 2 * PI * 1.25e12
GAMMA_0   = 2 * PI * 0.1e12
SSQ       = 0.57
F_UBI_RATIO = 0.6

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


def ramanujan_Rn(n: int, k: int = 3) -> float:
    """Ramanujan acceleration factor."""
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


# ── §1  PHONON-ENHANCED CROSS SECTION ─────────────────────────────────────

class LENRPhononCrossSection:
    """Phonon-enhanced LENR cross section.

    σ_rotor(ω) = σ₀ · exp(-(ω - ω_SCm)² / (2Γ²)) · S₂₆⁽³⁾

    At resonance (ω = ω_SCm): σ_rotor = σ₀ · S₂₆⁽³⁾
    The SCm phonon mode enhances the tunnelling probability by providing
    coherent lattice energy to overcome the Coulomb barrier.
    """

    def __init__(self, sigma_0: float = 1e-28, omega_scm: float = OMEGA_SCM,
                 gamma: float = GAMMA_0):
        self.sigma_0 = sigma_0
        self.omega_scm = omega_scm
        self.gamma = gamma

    def cross_section(self, omega: float) -> float:
        """Compute σ_rotor at frequency omega."""
        dw = omega - self.omega_scm
        gaussian = math.exp(-dw ** 2 / (2 * self.gamma ** 2))
        return self.sigma_0 * gaussian * S26_3RD

    def cross_section_density(self, omega: float, rho: float) -> float:
        """Compute σ_n(ω, ρ) — cross section with explicit frequency AND density.

        σ_n(ω, ρ) = σ₀ · exp(-(ω-ω_SCm)²/(2Γ²)) · S₂₆⁽³⁾ · (1 + [SSq]·ρ/(ρ+ρ_ref))

        The density factor provides a monotonic enhancement: at high ρ the
        cross section approaches σ_rotor·(1+[SSq]), while at ρ→0 it
        reduces to the frequency-only form.

        Args:
            omega: angular frequency (rad/s)
            rho: local lattice/material density (kg/m³)

        Returns:
            σ_n in m²
        """
        rho_ref = 1e6  # reference density scale (kg/m³)
        base = self.cross_section(omega)
        density_factor = 1.0 + SSQ * rho / (rho + rho_ref)
        return base * density_factor

    def enhancement_factor(self, omega: float) -> float:
        """Enhancement over bare σ₀."""
        return self.cross_section(omega) / self.sigma_0 if self.sigma_0 > 0 else 0.0

    def frequency_sweep(self, f_min_THz: float = 0.5, f_max_THz: float = 2.0,
                        n_points: int = 50) -> list:
        """Sweep frequency and compute cross section."""
        results = []
        for i in range(n_points + 1):
            f = f_min_THz + (f_max_THz - f_min_THz) * i / n_points
            omega = 2 * PI * f * 1e12
            sigma = self.cross_section(omega)
            results.append({
                "f_THz": f,
                "omega_rad_s": omega,
                "sigma_m2": sigma,
                "enhancement": sigma / self.sigma_0 if self.sigma_0 > 0 else 0,
            })
        return results

    def compute(self, dataset: dict) -> dict:
        """Full cross section calculation."""
        omega = float(dataset.get("omega", self.omega_scm))
        rho = float(dataset.get("rho_kg_m3", 0.0))
        sigma = self.cross_section(omega)
        sigma_rho = self.cross_section_density(omega, rho)
        sigma_res = self.cross_section(self.omega_scm)
        sweep = self.frequency_sweep()
        return {
            "sigma_rotor_m2": sigma,
            "sigma_n_density_m2": sigma_rho,
            "sigma_resonance_m2": sigma_res,
            "sigma_0_m2": self.sigma_0,
            "rho_kg_m3": rho,
            "density_enhancement": sigma_rho / sigma if sigma > 0 else 1.0,
            "enhancement_at_resonance": sigma_res / self.sigma_0,
            "sweep_len": len(sweep),
            "primary_equations": [
                f"σ_rotor(ω) = σ₀·exp(-(ω-ω_SCm)²/(2Γ²))·S₂₆⁽³⁾ = {sigma:.6e} m²",
                f"σ_n(ω,ρ) = σ_rotor·(1+[SSq]·ρ/(ρ+ρ_ref)) = {sigma_rho:.6e} m²",
                f"σ₀ = {self.sigma_0:.6e} m²",
                f"S₂₆⁽³⁾ = {S26_3RD:.6e}",
                f"ρ = {rho:.2e} kg/m³, density factor = {sigma_rho/sigma if sigma>0 else 1:.6f}",
            ],
        }


# ── §2  WIDOM-LARSEN UQFF ─────────────────────────────────────────────────

class WidomLarsenUQFF:
    """Widom-Larsen neutron production with UQFF phonon enhancement.

    The Widom-Larsen mechanism: e⁻ + p⁺ → n + ν_e
    requires electron mass renormalisation m*_e > m_n - m_p ≈ 1.293 MeV/c².

    SCm phonon field provides coherent energy to electrons on the surface:
    m*_e = m_e + ℏω_SCm·S₂₆⁽³⁾/(c²) · N_coherent

    where N_coherent is the number of phonon modes coupling coherently.
    """

    MASS_DIFF = (M_N - M_P) * C ** 2 / eV  # ≈ 1.293 MeV

    def __init__(self):
        self.delta_scm = HBAR * OMEGA_SCM

    def effective_electron_mass(self, N_coherent: int = 1) -> float:
        """Compute m*_e in kg."""
        return M_E + (self.delta_scm * S26_3RD * N_coherent) / C ** 2

    def threshold_N_coherent(self) -> int:
        """Minimum N_coherent for neutron production."""
        dm_required = (M_N - M_P) - M_E  # mass deficit in kg
        if dm_required <= 0:
            return 0
        per_phonon = (self.delta_scm * S26_3RD) / C ** 2
        if per_phonon <= 0:
            return -1
        return math.ceil(dm_required / per_phonon)

    def compute(self, dataset: dict) -> dict:
        """Full Widom-Larsen UQFF calculation."""
        N_coh = int(dataset.get("N_coherent", 1000))
        m_star = self.effective_electron_mass(N_coh)
        m_star_eV = m_star * C ** 2 / eV
        N_thresh = self.threshold_N_coherent()
        above_threshold = m_star > (M_N - M_P)

        return {
            "m_star_e_kg": m_star,
            "m_star_e_MeV": m_star_eV / 1e6,
            "m_e_MeV": M_E * C ** 2 / (eV * 1e6),
            "mass_diff_MeV": self.MASS_DIFF / 1e6,
            "N_coherent": N_coh,
            "N_threshold": N_thresh,
            "above_threshold": above_threshold,
            "primary_equations": [
                f"m*_e = m_e + ℏω_SCm·S₂₆⁽³⁾·N_coh/c² = {m_star_eV / 1e6:.6f} MeV/c²",
                f"Threshold (m_n - m_p) = {self.MASS_DIFF / 1e6:.6f} MeV/c²",
                f"N_coherent = {N_coh}, N_threshold = {N_thresh}",
                f"Above threshold: {above_threshold}",
            ],
        }


# ── §3  COP CALCULATOR ────────────────────────────────────────────────────

class LENRCOPCalculator:
    """Coefficient of Performance with SCm phonon coupling.

    COP = P_out / P_in
    P_out = P_nuclear + P_buoyancy
    P_buoyancy = (F_{UBi}/F_U) · S₂₆⁽³⁾ · P_phonon

    The UQFF buoyancy sector provides an additional energy pathway
    beyond nuclear reactions alone.
    """

    def compute(self, dataset: dict) -> dict:
        """COP calculation.

        Dataset keys:
            P_in_W      : input power (W)
            P_nuclear_W : nuclear output power (W)
            P_phonon_W  : phonon mode power (W, optional)
        """
        P_in = float(dataset.get("P_in_W", 100))
        P_nuc = float(dataset.get("P_nuclear_W", 50))
        P_phonon = float(dataset.get("P_phonon_W", 10))
        fubi = float(dataset.get("FUBi_ratio", F_UBI_RATIO))

        P_buoyancy = fubi * S26_3RD * P_phonon
        P_out = P_nuc + P_buoyancy
        COP = P_out / P_in if P_in > 0 else 0.0

        return {
            "COP": COP,
            "P_out_W": P_out,
            "P_in_W": P_in,
            "P_nuclear_W": P_nuc,
            "P_buoyancy_W": P_buoyancy,
            "primary_equations": [
                f"COP = P_out/P_in = {COP:.6f}",
                f"P_out = P_nuclear + P_buoyancy = {P_out:.6e} W",
                f"P_buoyancy = (F_{{UBi}}/F_U)·S₂₆⁽³⁾·P_phonon = {P_buoyancy:.6e} W",
            ],
        }


# ── §4  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("LENR Phonon-Rotor Coupling Extensions")
    print("=" * 72)

    xs = LENRPhononCrossSection()
    result = xs.compute({"omega": OMEGA_SCM})
    print("\nPhonon-enhanced cross section:")
    for eq in result["primary_equations"]:
        print(f"  {eq}")

    wl = WidomLarsenUQFF()
    result_wl = wl.compute({"N_coherent": 1000})
    print("\nWidom-Larsen UQFF:")
    for eq in result_wl["primary_equations"]:
        print(f"  {eq}")

    cop = LENRCOPCalculator()
    result_cop = cop.compute({"P_in_W": 100, "P_nuclear_W": 50, "P_phonon_W": 10})
    print("\nCOP calculation:")
    for eq in result_cop["primary_equations"]:
        print(f"  {eq}")

    print("\n✓ LENR rotor extensions complete")
