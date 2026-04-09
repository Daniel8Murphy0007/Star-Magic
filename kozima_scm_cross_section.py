"""
Kozima SCm-Modulated Neutron-Drop Cross-Section Calculator

Session 204 | Daniel Murphy
PURPOSE: Extends the existing NeutronCaptureRateCalculator (CP2 L29388) with
         proper SCm manifold modulation via the 26-level VDS enhancement term.

EXISTING (CP2):
  σ_n(ω) = σ₀ · (ω/ω_LENR)² · exp[-(ω - ω_LENR)² / (2·Δω²)]
  σ_n(ρ) = σ₀ · (ρ/ρ₀)             ← linear, no VDS

THIS MODULE ADDS:
  σ_n^SCm(ω,n) = σ₀ · exp[-(ω - ω_SCm)² / (2Γ²)] · (1 + [SSq]·n/26)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^
                  Pure Gaussian resonance at 1.25 THz   VDS 26-level factor

  σ_n^SCm(ρ) = σ₀ · (ρ_SCm/ρ_UA)^[SSq] · exp(-[SSq]·ρ/ρ_crit)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^
               Power-law VDS factor          Exponential cutoff at ρ_crit

  F_neutron^SCm = N_n · σ_n^SCm(ω) · Φ_phonon · (F_{U,Bi}/F_U - 1)
                  buoyancy-reversal coupled neutron production force

COMPARISON:  Both base (CP2) and SCm-modulated forms computed for cross-validation.

ARCHITECTURE: CondensedPhysics.py rules — pure calculator, no hardcoded systems.
"""

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
K_B = 1.381e-23

# UQFF calibrated constants
SSQ = 0.57
KAPPA = 5.787e-9
H_SCM = 0.99
BETA_I = 0.603
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37      # kg/m³ — SCm superconductive vacuum
RHO_UA = 7.09e-36       # kg/m³ — UA aether vacuum
RHO_VAC_SCM = 9.47e-27  # kg/m³ — VDS SCm density
RHO_VAC_UA = 5e-27      # kg/m³ — VDS UA density

# LENR parameters
F_LENR_THZ = 1.25e12    # Hz — SCm phonon resonance
OMEGA_SCM = 2 * PI * F_LENR_THZ  # rad/s — angular frequency
GAMMA_DEFAULT = 0.1e12 * 2 * PI  # rad/s — resonance width (~0.1 THz)
SIGMA_0 = 1e-4          # base cross-section (dimensionless scaled)
B_CRIT = 4.4e13         # T — magnetar critical field

# Default parameters
N_LEVELS = 26           # VDS dimensionality
RHO_CRIT_DEFAULT = 1e-20  # kg/m³ — SCm activation density threshold


# ── §1  SCm-MODULATED FREQUENCY CROSS-SECTION ────────────────────────────

class KozimaSCmCrossSection:
    """
    SCm-modulated Kozima neutron-drop cross-section calculator.

    Provides three cross-section models:
      1. Base (CP2):   σ_n(ω) = σ₀ · (ω/ω_LENR)² · exp[-(ω-ω_LENR)²/(2Δω²)]
      2. SCm Gaussian: σ_n^SCm(ω,n) = σ₀ · exp[-(ω-ω_SCm)²/(2Γ²)] · (1+[SSq]·n/26)
      3. Density:      σ_n^SCm(ρ) = σ₀ · (ρ_SCm/ρ_UA)^[SSq] · exp(-[SSq]·ρ/ρ_crit)

    All three are computed for cross-validation.
    """

    def __init__(self,
                 sigma_0: float = SIGMA_0,
                 omega_scm: float = OMEGA_SCM,
                 gamma: float = GAMMA_DEFAULT,
                 ssq: float = SSQ,
                 n_levels: int = N_LEVELS,
                 rho_scm: float = RHO_VAC_SCM,
                 rho_ua: float = RHO_VAC_UA,
                 rho_crit: float = RHO_CRIT_DEFAULT):
        self.sigma_0 = sigma_0
        self.omega_scm = omega_scm
        self.gamma = gamma
        self.ssq = ssq
        self.n_levels = n_levels
        self.rho_scm = rho_scm
        self.rho_ua = rho_ua
        self.rho_crit = rho_crit

    # ── 1a. Base CP2 cross-section (for comparison) ──────────────────────

    def compute_base_cross_section(self, omega: float,
                                   delta_omega: float = 0.1e12 * 2 * PI
                                   ) -> Dict[str, Any]:
        """
        Base NeutronCaptureRateCalculator form (CP2 L29406):
          σ_n(ω) = σ₀ · (ω/ω_LENR)² · exp[-(ω - ω_LENR)²/(2·Δω²)]
        """
        freq_ratio = (omega / self.omega_scm) ** 2
        exponent = -((omega - self.omega_scm) ** 2) / (2 * delta_omega ** 2)
        gaussian = math.exp(exponent)
        sigma_n = self.sigma_0 * freq_ratio * gaussian

        return {
            "model": "base_CP2",
            "sigma_n": sigma_n,
            "sigma_0": self.sigma_0,
            "omega": omega,
            "omega_LENR": self.omega_scm,
            "freq_ratio_sq": freq_ratio,
            "gaussian": gaussian,
            "equation": "σ_n(ω) = σ₀ · (ω/ω_LENR)² · exp[-(ω-ω_LENR)²/(2Δω²)]",
        }

    # ── 1b. SCm-modulated cross-section ──────────────────────────────────

    def compute_scm_cross_section(self, omega: float,
                                  n: int = 13) -> Dict[str, Any]:
        """
        SCm Gaussian resonance with 26-level VDS enhancement:
          σ_n^SCm(ω,n) = σ₀ · exp[-(ω - ω_SCm)²/(2Γ²)] · (1 + [SSq]·n/26)

        Parameters:
          omega: angular frequency (rad/s)
          n: VDS level index (0–26), default 13 (midpoint)
        """
        # Gaussian resonance (pure, no pre-factor)
        exponent = -((omega - self.omega_scm) ** 2) / (2 * self.gamma ** 2)
        gaussian = math.exp(exponent)

        # VDS 26-level enhancement
        vds_factor = 1.0 + (self.ssq * n) / self.n_levels

        sigma_scm = self.sigma_0 * gaussian * vds_factor

        return {
            "model": "SCm_modulated",
            "sigma_n_scm": sigma_scm,
            "sigma_0": self.sigma_0,
            "omega": omega,
            "omega_SCm": self.omega_scm,
            "f_SCm_THz": self.omega_scm / (2 * PI) / 1e12,
            "Gamma": self.gamma,
            "Gamma_THz": self.gamma / (2 * PI) / 1e12,
            "gaussian": gaussian,
            "vds_level_n": n,
            "vds_factor": vds_factor,
            "SSq": self.ssq,
            "equation": "σ_n^SCm(ω,n) = σ₀ · exp[-(ω-ω_SCm)²/(2Γ²)] · (1+[SSq]·n/26)",
        }

    # ── 1c. Density-scaled cross-section with VDS ────────────────────────

    def compute_density_cross_section(self, rho: float) -> Dict[str, Any]:
        """
        Density-scaled with VDS power-law and exponential cutoff:
          σ_n^SCm(ρ) = σ₀ · (ρ_SCm/ρ_UA)^[SSq] · exp(-[SSq]·ρ/ρ_crit)

        The [SSq] exponent (0.57) replaces the linear form in CP2.
        """
        # VDS power-law ratio
        ratio = self.rho_scm / self.rho_ua
        vds_power = ratio ** self.ssq

        # Exponential cutoff at critical density
        exp_cutoff = math.exp(-self.ssq * rho / self.rho_crit)

        sigma_rho = self.sigma_0 * vds_power * exp_cutoff

        # Also compute base (linear) for comparison
        sigma_base = self.sigma_0 * (rho / self.rho_ua)

        return {
            "model": "SCm_density",
            "sigma_n_scm_rho": sigma_rho,
            "sigma_n_base_rho": sigma_base,
            "sigma_0": self.sigma_0,
            "rho": rho,
            "rho_SCm": self.rho_scm,
            "rho_UA": self.rho_ua,
            "rho_crit": self.rho_crit,
            "ratio_SCm_UA": ratio,
            "SSq_exponent": self.ssq,
            "vds_power": vds_power,
            "exp_cutoff": exp_cutoff,
            "equation": "σ_n^SCm(ρ) = σ₀ · (ρ_SCm/ρ_UA)^[SSq] · exp(-[SSq]·ρ/ρ_crit)",
            "equation_base": "σ_n(ρ) = σ₀ · (ρ/ρ₀)  [CP2 linear]",
        }

    # ── 1d. Buoyancy-coupled neutron production force ────────────────────

    def compute_f_neutron_scm(self, omega: float, n: int = 13,
                              N_n: float = 1e28,
                              phi_phonon: float = 1e16,
                              F_U_Bi: float = 1.1,
                              F_U: float = 1.0) -> Dict[str, Any]:
        """
        SCm buoyancy-coupled neutron production force:
          F_neutron^SCm = N_n · σ_n^SCm(ω,n) · Φ_phonon · (F_{U,Bi}/F_U - 1)

        Parameters:
          N_n: neutron number density in reaction zone (m⁻³)
          phi_phonon: 1.25 THz QCL phonon fluence (photons/m²/s)
          F_U_Bi: buoyancy force (reversed)
          F_U: unified gravitational force
        """
        # Get SCm cross-section
        cs = self.compute_scm_cross_section(omega, n)
        sigma_scm = cs["sigma_n_scm"]

        # Buoyancy reversal factor
        if F_U == 0:
            buoyancy_reversal = 0.0
        else:
            buoyancy_reversal = (F_U_Bi / F_U) - 1.0

        # Neutron production force
        F_neutron = N_n * sigma_scm * phi_phonon * buoyancy_reversal

        return {
            "model": "F_neutron_SCm",
            "F_neutron_scm": F_neutron,
            "N_n": N_n,
            "sigma_n_scm": sigma_scm,
            "phi_phonon": phi_phonon,
            "F_U_Bi": F_U_Bi,
            "F_U": F_U,
            "buoyancy_reversal_factor": buoyancy_reversal,
            "buoyancy_reversed": buoyancy_reversal > 0,
            "vds_level": n,
            "cross_section_detail": cs,
            "equation": "F_neutron^SCm = N_n · σ_n^SCm(ω,n) · Φ_phonon · (F_{U,Bi}/F_U - 1)",
        }

    # ── 1e. Frequency sweep ──────────────────────────────────────────────

    def frequency_sweep(self, omega_min: float = 0.5e12 * 2 * PI,
                        omega_max: float = 2.0e12 * 2 * PI,
                        n_points: int = 100,
                        vds_level: int = 13) -> Dict[str, Any]:
        """
        Sweep frequency and compare base vs SCm cross-sections.
        Returns arrays for plotting.
        """
        d_omega = (omega_max - omega_min) / (n_points - 1)
        omegas = [omega_min + i * d_omega for i in range(n_points)]
        f_thz = [w / (2 * PI) / 1e12 for w in omegas]
        sigma_base = []
        sigma_scm = []

        for w in omegas:
            base = self.compute_base_cross_section(w)
            scm = self.compute_scm_cross_section(w, vds_level)
            sigma_base.append(base["sigma_n"])
            sigma_scm.append(scm["sigma_n_scm"])

        # Peak values
        peak_base = max(sigma_base)
        peak_scm = max(sigma_scm)
        peak_idx = sigma_scm.index(peak_scm)

        return {
            "n_points": n_points,
            "vds_level": vds_level,
            "f_THz": f_thz,
            "sigma_base": sigma_base,
            "sigma_scm": sigma_scm,
            "peak_base": peak_base,
            "peak_scm": peak_scm,
            "peak_f_THz": f_thz[peak_idx],
            "scm_enhancement_at_peak": peak_scm / peak_base if peak_base > 0 else float("inf"),
        }

    # ── 1f. VDS level scan ───────────────────────────────────────────────

    def vds_level_scan(self, omega: float = OMEGA_SCM) -> Dict[str, Any]:
        """
        Compute cross-section at resonance for each VDS level n=0..26.
        Shows how the 26D enhancement builds up.
        """
        results = []
        for n in range(self.n_levels + 1):
            cs = self.compute_scm_cross_section(omega, n)
            results.append({
                "n": n,
                "vds_factor": cs["vds_factor"],
                "sigma_scm": cs["sigma_n_scm"],
            })

        return {
            "omega": omega,
            "f_THz": omega / (2 * PI) / 1e12,
            "SSq": self.ssq,
            "n_levels": self.n_levels,
            "levels": results,
            "sigma_n0": results[0]["sigma_scm"],
            "sigma_n26": results[-1]["sigma_scm"],
            "total_enhancement": results[-1]["vds_factor"],
        }


# ── §2  MAIN ──────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("Kozima SCm-Modulated Neutron-Drop Cross-Section Calculator")
    print("=" * 72)

    calc = KozimaSCmCrossSection()

    # 1. On-resonance comparison (base vs SCm)
    omega_res = OMEGA_SCM
    base = calc.compute_base_cross_section(omega_res)
    scm = calc.compute_scm_cross_section(omega_res, n=13)

    print(f"\n── On-Resonance (f = {F_LENR_THZ/1e12:.2f} THz) ──")
    print(f"  Base (CP2):        σ_n = {base['sigma_n']:.6e}")
    print(f"  SCm (n=13):        σ_n = {scm['sigma_n_scm']:.6e}")
    print(f"  VDS factor (n=13): {scm['vds_factor']:.4f}")
    print(f"  Enhancement:       {scm['sigma_n_scm']/base['sigma_n']:.4f}x" if base['sigma_n'] > 0 else "  Enhancement: ∞")

    # 2. Off-resonance (0.5 THz)
    omega_off = 2 * PI * 0.5e12
    base_off = calc.compute_base_cross_section(omega_off)
    scm_off = calc.compute_scm_cross_section(omega_off, n=13)
    print(f"\n── Off-Resonance (f = 0.50 THz) ──")
    print(f"  Base (CP2):        σ_n = {base_off['sigma_n']:.6e}")
    print(f"  SCm (n=13):        σ_n = {scm_off['sigma_n_scm']:.6e}")

    # 3. Density-scaled comparison
    rho_test = 1e-22
    dens = calc.compute_density_cross_section(rho_test)
    print(f"\n── Density-Scaled (ρ = {rho_test:.1e} kg/m³) ──")
    print(f"  Base (linear):     σ_n = {dens['sigma_n_base_rho']:.6e}")
    print(f"  SCm (VDS power):   σ_n = {dens['sigma_n_scm_rho']:.6e}")
    print(f"  VDS power factor:  (ρ_SCm/ρ_UA)^{dens['SSq_exponent']} = {dens['vds_power']:.6e}")

    # 4. F_neutron with buoyancy coupling
    fn = calc.compute_f_neutron_scm(
        omega=omega_res, n=13,
        N_n=1e28, phi_phonon=1e16,
        F_U_Bi=1.1, F_U=1.0,
    )
    print(f"\n── Buoyancy-Coupled Neutron Force ──")
    print(f"  F_neutron^SCm:     {fn['F_neutron_scm']:.6e} N")
    print(f"  Buoyancy reversal: {fn['buoyancy_reversal_factor']:.4f}")
    print(f"  Reversed:          {fn['buoyancy_reversed']}")

    # 5. VDS level scan at resonance
    vds = calc.vds_level_scan()
    print(f"\n── VDS Level Scan (n=0..26 at resonance) ──")
    print(f"  n=0:  σ = {vds['sigma_n0']:.6e}  (baseline)")
    print(f"  n=26: σ = {vds['sigma_n26']:.6e}  (max VDS)")
    print(f"  Total enhancement: {vds['total_enhancement']:.4f}x")

    # 6. Frequency sweep
    sweep = calc.frequency_sweep(vds_level=13)
    print(f"\n── Frequency Sweep (0.5–2.0 THz) ──")
    print(f"  Peak (base):  σ = {sweep['peak_base']:.6e}")
    print(f"  Peak (SCm):   σ = {sweep['peak_scm']:.6e}  at {sweep['peak_f_THz']:.4f} THz")
    print(f"  SCm/base:     {sweep['scm_enhancement_at_peak']:.4f}x")

    # Export
    export = {
        "calculator": "KozimaSCmCrossSection",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "on_resonance": {"base": base, "scm": scm},
        "off_resonance": {"base": base_off, "scm": scm_off},
        "density_scaled": dens,
        "f_neutron_scm": fn,
        "vds_scan": vds,
        "sweep_summary": {
            "peak_base": sweep["peak_base"],
            "peak_scm": sweep["peak_scm"],
            "peak_f_THz": sweep["peak_f_THz"],
            "enhancement": sweep["scm_enhancement_at_peak"],
        },
    }
    out_path = "kozima_scm_cross_section_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(export, f, indent=2, default=str)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print("COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
