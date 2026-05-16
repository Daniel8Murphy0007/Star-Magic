"""
Buoyancy Lagrangian Euler-Lagrange Equation of Motion Calculator

Session 204 | Daniel Murphy
PURPOSE: The Lagrangian variational machinery exists for the magnetic sector:
           δS/δA_SCm = 0  →  ∇×B_SCm = μ₀ J_SCm  (Ampère in SCm medium)
         [uqff_lagrangian_derivation.py L466, lagrangian_re_runner.py L252]

         And the buoyancy Lagrangian is defined:
           S_buoy = -∫d⁴x β_i Σ_i Ug_i Ω_g (M/d_g)(1+ε_sw ρ_sw)[UA]cos(πt_n)
         [uqff_lagrangian_derivation.py L533]

         But the Euler-Lagrange equation δS/δφ_buoyancy = 0 is NOT derived —
         only δS/δΩ_g = 0 (which gives Ubi_i directly) exists.

         This module derives the full variational equation for the buoyancy
         scalar field φ_buoy(x,t), treating it as a dynamical degree of freedom
         in the 9-sector UQFF Lagrangian.

PHYSICS:
  The buoyancy sector Lagrangian density is:
    L_buoy = ½(∂_μ φ)² - V(φ) + J_buoy·φ

  where:
    V(φ) = β_i Σ_i Ug_i · Ω_g · M/d_g · [UA] · φ²/2
    J_buoy = (F_{U,Bi}/F_U - 1) · ρ_SCm · c² (source term from vacuum reversal)

  The Euler-Lagrange EOM:
    □φ + ∂V/∂φ = J_buoy
    □φ + β_i Σ_i Ug_i · Ω_g · M/d_g · [UA] · φ = J_buoy

  This is a Klein-Gordon equation with effective mass:
    m_eff² = β_i Σ_i Ug_i · Ω_g · M/(d_g · c² · ħ²) · [UA]

  Static solution (∂_t φ = 0):
    φ_static(r) = (J_buoy / m_eff²) · [1 - exp(-m_eff·r)]

  This gives the buoyancy reversal range — the distance over which
  SCm vacuum displacement produces observable weight anomalies.

ARCHITECTURE: Pure calculator. Receives dataset with Ug values, returns EOM.
"""

import json
import math
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
K_B = 1.381e-23

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99
U_UA = 1e-4

# Vacuum
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36

# Solar
M_SUN = 1.989e30
PC_M = 3.086e16  # parsec in meters


# ── §1  BUOYANCY LAGRANGIAN EOM ───────────────────────────────────────────

class BuoyancyLagrangianEOM:
    """
    Euler-Lagrange equation for the buoyancy scalar field φ_buoy(x,t).

    Lagrangian density:
      L = ½(∂_μ φ)² - V(φ) + J_buoy·φ

    Potential:
      V(φ) = ½ m_eff² φ²
      m_eff² = β_i · Σ Ug_i · Ω_g · M / (d_g · c² · ħ²) · [UA]

    EOM (Klein-Gordon):
      □φ + m_eff² φ = J_buoy    (□ = ∂²/∂t² - ∇²)

    Source term:
      J_buoy = (F_{U,Bi}/F_U - 1) · ρ_SCm · c²
    """

    def __init__(self,
                 beta_i: float = BETA_I,
                 u_ua: float = U_UA,
                 rho_scm: float = RHO_SCM):
        self.beta_i = beta_i
        self.u_ua = u_ua
        self.rho_scm = rho_scm

    # ── 1a. Effective mass ───────────────────────────────────────────────

    def compute_effective_mass(self,
                               Ug_values: List[float],
                               Omega_g: float = 1.0,
                               M_kg: float = 4.3e6 * M_SUN,
                               d_g_m: float = 2.44e20) -> Dict[str, Any]:
        """
        Compute the effective boson mass for the buoyancy field.

        m_eff² = β_i · Σ|Ug_i| · Ω_g · M / (d_g · c² · ħ²) · [UA]

        This determines the range of the buoyancy force:
          λ_buoy = ħ / (m_eff · c)  (Compton wavelength)
        """
        Ug_sum = sum(abs(u) for u in Ug_values)

        # Effective mass squared (in natural units, then convert)
        m_eff_sq_nat = (self.beta_i * Ug_sum * Omega_g * M_kg
                        / (d_g_m * C**2 * HBAR**2) * self.u_ua)

        m_eff = math.sqrt(abs(m_eff_sq_nat)) if m_eff_sq_nat != 0 else 0

        # Compton wavelength (buoyancy range)
        if m_eff > 0:
            lambda_buoy = HBAR / (m_eff * C)
        else:
            lambda_buoy = float("inf")

        # Energy scale
        E_eff = m_eff * C**2 if m_eff > 0 else 0

        return {
            "m_eff_sq": m_eff_sq_nat,
            "m_eff_kg": m_eff,
            "m_eff_eV": E_eff / 1.602e-19 if E_eff > 0 else 0,
            "lambda_buoy_m": lambda_buoy,
            "lambda_buoy_pc": lambda_buoy / PC_M if lambda_buoy < 1e30 else float("inf"),
            "Ug_sum": Ug_sum,
            "n_Ug_terms": len(Ug_values),
            "beta_i": self.beta_i,
            "Omega_g": Omega_g,
            "M_kg": M_kg,
            "d_g_m": d_g_m,
            "U_UA": self.u_ua,
            "equation": "m_eff² = β_i · Σ|Ug_i| · Ω_g · M/(d_g·c²·ħ²) · [UA]",
        }

    # ── 1b. Source term ──────────────────────────────────────────────────

    @staticmethod
    def near_horizon_heaviside(r: float,
                               r_s: float,
                               shell_width: float = 0.01,
                               amplification: float = 1.0e13) -> float:
        """
        Near-horizon Heaviside phase-transition amplifier f_H(r/r_s).

        The SCm vacuum undergoes a phase transition in a thin shell just
        outside the Schwarzschild horizon, amplifying the buoyancy source
        current by (1 + amplification · f_H) for r_s < r < r_s(1+shell_width).

        f_H = 1   for r_s ≤ r ≤ r_s(1+shell_width)
        f_H = 0   elsewhere

        Returns the multiplicative amplifier (1 + amplification · f_H).
        Item 3 / first_principles_derivation.py — see PAPER_1171 §6.
        """
        if r < r_s or r > r_s * (1.0 + shell_width):
            f_H = 0.0
        else:
            f_H = 1.0
        return 1.0 + amplification * f_H

    def compute_source_term(self,
                            F_U_Bi: float = 1.1,
                            F_U: float = 1.0,
                            r: Optional[float] = None,
                            r_s: Optional[float] = None,
                            shell_width: float = 0.01,
                            amplification: float = 1.0e13) -> Dict[str, Any]:
        """
        Buoyancy source current with optional near-horizon amplifier:
          J_buoy = (F_{U,Bi}/F_U - 1) · ρ_SCm · c² · (1 + 10¹³·f_H(r/r_s))

        f_H is the Heaviside near-horizon phase-transition function.
        When r and r_s are supplied, the amplifier is applied; otherwise
        the bare current is returned (backward compatible).
        """
        if F_U == 0:
            reversal = 0.0
        else:
            reversal = (F_U_Bi / F_U) - 1.0

        if r is not None and r_s is not None:
            amp = self.near_horizon_heaviside(r, r_s, shell_width, amplification)
        else:
            amp = 1.0

        J_buoy = reversal * self.rho_scm * C**2 * amp

        return {
            "J_buoy": J_buoy,
            "reversal_factor": reversal,
            "near_horizon_amplifier": amp,
            "buoyancy_reversed": reversal > 0,
            "F_U_Bi": F_U_Bi,
            "F_U": F_U,
            "rho_SCm": self.rho_scm,
            "equation": "J_buoy = (F_{U,Bi}/F_U - 1) · ρ_SCm · c² · (1 + 10¹³·f_H)",
        }

    # ── 1c. Full EOM derivation ──────────────────────────────────────────

    def derive_eom(self,
                   Ug_values: List[float],
                   Omega_g: float = 1.0,
                   M_kg: float = 4.3e6 * M_SUN,
                   d_g_m: float = 2.44e20,
                   F_U_Bi: float = 1.1,
                   F_U: float = 1.0,
                   epsilon_sw: float = 1e-5,
                   rho_sw: float = 1e-20,
                   t_n: float = 0.0) -> Dict[str, Any]:
        """
        Full Euler-Lagrange derivation for the buoyancy sector.

        Returns the complete derivation chain from Lagrangian → EOM → solution.
        """
        # Wind factor
        wind = 1.0 + epsilon_sw * rho_sw
        cos_tn = math.cos(PI * t_n)

        # Ubi values (existing CP2/Lagrangian formula)
        Ubi_values = []
        for Ug in Ug_values:
            ubi = -self.beta_i * Ug * Omega_g * M_kg / d_g_m * wind * self.u_ua * cos_tn
            Ubi_values.append(ubi)
        F_U_Bi_total = sum(Ubi_values)

        # Effective mass
        mass = self.compute_effective_mass(Ug_values, Omega_g, M_kg, d_g_m)

        # Source term
        source = self.compute_source_term(F_U_Bi, F_U)

        # Static solution
        m_eff = mass["m_eff_kg"]
        J_buoy = source["J_buoy"]

        if m_eff > 0 and mass["m_eff_sq"] != 0:
            phi_static_amplitude = J_buoy / mass["m_eff_sq"]
        else:
            phi_static_amplitude = 0.0

        # Derivation chain (human-readable)
        chain = [
            "══════════════════════════════════════════════════════════════",
            "BUOYANCY LAGRANGIAN EULER-LAGRANGE DERIVATION",
            "══════════════════════════════════════════════════════════════",
            "",
            "§1. Buoyancy Sector Lagrangian Density:",
            "  L_buoy = ½(∂_μ φ_buoy)² − V(φ_buoy) + J_buoy · φ_buoy",
            "",
            "§2. Potential (from 4-layer Ug summation):",
            "  V(φ) = ½ m_eff² φ²",
            f"  m_eff² = β_i·Σ|Ug_i|·Ω_g·M/(d_g·c²·ħ²)·[UA]",
            f"         = {self.beta_i}·{mass['Ug_sum']:.4e}·{Omega_g}·{M_kg:.4e}/({d_g_m:.4e}·c²·ħ²)·{self.u_ua}",
            f"         = {mass['m_eff_sq']:.6e}  kg²",
            f"  m_eff  = {m_eff:.6e}  kg  ({mass['m_eff_eV']:.4e} eV)",
            "",
            "§3. Source Current (vacuum displacement → buoyancy reversal):",
            f"  J_buoy = (F_{{U,Bi}}/F_U − 1) · ρ_SCm · c²",
            f"         = ({F_U_Bi}/{F_U} − 1) · {self.rho_scm:.4e} · c²",
            f"         = {J_buoy:.6e}  kg/m³·(m/s)²",
            f"  Buoyancy reversed: {source['buoyancy_reversed']}",
            "",
            "§4. Euler-Lagrange Equation:",
            "  δS/δφ_buoy = 0",
            "  ∂L/∂φ − ∂_μ(∂L/∂(∂_μφ)) = 0",
            "  −m_eff²·φ + J_buoy − □φ = 0",
            "",
            "  ══> □φ + m_eff²·φ = J_buoy  (Klein-Gordon with source)",
            "",
            "§5. Static Solution (∂_t φ = 0):",
            "  ∇²φ = m_eff²·φ − J_buoy",
            "  φ_static(r) = (J_buoy/m_eff²)·[1 − exp(−m_eff·r)]",
            f"  φ_∞ = J_buoy/m_eff² = {phi_static_amplitude:.6e}",
            "",
            "§6. Buoyancy Range (Compton wavelength):",
            f"  λ_buoy = ħ/(m_eff·c) = {mass['lambda_buoy_m']:.6e} m",
            f"         = {mass['lambda_buoy_pc']:.6e} pc" if mass['lambda_buoy_pc'] < 1e20 else "",
            "  This is the distance over which SCm vacuum displacement",
            "  produces measurable buoyancy anomalies.",
            "",
            "§7. Consistency Check (existing Ubi values from δS/δΩ_g = 0):",
        ]

        for i, (ug, ubi) in enumerate(zip(Ug_values, Ubi_values)):
            chain.append(f"  Ug{i+1} = {ug:.4e}  →  Ubi{i+1} = {ubi:.4e}")
        chain.append(f"  F_{{U,Bi}} (total) = {F_U_Bi_total:.4e}")
        chain.append("")
        chain.append("§8. Connection to Existing Variational Equations:")
        chain.append("  ✓ δS/δA_SCm = 0  →  ∇×B_SCm = μ₀ J_SCm  (magnetic, L466)")
        chain.append("  ✓ δS/δΩ_g  = 0  →  Ubi_i (buoyancy forces, L537)")
        chain.append("  ✓ δS/δφ̂    = 0  →  Um (helical string magnetism, L543)")
        chain.append("  ★ δS/δφ_buoy = 0  →  □φ + m_eff²φ = J_buoy  [THIS MODULE]")
        chain.append("══════════════════════════════════════════════════════════════")

        return {
            "lagrangian": "L_buoy = ½(∂_μ φ)² − ½m_eff²φ² + J_buoy·φ",
            "eom": "□φ + m_eff²·φ = J_buoy",
            "eom_type": "Klein-Gordon with source",
            "effective_mass": mass,
            "source_term": source,
            "phi_static_amplitude": phi_static_amplitude,
            "Ubi_values": Ubi_values,
            "F_U_Bi_total": F_U_Bi_total,
            "wind_factor": wind,
            "cos_pi_tn": cos_tn,
            "derivation_chain": chain,
            "variational_equations": {
                "magnetic": "δS/δA_SCm = 0 → ∇×B_SCm = μ₀ J_SCm",
                "buoyancy_Omega": "δS/δΩ_g = 0 → Ubi_i = -β_i Ug_i Ω_g M/d_g [UA]",
                "magnetism_hat": "δS/δφ̂ = 0 → Um = Σ μ_j/r_j (1-e^{-γt}) φ̂ N P E",
                "buoyancy_field": "δS/δφ_buoy = 0 → □φ + m_eff²φ = J_buoy  [NEW]",
            },
        }

    # ── 1d. Sgr A* default system ────────────────────────────────────────

    def derive_sgr_a_star(self) -> Dict[str, Any]:
        """
        Derive buoyancy EOM for Sgr A* (default astrophysical system).
        Uses canonical Ug placeholder values.
        """
        return self.derive_eom(
            Ug_values=[1e-10, 1e-10, 1e-10, 1e-10],
            Omega_g=1.0,
            M_kg=4.3e6 * M_SUN,
            d_g_m=2.44e20,
            F_U_Bi=1.1,
            F_U=1.0,
        )

    # ── 1e. Lab-scale system ─────────────────────────────────────────────

    def derive_lab_scale(self,
                         M_sample_kg: float = 0.5,
                         r_m: float = 0.1,
                         F_U_Bi: float = 1.05,
                         F_U: float = 1.0) -> Dict[str, Any]:
        """
        Derive buoyancy EOM for a lab-scale LENR system.
        """
        # Lab Ug estimates (gravitational acceleration components)
        g_local = G * 5.972e24 / (6.371e6) ** 2  # Earth surface
        Ug_lab = [g_local * 1e-6, g_local * 1e-7, g_local * 1e-8, g_local * 1e-9]

        return self.derive_eom(
            Ug_values=Ug_lab,
            Omega_g=1.0,
            M_kg=M_sample_kg,
            d_g_m=r_m,
            F_U_Bi=F_U_Bi,
            F_U=F_U,
        )


# ── §2  MAIN ──────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("Buoyancy Lagrangian Euler-Lagrange EOM Calculator")
    print("=" * 72)

    calc = BuoyancyLagrangianEOM()

    # 1. Sgr A* system
    print(f"\n── Sgr A* System (astrophysical scale) ──")
    sgr = calc.derive_sgr_a_star()

    for line in sgr["derivation_chain"]:
        print(f"  {line}")

    print(f"\n  EOM: {sgr['eom']}")
    print(f"  m_eff = {sgr['effective_mass']['m_eff_kg']:.6e} kg")
    print(f"  λ_buoy = {sgr['effective_mass']['lambda_buoy_m']:.6e} m")
    print(f"  φ_∞ = {sgr['phi_static_amplitude']:.6e}")
    print(f"  J_buoy = {sgr['source_term']['J_buoy']:.6e}")

    # 2. Lab-scale system
    print(f"\n── Lab-Scale System (Colman-Gillespie, 0.5 kg) ──")
    lab = calc.derive_lab_scale(M_sample_kg=0.5, r_m=0.1, F_U_Bi=1.05, F_U=1.0)

    print(f"  EOM: {lab['eom']}")
    print(f"  m_eff = {lab['effective_mass']['m_eff_kg']:.6e} kg")
    print(f"  λ_buoy = {lab['effective_mass']['lambda_buoy_m']:.6e} m")
    print(f"  φ_∞ = {lab['phi_static_amplitude']:.6e}")
    print(f"  J_buoy = {lab['source_term']['J_buoy']:.6e}")
    print(f"  Buoyancy reversed: {lab['source_term']['buoyancy_reversed']}")

    # 3. Variational equations summary
    print(f"\n── Complete Variational Equations (9-Sector Lagrangian) ──")
    for key, eq in sgr["variational_equations"].items():
        marker = "★ NEW" if "NEW" in eq else "  ✓  "
        print(f"  {marker} {eq}")

    # Export
    export = {
        "calculator": "BuoyancyLagrangianEOM",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "sgr_a_star": {
            "eom": sgr["eom"],
            "effective_mass": sgr["effective_mass"],
            "source_term": sgr["source_term"],
            "phi_static": sgr["phi_static_amplitude"],
            "Ubi_values": sgr["Ubi_values"],
            "variational_equations": sgr["variational_equations"],
        },
        "lab_scale": {
            "eom": lab["eom"],
            "effective_mass": lab["effective_mass"],
            "source_term": lab["source_term"],
            "phi_static": lab["phi_static_amplitude"],
        },
    }
    out_path = "buoyancy_lagrangian_eom_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(export, f, indent=2, default=str)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print("COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
