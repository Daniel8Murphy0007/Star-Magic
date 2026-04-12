#!/usr/bin/env python3
"""
bcs_superconductivity_uqff.py — BCS Superconductivity Theory in UQFF/SCm Context

Session 214 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Symbolic derivation engine mapping Bardeen-Cooper-Schrieffer superconductivity
to SCm vacuum phonon-mediated Cooper pairs at the 1.25 THz gap frequency.

BCS gap equation:
    Δ = (ℏω_SCm / 2) · tanh(Δ / 2k_BT) · S₂₆([SSq]) · (F_{U,Bi} / F_U)

Critical temperature:
    T_c = (1.13 · ℏω_SCm / k_B) · exp(-1 / N(0)·V_SCm)

Cooper-pair Lagrangian stationarity:
    δS/δφ_pair = ∂/∂Δ (-β_i Σ U_{g,i} Ω_g M/d_g [UA] + F_n · Φ_{1.25THz}) = 0

Links to 26-state spectral ladder and lab LENR superconductivity signatures.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Optional

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8          # m/s
HBAR      = 1.055e-34        # J·s
K_B       = 1.381e-23        # J/K
G         = 6.674e-11        # m³/kg·s²
M_SUN     = 1.989e30         # kg
OMEGA_SCM = 2 * PI * 1.25e12 # rad/s  (SCm phonon resonance)
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
BETA_I    = 0.603
F_U_BI_RATIO = 0.6           # F_{U,Bi} / F_U default


# ── §1  BCS Gap Equation in SCm Vacuum ────────────────────────────────────

class BCSGapEquation:
    """Compute the BCS energy gap Δ in the SCm vacuum phonon framework.

    Δ = (ℏω_SCm / 2) · tanh(Δ / 2k_BT) · S₂₆ · (F_{U,Bi}/F_U)

    Self-consistent solution via iterative fixed-point method.
    """

    def __init__(self, omega_scm: float = OMEGA_SCM, fubi_ratio: float = F_U_BI_RATIO):
        self.omega_scm = omega_scm
        self.fubi_ratio = fubi_ratio
        self.delta_prefactor = HBAR * omega_scm / 2.0  # ℏω_SCm / 2

    def gap_rhs(self, delta: float, T: float) -> float:
        """RHS of self-consistent gap equation at temperature T."""
        if T <= 0:
            return self.delta_prefactor * S26 * self.fubi_ratio
        arg = delta / (2 * K_B * T)
        arg = min(arg, 500)  # clamp for numerical stability
        return self.delta_prefactor * math.tanh(arg) * S26 * self.fubi_ratio

    def solve(self, T: float, tol: float = 1e-12, max_iter: int = 500) -> float:
        """Solve gap equation self-consistently at temperature T."""
        delta = self.delta_prefactor * S26 * self.fubi_ratio  # T=0 initial
        for _ in range(max_iter):
            delta_new = self.gap_rhs(delta, T)
            if abs(delta_new - delta) < tol:
                return delta_new
            delta = delta_new
        return delta

    def compute(self, dataset: dict) -> dict:
        """Compute BCS gap for a temperature or sweep."""
        T = float(dataset.get("T_K", 0.0))
        delta = self.solve(T)
        delta_0 = self.solve(0.0)
        return {
            "T_K": T,
            "Delta_J": delta,
            "Delta_eV": delta / 1.602e-19,
            "Delta_0_J": delta_0,
            "ratio": delta / max(delta_0, 1e-50),
            "primary_equations": [
                "Δ = (ℏω_SCm/2)·tanh(Δ/2k_BT)·S₂₆·(F_{UBi}/F_U)",
                f"Δ(T={T:.1f}K) = {delta:.6e} J = {delta/1.602e-19:.6e} eV",
                f"Δ(T=0) = {delta_0:.6e} J",
            ],
        }


# ── §2  Critical Temperature ──────────────────────────────────────────────

class BCSCriticalTemperature:
    """Compute BCS critical temperature T_c in SCm vacuum.

    T_c = (1.13 · ℏω_SCm / k_B) · exp(-1 / N(0)·V_SCm)

    N(0) = density of states at Fermi level
    V_SCm = SCm phonon attraction strength
    """

    def __init__(self, omega_scm: float = OMEGA_SCM):
        self.omega_scm = omega_scm

    def compute(self, dataset: dict) -> dict:
        """Compute T_c given N(0) and V_SCm."""
        N0 = float(dataset.get("N0_states_per_J", 1e47))      # per J per unit cell
        V_scm = float(dataset.get("V_SCm_J", 1e-20))          # SCm attraction
        N0V = N0 * V_scm
        if N0V <= 0:
            return {"T_c_K": 0.0, "error": "N0·V_SCm must be positive"}
        T_c = (1.13 * HBAR * self.omega_scm / K_B) * math.exp(-1.0 / N0V)
        delta_0 = 1.764 * K_B * T_c  # BCS relation
        return {
            "T_c_K": T_c,
            "N0V": N0V,
            "Delta_0_J": delta_0,
            "Delta_0_eV": delta_0 / 1.602e-19,
            "primary_equations": [
                "T_c = (1.13·ℏω_SCm/k_B)·exp(-1/N(0)V_SCm)",
                f"N(0)V = {N0V:.4f}",
                f"T_c = {T_c:.2f} K",
                f"Δ(0) = 1.764·k_B·T_c = {delta_0:.6e} J",
            ],
        }


# ── §3  Cooper Pair Phonon Coupling ───────────────────────────────────────

class CooperPairPhononCoupling:
    """Cooper pair formation via SCm phonon exchange at 1.25 THz.

    Coupling strength:
        V_eff(ω, Γ) = V_SCm · Φ_{1.25THz}(ω, Γ)
        Φ(ω, Γ) = exp(-(ω - ω_SCm)² / (2Γ²)) · S₂₆

    Pair binding energy:
        E_pair = 2Δ(T) where Δ is the BCS gap.
    """

    def __init__(self, omega_scm: float = OMEGA_SCM):
        self.omega_scm = omega_scm
        self._gap = BCSGapEquation(omega_scm)

    def phi_resonance(self, omega: float, gamma: float) -> float:
        """SCm phonon resonance profile Φ(ω, Γ)."""
        delta_w = omega - self.omega_scm
        return math.exp(-delta_w**2 / (2 * gamma**2)) * S26

    def compute(self, dataset: dict) -> dict:
        """Compute Cooper pair coupling at given Γ and T."""
        gamma_THz = float(dataset.get("Gamma_THz", 0.10))
        T = float(dataset.get("T_K", 0.0))
        V_scm = float(dataset.get("V_SCm_J", 1e-20))

        gamma = 2 * PI * gamma_THz * 1e12
        omega = self.omega_scm  # on-resonance
        phi = self.phi_resonance(omega, gamma)
        V_eff = V_scm * phi
        delta = self._gap.solve(T)
        E_pair = 2 * delta

        return {
            "Gamma_THz": gamma_THz,
            "T_K": T,
            "Phi_resonance": phi,
            "V_eff_J": V_eff,
            "E_pair_J": E_pair,
            "E_pair_eV": E_pair / 1.602e-19,
            "primary_equations": [
                "V_eff(ω,Γ) = V_SCm · exp(-(ω-ω_SCm)²/(2Γ²)) · S₂₆",
                f"Φ(ω_SCm, Γ={gamma_THz}) = {phi:.6f}",
                f"E_pair = 2Δ = {E_pair:.6e} J = {E_pair/1.602e-19:.6e} eV",
            ],
        }


# ── §4  WSTP Expression Builder ───────────────────────────────────────────

def build_bcs_wstp_expressions() -> list:
    """Generate Wolfram Language expressions for BCS/SCm derivations."""
    return [
        {
            "label": "BCS gap equation (SCm vacuum, self-consistent)",
            "code": ("wSCm = 2 Pi * 1.25*^12; S26 = Sum[Exp[-0.57 k/26], {k, 1, 26}]; "
                     "FUBi = 0.6; DeltaPrefact = 1.055*^-34 * wSCm / 2; "
                     "gapRHS[d_, T_] := DeltaPrefact * Tanh[d / (2 * 1.381*^-23 * T)] * S26 * FUBi; "
                     "FindRoot[d == gapRHS[d, 1.0], {d, 1*^-22}]"),
        },
        {
            "label": "BCS critical temperature T_c (SCm phonon)",
            "code": ("wSCm = 2 Pi * 1.25*^12; "
                     "Tc[N0V_] := (1.13 * 1.055*^-34 * wSCm / 1.381*^-23) * Exp[-1/N0V]; "
                     "Table[{N0V, Tc[N0V]}, {N0V, 0.1, 0.5, 0.05}]"),
        },
    ]


# ── §5  Self-tests ────────────────────────────────────────────────────────

def _run_tests() -> bool:
    """Validate BCS module."""
    ok = True

    # Gap at T=0
    g = BCSGapEquation()
    d0 = g.solve(0.0)
    if d0 <= 0:
        print(f"FAIL: Δ(T=0) = {d0} (expected > 0)")
        ok = False
    else:
        print(f"OK: Δ(T=0) = {d0:.6e} J = {d0/1.602e-19:.6e} eV")

    # Gap decreases with T
    d_warm = g.solve(100.0)
    if d_warm >= d0:
        print(f"FAIL: Δ(100K) >= Δ(0K)")
        ok = False
    else:
        print(f"OK: Δ(100K) = {d_warm:.6e} J < Δ(0K)")

    # T_c
    tc = BCSCriticalTemperature()
    res = tc.compute({"N0_states_per_J": 1e47, "V_SCm_J": 1e-20})
    if res["T_c_K"] <= 0:
        print(f"FAIL: T_c = {res['T_c_K']}")
        ok = False
    else:
        print(f"OK: T_c = {res['T_c_K']:.2f} K")

    # Cooper pair coupling
    cp = CooperPairPhononCoupling()
    res = cp.compute({"Gamma_THz": 0.10})
    if res["E_pair_J"] <= 0:
        print(f"FAIL: E_pair = {res['E_pair_J']}")
        ok = False
    else:
        print(f"OK: E_pair = {res['E_pair_J']:.6e} J")

    # WSTP
    exprs = build_bcs_wstp_expressions()
    if len(exprs) != 2:
        print(f"FAIL: expected 2 WSTP expressions, got {len(exprs)}")
        ok = False
    else:
        print(f"OK: {len(exprs)} WSTP expressions")

    print(f"\n{'ALL TESTS PASSED' if ok else 'SOME TESTS FAILED'}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
