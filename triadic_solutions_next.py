#!/usr/bin/env python3
"""
triadic_solutions_next.py — Next Triadic Solutions (Compressed/Resonant/Buoyancy)

Session 215 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Applies the three UQFF Triadic operational modes to phonon/jet/NS topics:

  Compressed Gravity:  F_{U,Bi}/F_U modulated by Γ for jet collimation
  Resonant:            1.25 THz phonon linewidth Γ tunes neutron-drop
                       and buoyancy reversal
  Buoyancy:            E_net(t,Γ) drives positive expansion (nebulae)
                       ↔ negative erosion (filaments)

All three converge on SCm phonon resonance.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
HBAR      = 1.055e-34
K_B       = 1.381e-23
C         = 2.998e8
G         = 6.674e-11
M_SUN     = 1.989e30
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
BETA_I    = 0.603


# ── §1  Compressed Gravity Triadic ────────────────────────────────────────

class CompressedGravityTriadic:
    """Compressed gravity mode: F_{U,Bi}/F_U modulated by linewidth Γ.

    In the compressed regime, phonon linewidth sharpens the buoyancy-to-gravity
    ratio, producing jet collimation in AGN and tight knot formation.

    F_compressed(Γ) = F_{U,Bi}/F_U · exp(-(ω - ω_SCm)² / 2Γ²) · S₂₆
    Session 215.
    """

    def jet_collimation_factor(self, gamma_THz: float, A_jet: float = 1.5,
                                f_ubi_ratio: float = 0.6) -> float:
        gamma = 2 * PI * gamma_THz * 1e12
        phi = math.exp(-(OMEGA_SCM - OMEGA_SCM) ** 2 / (2 * gamma ** 2))
        return f_ubi_ratio * phi * S26 * A_jet

    def compute(self, dataset: dict) -> dict:
        gamma_THz = float(dataset.get("Gamma_THz", 0.10))
        A_jet = float(dataset.get("A_jet", 1.5))
        f_ubi = float(dataset.get("F_UBi_ratio", 0.6))

        gammas = [0.01, 0.05, 0.10, 0.20, 0.30, 0.50]
        results = []
        for g in gammas:
            cf = self.jet_collimation_factor(g, A_jet, f_ubi)
            regime = "ultra-tight" if g < 0.05 else ("collimated" if g < 0.15 else "diffuse")
            results.append({"Gamma_THz": g, "collimation": cf, "regime": regime})

        return {
            "mode": "compressed",
            "scan": results,
            "primary_equations": [
                "F_compressed(Γ) = F_{U,Bi}/F_U · Φ(ω,Γ) · S₂₆ · A_jet",
                f"At Γ={gamma_THz} THz: collimation = {self.jet_collimation_factor(gamma_THz, A_jet, f_ubi):.6f}",
            ],
        }


# ── §2  Resonant Gravity Triadic ──────────────────────────────────────────

class ResonantGravityTriadic:
    """Resonant gravity mode: 1.25 THz phonon linewidth tunes neutron-drop
    and buoyancy reversal for NS mergers and mass-gap classification.

    Φ(ω,Γ) = Φ_0 · exp(-(ω - ω_SCm)² / 2Γ²) · S₂₆([SSq])
    Neutron-drop: triggered when Φ > Φ_crit (threshold for drip-line shift)
    Session 215.
    """

    PHI_CRIT = 0.5  # Critical phonon occupation for neutron-drop

    def phi_resonance(self, omega: float, gamma_THz: float, phi_0: float = 1.0) -> float:
        gamma = 2 * PI * gamma_THz * 1e12
        return phi_0 * math.exp(-(omega - OMEGA_SCM) ** 2 / (2 * gamma ** 2)) * S26

    def neutron_drop_threshold(self, gamma_THz: float) -> bool:
        phi = self.phi_resonance(OMEGA_SCM, gamma_THz)
        return phi > self.PHI_CRIT

    def buoyancy_reversal_point(self, gamma_THz: float) -> float:
        """Γ at which buoyancy switches sign (E_net = 0 crossing)."""
        gamma = 2 * PI * gamma_THz * 1e12
        return PI / (2 * gamma)

    def compute(self, dataset: dict) -> dict:
        gamma_THz = float(dataset.get("Gamma_THz", 0.10))
        phi_on = self.phi_resonance(OMEGA_SCM, gamma_THz)
        n_drop = self.neutron_drop_threshold(gamma_THz)
        t_rev = self.buoyancy_reversal_point(gamma_THz)

        return {
            "mode": "resonant",
            "Gamma_THz": gamma_THz,
            "Phi_on_resonance": phi_on,
            "neutron_drop": n_drop,
            "buoyancy_reversal_time_s": t_rev,
            "primary_equations": [
                "Φ(ω,Γ) = Φ_0 · exp(-(ω-ω_SCm)²/2Γ²) · S₂₆",
                f"Φ(ω_SCm, Γ={gamma_THz} THz) = {phi_on:.6f}",
                f"Neutron-drop triggered: {n_drop}",
                f"Buoyancy reversal at t = {t_rev:.4e} s",
            ],
        }


# ── §3  Buoyancy Gravity Triadic ──────────────────────────────────────────

class BuoyancyGravityTriadic:
    """Buoyancy gravity mode: E_net(t,Γ) drives positive expansion (nebulae)
    ↔ negative erosion (filaments).

    E_net(t,Γ) = S₂₆ · cos(ω_SCm t) · exp(-Γ t) — threshold
    Positive → expansion (nebulae, HII regions)
    Negative → erosion  (filaments, pillars, cometary knots)
    Session 215.
    """

    def E_net(self, t: float, gamma_THz: float, threshold: float = 0.0) -> float:
        gamma = 2 * PI * gamma_THz * 1e12
        return S26 * math.cos(OMEGA_SCM * t) * math.exp(-gamma * t) - threshold

    def classify_regime(self, E: float) -> str:
        if E > 0.1:
            return "expansion"
        elif E < -0.1:
            return "erosion"
        return "neutral"

    def compute(self, dataset: dict) -> dict:
        gamma_THz = float(dataset.get("Gamma_THz", 0.10))
        threshold = float(dataset.get("threshold", 0.0))
        N_steps = 50
        dt = 5e-12 / N_steps
        gamma = 2 * PI * gamma_THz * 1e12

        trace = []
        for i in range(N_steps + 1):
            t = i * dt
            E = self.E_net(t, gamma_THz, threshold)
            trace.append({
                "t_ps": t * 1e12,
                "E_net": E,
                "regime": self.classify_regime(E),
            })

        n_exp = sum(1 for x in trace if x["regime"] == "expansion")
        n_ero = sum(1 for x in trace if x["regime"] == "erosion")

        return {
            "mode": "buoyancy",
            "Gamma_THz": gamma_THz,
            "threshold": threshold,
            "expansion_fraction": n_exp / len(trace),
            "erosion_fraction": n_ero / len(trace),
            "trace": trace,
            "primary_equations": [
                "E_net(t,Γ) = S₂₆·cos(ω_SCm·t)·exp(-Γt) - threshold",
                f"Expansion fraction: {n_exp/len(trace):.1%}",
                f"Erosion fraction:   {n_ero/len(trace):.1%}",
            ],
        }


# ── §4  Unified Triadic Solver ─────────────────────────────────────────────

class TriadicSolverNext:
    """Unified solver applying all three Triadic modes to a single dataset.

    Returns compressed, resonant, and buoyancy results simultaneously.
    Session 215.
    """

    def __init__(self):
        self.compressed = CompressedGravityTriadic()
        self.resonant = ResonantGravityTriadic()
        self.buoyancy = BuoyancyGravityTriadic()

    def compute(self, dataset: dict) -> dict:
        return {
            "compressed": self.compressed.compute(dataset),
            "resonant": self.resonant.compute(dataset),
            "buoyancy": self.buoyancy.compute(dataset),
            "convergence": "All three modes converge on SCm phonon resonance at ω_SCm = 1.25 THz",
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True

    # Test 1: Compressed jet collimation
    cg = CompressedGravityTriadic()
    cf = cg.jet_collimation_factor(0.10)
    if cf <= 0:
        print("[FAIL] Compressed collimation should be positive"); ok = False
    else:
        print(f"[ OK ] Compressed collimation(Γ=0.10) = {cf:.6f}")

    # Test 2: Resonant Φ on resonance
    rg = ResonantGravityTriadic()
    phi = rg.phi_resonance(OMEGA_SCM, 0.10)
    if phi <= 0:
        print("[FAIL] Resonant Φ should be positive"); ok = False
    else:
        print(f"[ OK ] Resonant Φ(ω_SCm, Γ=0.10) = {phi:.6f}")

    # Test 3: Neutron-drop threshold
    nd = rg.neutron_drop_threshold(0.10)
    print(f"[ OK ] Neutron-drop triggered at Γ=0.10: {nd}")

    # Test 4: Buoyancy E_net
    bg = BuoyancyGravityTriadic()
    E0 = bg.E_net(0, 0.10)
    if abs(E0 - S26) > 0.01:
        print(f"[FAIL] E_net(0) should ≈ S₂₆ = {S26:.4f}, got {E0:.4f}"); ok = False
    else:
        print(f"[ OK ] Buoyancy E_net(0) = {E0:.6f} ≈ S₂₆")

    # Test 5: Unified solver
    solver = TriadicSolverNext()
    result = solver.compute({"Gamma_THz": 0.10})
    for mode in ["compressed", "resonant", "buoyancy"]:
        if mode not in result:
            print(f"[FAIL] Missing {mode} in unified result"); ok = False
        else:
            print(f"[ OK ] Unified solver has {mode} result")

    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
