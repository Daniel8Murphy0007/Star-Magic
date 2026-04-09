"""
SCm Activation Function Calculator — Linear vs Gaussian B-Field Suppression

Session 204 | Daniel Murphy
PURPOSE: The codebase uses linear magnetic suppression:
           β = exp(-B/B_crit)              [CP3 L1802, C++ COMPRESSED_RESONANCE]
         This module adds the Gaussian (squared) form which is physically more
         appropriate for smooth SCm manifold collapse:
           A_SCm(B) = exp[-(B/B_crit)²]   [Gaussian suppression]

         Both forms are computed for comparison along with:
         - Transition field B_trans where the two forms diverge by >10%
         - Critical field analysis for LENR lab systems
         - SCm collapse threshold for cosmogenesis

PHYSICS:
  The linear form exp(-B/B_crit) is a simple exponential decay.
  The Gaussian form exp[-(B/B_crit)²] has:
    - Flatter response at low B (preserving SCm coherence)
    - Sharper cutoff near B_crit (cleaner transition)
    - Matches BCS superconductor gap equation behavior
  At B = B_crit: linear → e⁻¹ ≈ 0.368, Gaussian → e⁻¹ ≈ 0.368 (same!)
  At B = 0.5·B_crit: linear → 0.607, Gaussian → 0.779 (Gaussian 28% higher)
  At B = 2·B_crit: linear → 0.135, Gaussian → 0.018 (Gaussian 87% lower)

ARCHITECTURE: Pure calculator, no hardcoded systems. Receives dataset, returns equations.
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
K_B = 1.381e-23
MU_0 = 4 * PI * 1e-7

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36

# Critical fields
B_CRIT_MAGNETAR = 4.4e13   # T — magnetar Schwinger limit
B_CRIT_LAB = 1.5            # T — typical lab NdFeB + solenoid
B_CRIT_BCS = 10.0           # T — typical Type-II SC upper critical field

# LENR
OMEGA_SCM = 2 * PI * 1.25e12  # rad/s — 1.25 THz SCm phonon


# ── §1  ACTIVATION FUNCTION MODELS ───────────────────────────────────────

class SCmActivationFunction:
    """
    SCm manifold activation/suppression as a function of magnetic field B.

    Two models:
      Linear (existing):   β(B) = exp(-B/B_crit)
      Gaussian (new):      A_SCm(B) = exp[-(B/B_crit)²]

    Also provides:
      Combined:  α(B) = w·β(B) + (1-w)·A_SCm(B)  for smooth interpolation
      Derivative: dA/dB for sensitivity analysis
      B_transition: field where |Gaussian - Linear| > threshold
    """

    def __init__(self, B_crit: float = B_CRIT_MAGNETAR):
        self.B_crit = B_crit

    # ── 1a. Linear (existing CP3 form) ───────────────────────────────────

    def beta_linear(self, B: float) -> Dict[str, Any]:
        """
        β(B) = exp(-B/B_crit)
        Existing form from CP3 L1802 and C++ COMPRESSED_RESONANCE modules.
        """
        x = B / self.B_crit
        beta = math.exp(-x)
        dbeta_dB = -beta / self.B_crit  # derivative

        return {
            "model": "linear",
            "beta": beta,
            "B": B,
            "B_crit": self.B_crit,
            "B_over_Bcrit": x,
            "d_beta_dB": dbeta_dB,
            "equation": "β(B) = exp(-B/B_crit)",
        }

    # ── 1b. Gaussian (new SCm form) ──────────────────────────────────────

    def a_scm_gaussian(self, B: float) -> Dict[str, Any]:
        """
        A_SCm(B) = exp[-(B/B_crit)²]
        Gaussian suppression — flatter at low B, sharper cutoff near B_crit.
        """
        x = B / self.B_crit
        a_scm = math.exp(-(x ** 2))
        da_dB = -2 * x / self.B_crit * a_scm  # derivative

        return {
            "model": "gaussian",
            "A_SCm": a_scm,
            "B": B,
            "B_crit": self.B_crit,
            "B_over_Bcrit": x,
            "d_A_dB": da_dB,
            "equation": "A_SCm(B) = exp[-(B/B_crit)²]",
        }

    # ── 1c. Comparison ───────────────────────────────────────────────────

    def compare(self, B: float) -> Dict[str, Any]:
        """
        Compute both forms at field B and quantify the difference.
        """
        lin = self.beta_linear(B)
        gau = self.a_scm_gaussian(B)

        beta = lin["beta"]
        a_scm = gau["A_SCm"]
        diff = abs(a_scm - beta)
        rel_diff = diff / beta if beta > 0 else float("inf")

        # Which form dominates
        if a_scm > beta:
            dominant = "gaussian_higher"
        elif a_scm < beta:
            dominant = "linear_higher"
        else:
            dominant = "equal"

        return {
            "B": B,
            "B_crit": self.B_crit,
            "B_over_Bcrit": B / self.B_crit,
            "beta_linear": beta,
            "A_SCm_gaussian": a_scm,
            "absolute_difference": diff,
            "relative_difference": rel_diff,
            "dominant": dominant,
            "linear_detail": lin,
            "gaussian_detail": gau,
        }

    # ── 1d. Combined / interpolated ──────────────────────────────────────

    def combined(self, B: float, weight_gaussian: float = 0.5) -> Dict[str, Any]:
        """
        Weighted combination: α(B) = w·A_SCm(B) + (1-w)·β(B)
        """
        lin = self.beta_linear(B)
        gau = self.a_scm_gaussian(B)
        w = weight_gaussian

        alpha = w * gau["A_SCm"] + (1 - w) * lin["beta"]

        return {
            "model": "combined",
            "alpha": alpha,
            "weight_gaussian": w,
            "beta_linear": lin["beta"],
            "A_SCm_gaussian": gau["A_SCm"],
            "B": B,
            "B_crit": self.B_crit,
            "equation": "α(B) = w·A_SCm(B) + (1-w)·β(B)",
        }

    # ── 1e. Transition field finder ──────────────────────────────────────

    def find_transition_field(self, threshold: float = 0.10,
                              n_points: int = 10000) -> Dict[str, Any]:
        """
        Find the field B where |Gaussian - Linear| / Linear > threshold.
        Scans from B=0 to B=3·B_crit.
        """
        B_max = 3 * self.B_crit
        dB = B_max / n_points

        B_trans = None
        for i in range(1, n_points + 1):
            B = i * dB
            x = B / self.B_crit
            beta = math.exp(-x)
            a_scm = math.exp(-(x ** 2))
            if beta > 0:
                rel = abs(a_scm - beta) / beta
                if rel > threshold and B_trans is None:
                    B_trans = B

        return {
            "B_transition": B_trans,
            "B_trans_over_Bcrit": B_trans / self.B_crit if B_trans else None,
            "threshold": threshold,
            "B_crit": self.B_crit,
        }

    # ── 1f. Field sweep ──────────────────────────────────────────────────

    def field_sweep(self, B_min: float = 0.0, B_max: float = None,
                    n_points: int = 100) -> Dict[str, Any]:
        """
        Sweep B from B_min to B_max, return both activation curves.
        """
        if B_max is None:
            B_max = 3 * self.B_crit

        dB = (B_max - B_min) / max(n_points - 1, 1)
        B_values = [B_min + i * dB for i in range(n_points)]
        beta_values = []
        a_scm_values = []
        diff_values = []

        for B in B_values:
            x = B / self.B_crit
            beta = math.exp(-x)
            a_scm = math.exp(-(x ** 2))
            beta_values.append(beta)
            a_scm_values.append(a_scm)
            diff_values.append(a_scm - beta)

        return {
            "n_points": n_points,
            "B_crit": self.B_crit,
            "B": B_values,
            "B_over_Bcrit": [b / self.B_crit for b in B_values],
            "beta_linear": beta_values,
            "A_SCm_gaussian": a_scm_values,
            "difference": diff_values,
        }

    # ── 1g. LENR lab system evaluation ───────────────────────────────────

    def evaluate_lenr_systems(self) -> List[Dict[str, Any]]:
        """
        Evaluate activation function for the 6 canonical LENR lab systems.
        Uses their actual B fields against magnetar B_crit.
        """
        systems = [
            ("Micro-Plasmoid Reversal", 0.5, "PAPER_859"),
            ("Water Reactor H2O2", 0.01, "PAPER_863"),
            ("Colman-Gillespie Pulsed Motor", 1.2, "PAPER_835"),
            ("Kozima Neutron Drop", 0.0, "PAPER_840"),
            ("LRC Spark-Gap Monopole", 0.1, "PAPER_864"),
            ("Caduceus Twin-Helix Motor", 0.8, "PAPER_866"),
        ]

        results = []
        for name, B, paper in systems:
            comp = self.compare(B)
            results.append({
                "system": name,
                "paper": paper,
                "B_T": B,
                "beta_linear": comp["beta_linear"],
                "A_SCm_gaussian": comp["A_SCm_gaussian"],
                "relative_difference": comp["relative_difference"],
                "note": "Lab fields B << B_crit: both ≈ 1.0 (SCm fully active)",
            })

        return results


# ── §2  MAIN ──────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("SCm Activation Function — Linear vs Gaussian B-Field Suppression")
    print("=" * 72)

    calc = SCmActivationFunction(B_crit=B_CRIT_MAGNETAR)

    # 1. Key comparison points
    test_fields = [
        ("Zero field", 0.0),
        ("Lab (1.5 T)", 1.5),
        ("SC magnet (10 T)", 10.0),
        ("Half B_crit", 0.5 * B_CRIT_MAGNETAR),
        ("At B_crit", B_CRIT_MAGNETAR),
        ("2x B_crit", 2.0 * B_CRIT_MAGNETAR),
    ]

    print(f"\n  B_crit = {B_CRIT_MAGNETAR:.2e} T (magnetar Schwinger limit)")
    print(f"\n{'Label':25s} {'B/B_crit':>10s} {'β(linear)':>12s} {'A_SCm(gauss)':>14s} {'Rel Diff':>10s}")
    print("-" * 72)

    for label, B in test_fields:
        comp = calc.compare(B)
        print(f"  {label:23s} {comp['B_over_Bcrit']:10.4f} "
              f"{comp['beta_linear']:12.6f} {comp['A_SCm_gaussian']:14.6f} "
              f"{comp['relative_difference']:10.4f}")

    # 2. Transition field
    trans = calc.find_transition_field(threshold=0.10)
    print(f"\n── Transition Field (>10% divergence) ──")
    if trans["B_transition"]:
        print(f"  B_trans = {trans['B_transition']:.4e} T  ({trans['B_trans_over_Bcrit']:.4f} × B_crit)")
    else:
        print(f"  Not reached in scan range")

    # 3. LENR lab systems
    print(f"\n── LENR Lab Systems (B << B_crit) ──")
    lenr = calc.evaluate_lenr_systems()
    for sys in lenr:
        print(f"  {sys['system']:35s}  B={sys['B_T']:.2f} T  "
              f"β={sys['beta_linear']:.10f}  A_SCm={sys['A_SCm_gaussian']:.10f}  "
              f"diff={sys['relative_difference']:.2e}")

    # 4. Physics insight
    print(f"\n── Physics Interpretation ──")
    print(f"  At lab fields (0–10 T):  Both ≈ 1.0000 (SCm fully active)")
    print(f"  At 0.5·B_crit:  Linear = 0.607, Gaussian = 0.779 (28% difference)")
    print(f"  At B_crit:      Both = e⁻¹ ≈ 0.368 (crossing point)")
    print(f"  At 2·B_crit:    Linear = 0.135, Gaussian = 0.018 (87% difference)")
    print(f"  → Gaussian preserves SCm at low B, collapses faster at high B")
    print(f"  → Matches BCS superconductor gap behavior near T_c")

    # Export
    export = {
        "calculator": "SCmActivationFunction",
        "B_crit": B_CRIT_MAGNETAR,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "comparisons": [calc.compare(B) for _, B in test_fields],
        "transition": trans,
        "lenr_systems": lenr,
    }
    out_path = "scm_activation_results.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(export, f, indent=2, default=str)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print("COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
