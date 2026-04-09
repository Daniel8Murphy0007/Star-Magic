"""
Mock Theta Functions f₂₆(q), φ₂₆(q), ψ₂₆(q) — Proper q-Series with (q;q)ₙ

Session 204 | Daniel Murphy
PURPOSE: The codebase has ONE mock-theta-related calculator:

         1. CondensedPhysics3.py L2734 (RamanujanPolynomialsQ26Calculator):
              Q_26 = Σ_{n=1}^{26} c_n exp(-nπ√n)
            This is a SIMPLIFIED exponential sum labelled "mock theta-like",
            NOT a proper Ramanujan mock theta function. It lacks the
            essential q-Pochhammer symbol (q;q)_n = Π_{k=0}^{n-1}(1-q^{k+1}).

         2. qcalcgeom_helpers.py L506:
              def pochhammer(k, n=26) → rising factorial (x)_n = x(x+1)...(x+n-1)
            This is the RISING FACTORIAL Pochhammer, NOT the q-Pochhammer.

         This module implements:
           1. q-Pochhammer symbol (q;q)_n = Π_{k=0}^{n-1}(1 - q^{k+1})
           2. Third-order mock theta f₂₆(q) = Σ_{n=0}^{25} q^{n²}/(-q;q)_n²
           3. Generalized φ₂₆(q) = Σ_{n=0}^{25} q^{n²}/(−q²;q²)_n
           4. Generalized ψ₂₆(q) = Σ_{n=1}^{26} q^{n²}/(q;q²)_n
           5. UQFF coupling via q = [SSq]·exp(−κt) thermal parameter
           6. VDS-weighted 26-state mock theta amplitudes

         Mathematical framework:
           The Ramanujan mock theta functions are q-hypergeometric series
           that transform like modular forms near rational cusps but lack
           a modular complement. Zwegers (2002) showed they are holomorphic
           parts of harmonic Maass forms. The 26-state UQFF coupling maps
           each layer i → q_i = [SSq]·exp(−κ·t_i) where t_i = i·T/26.

ARCHITECTURE: Pure calculator. No hardcoded astrophysical systems.
"""

import json
import math
from dataclasses import dataclass, field
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
M_SUN = 1.989e30

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36

# 26-state parameters
N_LEVELS = 26


# ── §1  q-POCHHAMMER SYMBOL ──────────────────────────────────────────────

def q_pochhammer(a: float, q: float, n: int) -> float:
    """
    q-Pochhammer symbol (a; q)_n = Π_{k=0}^{n-1} (1 - a·q^k).

    For n=0, returns 1 (empty product).
    For n<0, undefined — raises ValueError.

    This is the INFINITE PRODUCT building block for q-series,
    distinct from the rising factorial Pochhammer (x)_n.

    Standard special cases:
      (q; q)_n = Π_{k=0}^{n-1}(1 - q^{k+1})  [set a=q]
      (q; q)_∞ = Euler function               [n → ∞]
    """
    if n < 0:
        raise ValueError(f"q-Pochhammer requires n≥0, got n={n}")
    if n == 0:
        return 1.0
    product = 1.0
    for k in range(n):
        factor = 1.0 - a * q**k
        product *= factor
        if abs(product) < 1e-300:
            return 0.0
    return product


def q_pochhammer_neg(q: float, n: int) -> float:
    """(-q; q)_n = Π_{k=0}^{n-1}(1 + q^{k+1}), appears in f(q)."""
    return q_pochhammer(-q, q, n)


# ── §2  MOCK THETA FUNCTIONS ─────────────────────────────────────────────

class MockThetaQ26:
    """
    Ramanujan third-order mock theta functions truncated at 26 states.

    f(q) = Σ_{n=0}^{∞} q^{n²} / (-q;q)_n²      [third-order, type f]
    φ(q) = Σ_{n=0}^{∞} q^{n²} / (-q²;q²)_n      [third-order, type φ]
    ψ(q) = Σ_{n=1}^{∞} q^{n²} / (q;q²)_n        [third-order, type ψ]

    The 26-state UQFF truncation:
      f₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q;q)_n²
      φ₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q²;q²)_n
      ψ₂₆(q) = Σ_{n=1}^{26} q^{n²} / (q;q²)_n

    UQFF coupling: q = [SSq]·exp(−κ·t) maps the unit disk parameter
    to the superconducting manifold decay.

    Reference: Ramanujan's "Lost Notebook" (1920), rediscovered by
    Andrews (1976); Zwegers thesis (2002) on harmonic Maass forms.
    """

    @staticmethod
    def f26(q: float) -> Dict[str, Any]:
        """
        f₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q;q)_n².

        Third-order mock theta function of the first kind.
        """
        if abs(q) >= 1.0:
            return {"error": "|q| must be < 1", "q": q}

        terms = []
        total = 0.0
        for n in range(N_LEVELS):
            numerator = q ** (n * n)
            denominator = q_pochhammer_neg(q, n) ** 2
            if abs(denominator) < 1e-300:
                break
            term = numerator / denominator
            total += term
            terms.append({
                "n": n,
                "q^n²": numerator,
                "(-q;q)_n²": denominator,
                "term": term,
                "partial_sum": total,
            })

        return {
            "f26_q": total,
            "q": q,
            "n_terms": len(terms),
            "terms": terms[:6],  # first 6 for display
            "equation": "f₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q;q)_n²",
            "convergence": "mock modular — transforms under τ→-1/τ up to shadow",
        }

    @staticmethod
    def phi26(q: float) -> Dict[str, Any]:
        """
        φ₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q²;q²)_n.

        Third-order mock theta function of the second kind.
        """
        if abs(q) >= 1.0:
            return {"error": "|q| must be < 1", "q": q}

        terms = []
        total = 0.0
        q2 = q * q
        for n in range(N_LEVELS):
            numerator = q ** (n * n)
            denominator = q_pochhammer(-q2, q2, n)
            if abs(denominator) < 1e-300:
                break
            term = numerator / denominator
            total += term
            terms.append({
                "n": n,
                "q^n²": numerator,
                "(-q²;q²)_n": denominator,
                "term": term,
                "partial_sum": total,
            })

        return {
            "phi26_q": total,
            "q": q,
            "n_terms": len(terms),
            "terms": terms[:6],
            "equation": "φ₂₆(q) = Σ_{n=0}^{25} q^{n²} / (-q²;q²)_n",
        }

    @staticmethod
    def psi26(q: float) -> Dict[str, Any]:
        """
        ψ₂₆(q) = Σ_{n=1}^{26} q^{n²} / (q;q²)_n.

        Third-order mock theta function of the third kind.
        Note: sum starts at n=1.
        """
        if abs(q) >= 1.0:
            return {"error": "|q| must be < 1", "q": q}

        terms = []
        total = 0.0
        q2 = q * q
        for n in range(1, N_LEVELS + 1):
            numerator = q ** (n * n)
            denominator = q_pochhammer(q, q2, n)
            if abs(denominator) < 1e-300:
                break
            term = numerator / denominator
            total += term
            terms.append({
                "n": n,
                "q^n²": numerator,
                "(q;q²)_n": denominator,
                "term": term,
                "partial_sum": total,
            })

        return {
            "psi26_q": total,
            "q": q,
            "n_terms": len(terms),
            "terms": terms[:6],
            "equation": "ψ₂₆(q) = Σ_{n=1}^{26} q^{n²} / (q;q²)_n",
        }


# ── §3  UQFF COUPLING ────────────────────────────────────────────────────

class MockThetaUQFF:
    """
    UQFF 26-state mock theta coupling.

    The thermal q parameter is coupled to the UQFF superconducting
    manifold decay: q_i = [SSq] · exp(−κ · t_i) for each layer i.

    The VDS (Vacuum Density Series) amplitude at level i is:
      A_i = (2π)^{i/6} · ρ_SCm/ρ_UA

    The coupled mock theta amplitude:
      Θ₂₆ = Σ_{i=1}^{26} A_i · [f₂₆(q_i) + φ₂₆(q_i) + ψ₂₆(q_i)]

    where the three mock theta species encode different sectors:
      f₂₆ → magnetic channel (Ug1 sector)
      φ₂₆ → reactive channel (Ug2 sector)
      ψ₂₆ → rotational channel (Ug3 sector)
    """

    @staticmethod
    def compute_q(t: float, ssq: float = SSQ, kappa: float = KAPPA) -> float:
        """q(t) = [SSq]·exp(−κ·t), the thermal q-series parameter."""
        return ssq * math.exp(-kappa * t)

    @staticmethod
    def vds_amplitude(level: int) -> float:
        """A_i = (2π)^{i/6} · ρ_SCm/ρ_UA."""
        return (2 * PI) ** (level / 6.0) * (RHO_SCM / RHO_UA)

    @staticmethod
    def coupled_theta26(dataset: dict) -> Dict[str, Any]:
        """
        Full 26-state coupled mock theta amplitude.

        Parameters (from dataset dict):
          t    : observation time (s), default 0
          ssq  : superconducting quantum strength, default 0.57
          kappa: decay constant (1/s), default 5.787e-9
        """
        t = dataset.get("t", 0.0)
        ssq = dataset.get("ssq", SSQ)
        kappa = dataset.get("kappa", KAPPA)

        theta_total = 0.0
        layer_results = []

        for i in range(1, N_LEVELS + 1):
            t_i = i * t / N_LEVELS if t > 0 else 0.0
            q_i = MockThetaUQFF.compute_q(t_i, ssq, kappa)
            A_i = MockThetaUQFF.vds_amplitude(i)

            f_val = MockThetaQ26.f26(q_i)
            phi_val = MockThetaQ26.phi26(q_i)
            psi_val = MockThetaQ26.psi26(q_i)

            f_num = f_val.get("f26_q", 0.0)
            phi_num = phi_val.get("phi26_q", 0.0)
            psi_num = psi_val.get("psi26_q", 0.0)

            layer_theta = A_i * (f_num + phi_num + psi_num)
            theta_total += layer_theta

            layer_results.append({
                "level": i,
                "q_i": q_i,
                "A_i": A_i,
                "f26": f_num,
                "phi26": phi_num,
                "psi26": psi_num,
                "layer_theta": layer_theta,
            })

        return {
            "Theta_26": theta_total,
            "n_layers": N_LEVELS,
            "q_0": ssq,
            "t": t,
            "layers": layer_results[:5],  # first 5 for display
            "equation": (
                "Θ₂₆ = Σ_{i=1}^{26} A_i·[f₂₆(q_i) + φ₂₆(q_i) + ψ₂₆(q_i)] "
                "where q_i = [SSq]·exp(−κt_i), A_i = (2π)^{i/6}·ρ_SCm/ρ_UA"
            ),
        }


# ── §4  SELF-TEST ─────────────────────────────────────────────────────────

def self_test() -> Dict[str, Any]:
    """Validate mock theta functions at q = SSQ = 0.57."""
    results = {}
    q = SSQ

    # Test q-Pochhammer basic properties
    # (q;q)_0 = 1 (empty product)
    assert q_pochhammer(q, q, 0) == 1.0, "(q;q)_0 must be 1"
    # (q;q)_1 = 1 - q
    val_1 = q_pochhammer(q, q, 1)
    assert abs(val_1 - (1 - q)) < 1e-14, f"(q;q)_1={val_1} ≠ 1-q={1-q}"
    # (q;q)_2 = (1-q)(1-q²)
    val_2 = q_pochhammer(q, q, 2)
    expected_2 = (1 - q) * (1 - q**2)
    assert abs(val_2 - expected_2) < 1e-14, f"(q;q)_2={val_2} ≠ {expected_2}"
    results["q_pochhammer"] = "PASS"

    # Test f₂₆(q) at q=0.57
    f_result = MockThetaQ26.f26(q)
    f_val = f_result["f26_q"]
    # f(q) must be > 1 (first term is q^0/1 = 1, remaining positive)
    assert f_val > 1.0, f"f₂₆(0.57) = {f_val} should be > 1"
    assert f_val < 1e10, f"f₂₆(0.57) = {f_val} diverging"
    results["f26"] = {"value": f_val, "status": "PASS"}

    # Test φ₂₆(q)
    phi_result = MockThetaQ26.phi26(q)
    phi_val = phi_result["phi26_q"]
    assert phi_val > 1.0, f"φ₂₆(0.57) = {phi_val} should be > 1"
    assert phi_val < 1e10, f"φ₂₆(0.57) = {phi_val} diverging"
    results["phi26"] = {"value": phi_val, "status": "PASS"}

    # Test ψ₂₆(q)
    psi_result = MockThetaQ26.psi26(q)
    psi_val = psi_result["psi26_q"]
    assert psi_val > 0.0, f"ψ₂₆(0.57) = {psi_val} should be > 0"
    assert psi_val < 1e10, f"ψ₂₆(0.57) = {psi_val} diverging"
    results["psi26"] = {"value": psi_val, "status": "PASS"}

    # Cross-check: f(0) = 1, φ(0) = 1, ψ(0) = 0 (sum from n=1)
    f0 = MockThetaQ26.f26(0.0)["f26_q"]
    phi0 = MockThetaQ26.phi26(0.0)["phi26_q"]
    psi0 = MockThetaQ26.psi26(0.0)["psi26_q"]
    assert abs(f0 - 1.0) < 1e-14, f"f₂₆(0)={f0} ≠ 1"
    assert abs(phi0 - 1.0) < 1e-14, f"φ₂₆(0)={phi0} ≠ 1"
    assert abs(psi0) < 1e-14, f"ψ₂₆(0)={psi0} ≠ 0"
    results["boundary_q0"] = "PASS"

    # Test UQFF coupled theta (t=0 → q=SSQ at all layers)
    coupled = MockThetaUQFF.coupled_theta26({"t": 0.0})
    theta = coupled["Theta_26"]
    assert theta != 0.0, "Θ₂₆ should be non-zero at t=0"
    results["coupled_theta26"] = {"Theta_26": theta, "status": "PASS"}

    # All three mock theta values should be positive and finite
    assert f_val > 0 and phi_val > 0 and psi_val > 0, \
        f"All mock thetas should be positive: f={f_val}, φ={phi_val}, ψ={psi_val}"
    results["all_positive"] = "PASS"

    results["overall"] = "PASS"
    return results


if __name__ == "__main__":
    print("=" * 72)
    print("Mock Theta Q₂₆ — Self-Test")
    print("=" * 72)
    r = self_test()
    print(json.dumps(r, indent=2, default=str))
    print("=" * 72)
    print(f"Overall: {r['overall']}")
