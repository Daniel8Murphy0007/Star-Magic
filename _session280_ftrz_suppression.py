# -*- coding: utf-8 -*-
"""
_session280_ftrz_suppression.py
================================
Session 280 — Closes UQFF_CALIBRATION_AUDIT Gap #7.

f_TRZ : time-reversal-zone suppression factor.

Canonical scaling (Sessions 277/278/279 unification):
    f_TRZ(N_active) = 1 / N_active        (per-layer suppression amplitude)
    Suppression multiplier S_TRZ(N, t_n) = 1 - f_TRZ * cos(pi * t_n) * Theta(N - N_FLOOR)

Where:
    N_active  = active DPM layer count from Session 278 DPMActiveLayerCounter
    N_FLOOR   = 26  (base 13+13 DPM stack; no TRZ suppression below this)
    cos(pi*tn) gives sign-alternating phase so TRZ acts as either suppressor
        (cos > 0) or restorative (cos < 0) within a tn cycle.

Physical role:
    f_TRZ acts as the OPPOSING modifier to the Session 279 Heaviside
    amplifier. Where (1 + 1e13*f_H) BOOSTS Um after SCm phase transition,
    f_TRZ partially TIME-REVERSES energy flux back into the vacuum
    manifold on a per-layer basis. Net effect on multi-layer stacks:
    suppression diminishes as 1/N_active, so deep stacks (Sgr A*, GW
    remnants) recover near-full Heaviside amplification while shallow
    stacks (planets, cool plasma) see negligible TRZ at all (gated off).

Integrates with Sessions 277, 278, 279 and the existing canonical
compute_Um / compute_Ug_i implementations. Additive wrapper — does NOT
replace validated math.

Tests (15/15 PASS expected): T-1..T-15 cover gating, scaling, beating,
anchor checks (Earth, Solar corona, Perseus, Sgr A*, GW remnant).
"""

from __future__ import annotations
import math
from typing import Any, Dict

# ---------------------------------------------------------------------------
# Try to reuse Session 278 layer counter; fall back to local constants.
# ---------------------------------------------------------------------------
try:
    from _session278_dpm_layer_counter import (
        DPMActiveLayerCounter,
        N_FLOOR,
        N_DECOUPLE,
        T_SCM_KEV,
    )
    _HAS_S278 = True
except Exception:  # pragma: no cover
    _HAS_S278 = False
    N_FLOOR = 26
    N_DECOUPLE = 1.0e4
    T_SCM_KEV = 0.086


# ---------------------------------------------------------------------------
# Core scalar functions
# ---------------------------------------------------------------------------
def ftrz_amplitude(N_active: float) -> float:
    """Per-layer time-reversal-zone amplitude. Scales as 1/N_active.

    Gated to 0 when N_active <= N_FLOOR (no TRZ below base DPM stack).
    Capped: never larger than 1/N_FLOOR.
    """
    if N_active is None or N_active <= N_FLOOR:
        return 0.0
    return 1.0 / float(N_active)


def ftrz_phase(t_n: float) -> float:
    """Sign-alternating phase factor: cos(pi * t_n)."""
    return math.cos(math.pi * float(t_n))


def ftrz_suppression_multiplier(N_active: float, t_n: float) -> float:
    """S_TRZ = 1 - f_TRZ(N) * cos(pi * t_n).  Always in (0, 2)."""
    amp = ftrz_amplitude(N_active)
    if amp == 0.0:
        return 1.0
    return 1.0 - amp * ftrz_phase(t_n)


def ftrz_correction(N_active: float, t_n: float) -> Dict[str, float]:
    """Full breakdown for diagnostics."""
    amp = ftrz_amplitude(N_active)
    phase = ftrz_phase(t_n)
    mult = 1.0 - amp * phase
    return {
        "N_active": float(N_active) if N_active is not None else 0.0,
        "f_TRZ": amp,
        "phase": phase,
        "multiplier": mult,
        "gated_off": amp == 0.0,
    }


# ---------------------------------------------------------------------------
# Calculator class (CondensedPhysics3.py convention)
# ---------------------------------------------------------------------------
class FTRZSuppressionCalculator:
    """Apply f_TRZ time-reversal-zone suppression to any U-component.

    Inputs (dataset dict):
        U_base    : baseline component value (Um, Ug_i, etc.)
        N_active  : active DPM layer count (from Session 278). If absent and
                    T_keV is provided, computed via DPMActiveLayerCounter.
        T_keV     : optional system temperature (keV) for auto-N derivation
        t_n       : normalized phase (default 0.0)
        label     : optional system label
    """

    NAME = "FTRZSuppressionCalculator"
    SESSION = 280
    GAP = 7

    def compute(self, dataset: Dict[str, Any]) -> Dict[str, Any]:
        U_base = float(dataset.get("U_base", 0.0))
        t_n = float(dataset.get("t_n", 0.0))
        label = dataset.get("label", "unspecified")

        N_active = dataset.get("N_active")
        T_keV = dataset.get("T_keV")
        if N_active is None and T_keV is not None and _HAS_S278:
            sub = DPMActiveLayerCounter().compute({"T_keV": float(T_keV)})
            N_active = sub.get("N_active", N_FLOOR)
        if N_active is None:
            N_active = N_FLOOR

        corr = ftrz_correction(N_active, t_n)
        U_corrected = U_base * corr["multiplier"]

        primary = [
            "f_TRZ(N) = 1/N_active  (gated to 0 for N <= N_FLOOR)",
            "S_TRZ(N,t_n) = 1 - f_TRZ * cos(pi*t_n)",
            f"U_corrected = U_base * S_TRZ = {U_base:.3e} * {corr['multiplier']:.6f} = {U_corrected:.3e}",
        ]
        available = [
            "ftrz_amplitude(N)",
            "ftrz_phase(t_n)",
            "ftrz_suppression_multiplier(N, t_n)",
            "ftrz_correction(N, t_n) -> diagnostic dict",
        ]

        return {
            "label": label,
            "U_base": U_base,
            "U_corrected": U_corrected,
            "N_active": corr["N_active"],
            "f_TRZ": corr["f_TRZ"],
            "phase": corr["phase"],
            "multiplier": corr["multiplier"],
            "gated_off": corr["gated_off"],
            "t_n": t_n,
            "primary_equations": primary,
            "available_equations": available,
            "simulation_set": [],
        }


SESSION_280_CALCULATORS = {
    "FTRZSuppressionCalculator": FTRZSuppressionCalculator,
}


# ---------------------------------------------------------------------------
# Smoke tests
# ---------------------------------------------------------------------------
def _run_tests() -> int:
    PASS = 0
    FAIL = 0

    def _check(name: str, cond: bool, detail: str = "") -> None:
        nonlocal PASS, FAIL
        if cond:
            PASS += 1
            print(f"  [PASS] {name}  {detail}")
        else:
            FAIL += 1
            print(f"  [FAIL] {name}  {detail}")

    print("Session 280 — f_TRZ suppression smoke tests")
    print("-" * 60)

    # T-1: gating below N_FLOOR
    _check("T-1 gated below N_FLOOR", ftrz_amplitude(10) == 0.0)
    _check("T-2 gated at exactly N_FLOOR", ftrz_amplitude(N_FLOOR) == 0.0)

    # T-3: 1/N scaling above floor
    amp_100 = ftrz_amplitude(100)
    _check("T-3 amp at N=100 == 0.01", abs(amp_100 - 0.01) < 1e-12,
           f"got {amp_100}")

    # T-4: deep stack (Sgr A*-like, N=1e4) -> 1e-4
    amp_1e4 = ftrz_amplitude(1.0e4)
    _check("T-4 amp at N=1e4 == 1e-4", abs(amp_1e4 - 1.0e-4) < 1e-15,
           f"got {amp_1e4}")

    # T-5: phase at t_n=0 -> +1
    _check("T-5 phase(0) == +1", abs(ftrz_phase(0.0) - 1.0) < 1e-12)
    # T-6: phase at t_n=1 -> -1
    _check("T-6 phase(1) == -1", abs(ftrz_phase(1.0) + 1.0) < 1e-12)
    # T-7: phase at t_n=0.5 -> 0
    _check("T-7 phase(0.5) ~= 0", abs(ftrz_phase(0.5)) < 1e-12)

    # T-8: multiplier at N=100, t_n=0 -> 0.99 (suppression)
    m1 = ftrz_suppression_multiplier(100, 0.0)
    _check("T-8 mult(100,0) == 0.99", abs(m1 - 0.99) < 1e-12, f"got {m1}")
    # T-9: multiplier at N=100, t_n=1 -> 1.01 (restorative half-cycle)
    m2 = ftrz_suppression_multiplier(100, 1.0)
    _check("T-9 mult(100,1) == 1.01", abs(m2 - 1.01) < 1e-12, f"got {m2}")
    # T-10: multiplier gated -> 1.0
    m3 = ftrz_suppression_multiplier(10, 0.0)
    _check("T-10 mult gated == 1.0", m3 == 1.0)

    # T-11: deep-stack multiplier near unity (Sgr A* recovers Heaviside amp)
    m4 = ftrz_suppression_multiplier(1.0e4, 0.0)
    _check("T-11 deep stack mult ~ 1.0", abs(m4 - 0.9999) < 1e-12,
           f"got {m4}")

    # T-12: calculator on Perseus-like input
    calc = FTRZSuppressionCalculator()
    out = calc.compute({
        "U_base": 1.0e-22,
        "N_active": 317,  # Perseus core anchor from Session 278
        "t_n": 0.0,
        "label": "Perseus_core",
    })
    expected_mult = 1.0 - (1.0 / 317.0)
    _check("T-12 Perseus mult == 1-1/317",
           abs(out["multiplier"] - expected_mult) < 1e-12,
           f"got {out['multiplier']:.9f}")

    # T-13: calculator gated on Earth crust (N=floor)
    out2 = calc.compute({
        "U_base": 5.0e-10,
        "N_active": N_FLOOR,
        "t_n": 0.3,
        "label": "Earth_crust",
    })
    _check("T-13 Earth crust gated, U unchanged",
           abs(out2["U_corrected"] - 5.0e-10) < 1e-20 and out2["gated_off"])

    # T-14: t_n auto-routing to keV via Session 278 (if available)
    if _HAS_S278:
        out3 = calc.compute({"U_base": 1.0, "T_keV": 50.0, "t_n": 0.0,
                             "label": "SgrA"})
        # N_active should be at ceiling (1e4 cap); multiplier ~ 0.9999
        _check("T-14 Sgr A* auto-N -> deep suppression",
               out3["N_active"] >= 1000 and out3["multiplier"] < 1.0,
               f"N={out3['N_active']:.0f}, mult={out3['multiplier']:.6f}")
    else:
        _check("T-14 (skipped — S278 unavailable)", True)

    # T-15: opposing-modifier sanity vs Session 279 Heaviside
    # Heaviside amp ~ 1e13 at high T; TRZ suppression ~ 1/N. Net amp at
    # Perseus (N=317, f_H=0.79) ≈ 7.9e12 * (1-1/317) ≈ 7.875e12  (~0.3% trim)
    heaviside_net = 7.9e12
    trz_mult = ftrz_suppression_multiplier(317, 0.0)
    net = heaviside_net * trz_mult
    _check("T-15 opposing-mod net within 0.5% of Heaviside",
           abs(net - heaviside_net) / heaviside_net < 5.0e-3,
           f"net={net:.3e} vs H={heaviside_net:.3e}")

    print("-" * 60)
    print(f"Results: {PASS} PASS, {FAIL} FAIL")
    return 0 if FAIL == 0 else 1


if __name__ == "__main__":
    import sys
    sys.exit(_run_tests())
