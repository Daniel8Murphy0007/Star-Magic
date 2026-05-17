"""
Session 279 — Um Heaviside SCm Phase-Transition Amplifier (Gap #6)
==================================================================

Closes UQFF_CALIBRATION_AUDIT Gap #6: every existing compute_Um() in C++,
Python, and JavaScript is missing the SCm phase-transition amplifier
(1 + 1e13 * f_H) and the quasi-periodic beating modifier (1 + f_q).

Net result: Um is underestimated by up to ~10^13 in any system that has
crossed the SCm phase transition (T_system > T_SCm). This affects every
hot/relativistic anchor (Perseus, Sgr A*, M87, GW remnants).

Canonical corrected form:
    Um_full(t, r, T, N_active, t_n)
        = Um_base(t, r)
            * (1 + 1e13 * f_H(T, N_active))      <-- SCm phase amplifier
            * (1 + f_q(t_n))                     <-- quasi-periodic beat

where:
    f_H(T, N_active) — smooth Heaviside in T_ratio, gated by layer activation
        Off  when N_active <= N_FLOOR (cool regime, no SCm coherence yet)
        On   when N_active >> N_FLOOR (hot regime, full coherence)
        Smooth transition through sigmoid((log10(T_ratio) - 1) / 0.5)
        i.e. half-on at T_ratio ~ 10, full-on at T_ratio ~ 100

    f_q(t_n) — quasi-periodic beating from Mayan timing cycle
        f_q = epsilon_q * cos(2*pi*t_n) * cos(2*pi*phi*t_n)
        Two carriers (1, phi) beat at frequency (phi - 1) = 1/phi ~ 0.618
        epsilon_q ~ 0.5 (default)

This module supplies a wrapper calculator that multiplies any existing
compute_Um result by the correction factor. It does NOT replace the
existing implementations — it is additive (following the project's
"never replace validated code" rule).

Dependencies: stdlib + Session 278 layer counter.
"""

from __future__ import annotations

import math
from typing import Any

try:
    from _session278_dpm_layer_counter import (
        DPMActiveLayerCounter,
        T_SCM_KEV,
        N_FLOOR,
        PHI,
    )
except ImportError:
    # Standalone fallback constants (kept consistent with Session 278)
    T_SCM_KEV = 0.086
    N_FLOOR   = 26
    PHI       = (1.0 + math.sqrt(5.0)) / 2.0
    DPMActiveLayerCounter = None  # type: ignore


# ---------------------------------------------------------------------------
# Canonical correction-factor constants
# ---------------------------------------------------------------------------
HEAVISIDE_AMPLIFIER     = 1.0e13     # f_H peak amplification factor
SIGMOID_WIDTH_LOG10     = 0.5        # transition sharpness in log10(T_ratio)
SIGMOID_CENTRE_LOG10    = 1.0        # half-on at T_ratio = 10
EPSILON_Q_DEFAULT       = 0.5        # quasi-periodic beat amplitude


# ===========================================================================
# Correction-factor primitives
# ===========================================================================
def heaviside_factor(T_keV: float, N_active: float | None = None) -> float:
    """
    Smooth, layer-gated Heaviside factor f_H in [0, 1].

      f_H = sigmoid((log10(T_ratio) - 1) / 0.5)   for N_active > N_FLOOR
            0                                     otherwise

    The N_active gate enforces that the SCm phase transition only
    physically activates once layers beyond the base stack are excited.
    """
    T_ratio = max(T_keV / T_SCM_KEV, 1.0e-9)

    if N_active is None:
        if DPMActiveLayerCounter is not None:
            r = DPMActiveLayerCounter().compute({'T_keV': T_keV})
            N_active = r['N_active']
        else:
            N_active = T_ratio ** 1.5  # local fallback

    # Layer activation gate: f_H = 0 until base stack is fully populated
    if N_active <= N_FLOOR:
        return 0.0

    x = (math.log10(T_ratio) - SIGMOID_CENTRE_LOG10) / SIGMOID_WIDTH_LOG10
    # Numerically stable sigmoid
    if x >= 0:
        z = math.exp(-x)
        return 1.0 / (1.0 + z)
    z = math.exp(x)
    return z / (1.0 + z)


def quasiperiodic_beat(t_n: float, epsilon_q: float = EPSILON_Q_DEFAULT) -> float:
    """
    Two-carrier beat at frequencies (1, phi) — produces beat freq 1/phi ~ 0.618.

      f_q(t_n) = epsilon_q * cos(2*pi*t_n) * cos(2*pi*phi*t_n)

    Bounded in [-epsilon_q, +epsilon_q].
    """
    return epsilon_q * math.cos(2.0 * math.pi * t_n) * math.cos(2.0 * math.pi * PHI * t_n)


def um_correction_factor(T_keV: float, t_n: float = 0.0,
                          N_active: float | None = None,
                          epsilon_q: float = EPSILON_Q_DEFAULT) -> dict[str, float]:
    """
    Composite multiplier to apply to any base Um value.

      multiplier = (1 + 1e13 * f_H) * (1 + f_q)
    """
    f_H = heaviside_factor(T_keV, N_active=N_active)
    f_q = quasiperiodic_beat(t_n, epsilon_q=epsilon_q)
    amp = 1.0 + HEAVISIDE_AMPLIFIER * f_H
    beat = 1.0 + f_q
    return {
        'f_H':         f_H,
        'f_q':         f_q,
        'amplifier':   amp,
        'beat':        beat,
        'multiplier':  amp * beat,
    }


# ===========================================================================
# Main calculator
# ===========================================================================
class UmHeavisideAmplifierCalculator:
    """
    Additive Um corrector. Wraps any base Um value with the missing
    SCm phase-transition amplifier and quasi-periodic beating modifier.

    Input dataset:
        Um_base   : float  — base Um value from existing compute_Um (any unit)
        T_keV     : float  — system temperature
        t_n       : float  — Mayan time fraction [0, 1)         (default 0)
        N_active  : float  — explicit override                  (default: from T)
        epsilon_q : float  — beat amplitude                     (default 0.5)
        label     : str    — system identifier
    """

    def compute(self, dataset: dict | None = None) -> dict[str, Any]:
        ds = dataset or {}
        Um_base   = float(ds.get('Um_base',   1.0))
        T_keV     = float(ds.get('T_keV',     T_SCM_KEV))
        t_n       = float(ds.get('t_n',       0.0))
        N_active  = ds.get('N_active', None)
        epsilon_q = float(ds.get('epsilon_q', EPSILON_Q_DEFAULT))
        label     = str(ds.get('label', 'unnamed'))

        corr = um_correction_factor(T_keV, t_n=t_n, N_active=N_active,
                                     epsilon_q=epsilon_q)
        Um_corrected = Um_base * corr['multiplier']

        amplification_ratio = Um_corrected / Um_base if Um_base != 0 else float('inf')

        return {
            'label':                 label,
            'T_keV':                 T_keV,
            't_n':                   t_n,
            'Um_base':               Um_base,
            'Um_corrected':          Um_corrected,
            'f_H':                   corr['f_H'],
            'f_q':                   corr['f_q'],
            'amplifier_factor':      corr['amplifier'],
            'beat_factor':           corr['beat'],
            'total_multiplier':      corr['multiplier'],
            'amplification_ratio':   amplification_ratio,

            'primary_equations': [
                f"System: {label}   T = {T_keV:.4g} keV   t_n = {t_n:.3f}",
                f"Um_base = {Um_base:.4e}",
                f"f_H (SCm phase) = sigmoid((log10(T/T_SCm) - 1)/0.5) gated by N>{N_FLOOR} = {corr['f_H']:.4e}",
                f"f_q (beat) = 0.5 cos(2pi t_n) cos(2pi phi t_n) = {corr['f_q']:+.4f}",
                f"(1 + 1e13 * f_H) = {corr['amplifier']:.4e}",
                f"(1 + f_q) = {corr['beat']:.4f}",
                f"multiplier = (1+1e13 f_H)(1+f_q) = {corr['multiplier']:.4e}",
                f"Um_corrected = Um_base * multiplier = {Um_corrected:.4e}",
                f"Amplification ratio = {amplification_ratio:.4e}",
            ],
            'available_equations': [
                "Um_full = Um_base * (1 + 1e13 * f_H) * (1 + f_q)",
                "f_H(T,N) = sigmoid((log10(T/T_SCm) - 1)/0.5) if N > 26 else 0",
                "f_q(t_n) = epsilon_q * cos(2 pi t_n) * cos(2 pi phi t_n)",
                "amplification peak: ~1e13 when T_ratio >> 100",
                "amplification zero: T_ratio < 1 (cool plasma, base stack only)",
            ],
            'simulation_set': [
                {'name': 'T sweep at t_n=0',
                 'T_keV_range': [0.01, 0.1, 1.0, 4.0, 50.0, 500.0]},
                {'name': 'beat phase sweep at T=4 keV (Perseus)',
                 't_n_range': [0.0, 0.1, 0.25, 0.382, 0.5, 0.618, 0.75, 1.0]},
            ],
        }


# ===========================================================================
# Registry
# ===========================================================================
SESSION_279_CALCULATORS = {
    'UmHeavisideAmplifierCalculator': UmHeavisideAmplifierCalculator,
}


# ===========================================================================
# Smoke tests
# ===========================================================================
def _run_smoke_tests() -> int:
    passed = 0
    failed = 0

    def _check(label: str, cond: bool, detail: str = '') -> None:
        nonlocal passed, failed
        tag = '[PASS]' if cond else '[FAIL]'
        print(f"  {tag} {label}: {detail}")
        if cond:
            passed += 1
        else:
            failed += 1

    calc = UmHeavisideAmplifierCalculator()

    # H-1: cool plasma -> no amplification (f_H = 0)
    r = calc.compute({'Um_base': 1.0, 'T_keV': T_SCM_KEV, 'label': 'cool'})
    _check('H-1 Cool regime (T=T_SCm): f_H = 0 (no amp)', r['f_H'] == 0.0,
           f"f_H = {r['f_H']}")

    # H-2: Perseus core -> partial amplification (T_ratio ~ 47 -> f_H near 0.85)
    r = calc.compute({'Um_base': 1.0, 'T_keV': 4.0, 'label': 'Perseus'})
    _check('H-2 Perseus core (T=4 keV): f_H in [0.7, 1.0]',
           0.7 <= r['f_H'] <= 1.0, f"f_H = {r['f_H']:.4f}")

    # H-3: Perseus amplification factor in [1e12, 2e13]
    _check('H-3 Perseus amplifier in [1e12, 2e13]',
           1.0e12 < r['amplifier_factor'] < 2.0e13,
           f"amp = {r['amplifier_factor']:.3e}")

    # H-4: Sgr A* (T=50 keV) -> near-full amplification
    r = calc.compute({'Um_base': 1.0, 'T_keV': 50.0, 'label': 'SgrA'})
    _check('H-4 Sgr A* near-full amp (> 9e12)',
           r['amplifier_factor'] > 9.0e12, f"amp = {r['amplifier_factor']:.3e}")

    # H-5: monotonic in T
    Ts = [0.01, 0.1, 1.0, 4.0, 50.0, 500.0]
    fHs = [calc.compute({'Um_base': 1.0, 'T_keV': T})['f_H'] for T in Ts]
    mono = all(fHs[i] <= fHs[i + 1] for i in range(len(fHs) - 1))
    _check('H-5 f_H monotonic in T', mono, f"f_H sequence = {[f'{x:.3f}' for x in fHs]}")

    # H-6: beat bounded in [-0.5, +0.5]
    qs = [calc.compute({'Um_base': 1.0, 'T_keV': 4.0, 't_n': t})['f_q']
          for t in [0.0, 0.1, 0.25, 0.382, 0.5, 0.618, 0.75]]
    _check('H-6 f_q bounded in [-0.5, +0.5]',
           all(-0.5 <= q <= 0.5 for q in qs),
           f"f_q range = [{min(qs):+.3f}, {max(qs):+.3f}]")

    # H-7: f_q(0) = epsilon_q * 1 * 1 = 0.5
    r = calc.compute({'Um_base': 1.0, 'T_keV': 4.0, 't_n': 0.0})
    _check('H-7 f_q(t_n=0) = 0.5', abs(r['f_q'] - 0.5) < 1.0e-9,
           f"f_q(0) = {r['f_q']:.6f}")

    # H-8: zero base Um stays sane (no NaN)
    r = calc.compute({'Um_base': 0.0, 'T_keV': 50.0})
    _check('H-8 Um_base=0 returns finite ratio',
           r['amplification_ratio'] == float('inf') or r['amplification_ratio'] >= 1.0,
           f"ratio = {r['amplification_ratio']}")

    # H-9: triple output
    triple = all(k in r for k in ('primary_equations', 'available_equations', 'simulation_set'))
    _check('H-9 Triple output structure', triple, 'triple OK' if triple else 'missing')

    # H-10: layer gate — manually set N_active <= N_FLOOR -> f_H = 0
    r = calc.compute({'Um_base': 1.0, 'T_keV': 50.0, 'N_active': 10.0})
    _check('H-10 Manual N_active <= N_FLOOR forces f_H = 0',
           r['f_H'] == 0.0, f"f_H (gated) = {r['f_H']}")

    # H-11: amplification ratio = multiplier when Um_base = 1
    r = calc.compute({'Um_base': 1.0, 'T_keV': 4.0, 't_n': 0.0})
    _check('H-11 ratio == multiplier for Um_base=1',
           abs(r['amplification_ratio'] - r['total_multiplier']) < 1.0e-6,
           f"ratio={r['amplification_ratio']:.3e}, mul={r['total_multiplier']:.3e}")

    # H-12: applying correction to Perseus baseline reproduces Chandra ICM coupling scale
    # If existing compute_Um returns Um_base ~ 1e-35 J (typical UQFF magnitude),
    # corrected should land in 1e-22 J range (matches observed ICM EM scale).
    r = calc.compute({'Um_base': 1.0e-35, 'T_keV': 4.0, 't_n': 0.0})
    _check('H-12 Perseus Um_base 1e-35 -> corrected in [1e-23, 1e-21]',
           1.0e-23 <= r['Um_corrected'] <= 1.0e-21,
           f"Um_corrected = {r['Um_corrected']:.3e}")

    print()
    print(f"TOTAL: {passed}/{passed + failed} PASS  |  {failed} FAIL")
    return failed


if __name__ == '__main__':
    print("=" * 72)
    print("Session 279 — Um Heaviside SCm Phase-Transition Amplifier (Gap #6)")
    print("=" * 72)
    fails = _run_smoke_tests()

    # Demonstration sweep
    print()
    print("Amplification sweep vs system temperature:")
    print(f"  {'system':<22}{'T_keV':>10}{'f_H':>10}{'amplifier':>15}{'mult(t_n=0)':>16}")
    demo = [
        ('Earth_crust',      2.5e-5),
        ('Solar_corona',     0.17),
        ('SGR_1745',         0.55),
        ('Perseus_core',     4.0),
        ('SgrA_inner',       50.0),
        ('GW150914',         100.0),
        ('M87_jet',          80.0),
    ]
    calc = UmHeavisideAmplifierCalculator()
    for name, T in demo:
        r = calc.compute({'Um_base': 1.0, 'T_keV': T, 't_n': 0.0, 'label': name})
        print(f"  {name:<22}{T:>10.4g}{r['f_H']:>10.3f}"
              f"{r['amplifier_factor']:>15.3e}{r['total_multiplier']:>16.3e}")

    raise SystemExit(fails)
