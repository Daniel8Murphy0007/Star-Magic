#!/usr/bin/env python3
"""
scm_lenr_cop_gamma.py — LENR COP Parametric Engine as f(Gamma)

Session 225 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Maps phonon linewidth Gamma to LENR Coefficient of Performance (COP):

  COP(Gamma) = P_out(Gamma) / P_in(Gamma)

where:
  P_in(Gamma)  = hbar * omega_SCm * Phi(Gamma) * A_cell    [input phonon power]
  P_out(Gamma) = E_dd * R_nd(Gamma) * V_active              [excess heat power]

  R_nd(Gamma)  = sigma_n * n_H * phi_0 * Phi(Gamma) * exp(-E_a / k_B T)
                 × SCm micro-plasmoid formation probability

with Phi(Gamma) = exp(-(Gamma - Gamma_0)^2 / (2 sigma_G^2)) × S_26^{(3)}

Classes:
  NeutronDropRate         — R_nd(Gamma) neutron-assisted transmutation rate
  MicroPlasmoidThreshold  — Gamma threshold for micro-plasmoid ignition
  LENRPowerBalance        — P_out and P_in as f(Gamma)
  COPParametricEngine     — Full COP(Gamma) sweep with Kozima/F-P validation

References:
  - scm_phonon_linewidth.py: Gamma sweep framework
  - positive_et_expansion.py: E_net factor
  - Kozima neutron-drop model (H. Kozima, Cold Fusion, 1998)
  - Fleischmann-Pons excess heat (M. Fleischmann & S. Pons, 1989)
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, Any, List

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
HBAR      = 1.055e-34      # J·s
K_B       = 1.381e-23      # J/K
C         = 2.998e8        # m/s
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0

# SCm phonon
OMEGA_SCM = 2 * PI * 1.25e12    # rad/s (1.25 THz)
GAMMA_0   = 2 * PI * 0.1e12     # optimal linewidth (rad/s)
SIGMA_G   = 0.08 * 2 * PI * 1e12

# LENR cell parameters
SIGMA_N   = 1e-4          # neutron capture cross section (barns → effective m²)
N_H       = 6e28           # hydrogen density in Pd lattice (m⁻³)
PHI_0     = 1e20           # background phonon fluence (m⁻² s⁻¹)
E_A       = 0.3 * 1.602e-19  # activation energy ~0.3 eV (J)
T_CELL    = 350.0          # cell temperature (K)
E_DD      = 23.8e6 * 1.602e-19  # D-D fusion energy ~23.8 MeV (J)
V_ACTIVE  = 1e-6           # active volume (1 cm³ → m³)
A_CELL    = 1e-4           # cell area (1 cm² → m²)

# F_NEUTRON constant
F_NEUTRON = 1e-10


def _ramanujan_Rn(n: int, k: int = 3) -> float:
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


S26_3RD = sum(SSQ**n / n**26 * _ramanujan_Rn(n, 3) for n in range(1, 27))


def _phonon_fluence(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# ═══════════════════════════════════════════════════════════════════════════════
# §1  NEUTRON-DROP RATE  R_nd(Gamma)
# ═══════════════════════════════════════════════════════════════════════════════

class NeutronDropRate:
    """SCm neutron-assisted transmutation rate.

    R_nd(Gamma) = sigma_n · n_H · phi_0 · Phi(Gamma) · exp(-E_a / k_B T)
                  × P_plasmoid(Gamma)

    where P_plasmoid is the micro-plasmoid formation probability.
    """

    def compute(self, dataset: dict) -> dict:
        gamma = float(dataset.get('gamma', GAMMA_0))
        T = float(dataset.get('T_cell', T_CELL))
        sigma_n = float(dataset.get('sigma_n', SIGMA_N))
        n_H = float(dataset.get('n_H', N_H))
        phi_0 = float(dataset.get('phi_0', PHI_0))
        E_a = float(dataset.get('E_a', E_A))

        Phi = _phonon_fluence(gamma)
        boltzmann = math.exp(-E_a / (K_B * T))

        # Micro-plasmoid formation: peaks at Gamma_0
        P_plasmoid = Phi / S26_3RD if S26_3RD > 0 else 0  # normalized [0,1]

        R_nd = sigma_n * n_H * phi_0 * Phi * boltzmann * P_plasmoid

        return {
            'gamma_rad': gamma,
            'gamma_THz': gamma / (2 * PI * 1e12),
            'Phi': Phi,
            'boltzmann_factor': boltzmann,
            'P_plasmoid': P_plasmoid,
            'R_nd': R_nd,
            'primary_equations': [
                'Neutron-Drop Rate:',
                '  R_nd = sigma_n · n_H · phi_0 · Phi(Gamma) · exp(-E_a/k_BT) · P_pl',
                f'  Gamma = {gamma/(2*PI*1e12):.4f} THz',
                f'  Phi(Gamma) = {Phi:.6e}',
                f'  Boltzmann = exp(-{E_a/1.602e-19:.2f} eV / {K_B*T/1.602e-19:.4f} eV) = {boltzmann:.6e}',
                f'  P_plasmoid = {P_plasmoid:.6f}',
                f'  R_nd = {R_nd:.6e} reactions/m^3/s',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §2  MICRO-PLASMOID THRESHOLD
# ═══════════════════════════════════════════════════════════════════════════════

class MicroPlasmoidThreshold:
    """Determine the linewidth Gamma threshold for micro-plasmoid ignition.

    A micro-plasmoid forms when Phi(Gamma) exceeds a critical threshold:
      Phi(Gamma) > Phi_crit = F_neutron / (sigma_n · n_H · V_active)

    This defines a range [Gamma_min, Gamma_max] around Gamma_0.
    """

    def compute(self, dataset: dict) -> dict:
        sigma_n = float(dataset.get('sigma_n', SIGMA_N))
        n_H = float(dataset.get('n_H', N_H))
        V_active = float(dataset.get('V_active', V_ACTIVE))

        Phi_crit = F_NEUTRON / (sigma_n * n_H * V_active) if sigma_n * n_H * V_active > 0 else float('inf')

        # Phi(Gamma) = exp(-(G-G0)^2/(2 sig^2)) × S26
        # Phi > Phi_crit  →  exp(...) > Phi_crit / S26
        ratio = Phi_crit / S26_3RD if S26_3RD > 0 else float('inf')

        if 0 < ratio < 1:
            delta = math.sqrt(-2 * SIGMA_G**2 * math.log(ratio))
            gamma_min = GAMMA_0 - delta
            gamma_max = GAMMA_0 + delta
            width_THz = 2 * delta / (2 * PI * 1e12)
            ignition = True
        else:
            gamma_min = gamma_max = GAMMA_0
            width_THz = 0.0
            ignition = ratio <= 1  # always on if ratio <= 0

        return {
            'Phi_crit': Phi_crit,
            'Phi_peak': S26_3RD,
            'ratio': ratio,
            'ignition_possible': ignition,
            'gamma_min_THz': gamma_min / (2 * PI * 1e12),
            'gamma_max_THz': gamma_max / (2 * PI * 1e12),
            'ignition_width_THz': width_THz,
            'primary_equations': [
                'Micro-Plasmoid Threshold:',
                '  Phi(Gamma) > Phi_crit = F_neutron / (sigma_n · n_H · V)',
                f'  Phi_crit = {Phi_crit:.6e}',
                f'  Phi_peak = S_26^(3) = {S26_3RD:.6e}',
                f'  Ratio = {ratio:.6e}',
                f'  Ignition possible: {ignition}',
                f'  Gamma range: [{gamma_min/(2*PI*1e12):.4f}, {gamma_max/(2*PI*1e12):.4f}] THz',
                f'  Width: {width_THz:.4f} THz',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §3  LENR POWER BALANCE
# ═══════════════════════════════════════════════════════════════════════════════

class LENRPowerBalance:
    """Input and output power as functions of Gamma.

    P_in(Gamma) = hbar · omega_SCm · phi_0 · Phi(Gamma) · A_cell
    P_out(Gamma) = E_dd · R_nd(Gamma) · V_active
    """

    def compute(self, dataset: dict) -> dict:
        gamma = float(dataset.get('gamma', GAMMA_0))
        E_dd = float(dataset.get('E_dd', E_DD))
        V_active = float(dataset.get('V_active', V_ACTIVE))
        A_cell = float(dataset.get('A_cell', A_CELL))
        phi_0 = float(dataset.get('phi_0', PHI_0))

        Phi = _phonon_fluence(gamma)

        # Input power: phonon pumping
        P_in = HBAR * OMEGA_SCM * phi_0 * Phi * A_cell

        # Output power: fusion reactions
        rnd_result = NeutronDropRate().compute(dataset)
        R_nd = rnd_result['R_nd']
        P_out = E_dd * R_nd * V_active

        COP = P_out / P_in if P_in > 0 else 0.0

        return {
            'gamma_THz': gamma / (2 * PI * 1e12),
            'Phi': Phi,
            'P_in_W': P_in,
            'P_out_W': P_out,
            'COP': COP,
            'R_nd': R_nd,
            'primary_equations': [
                'LENR Power Balance:',
                '  P_in = hbar · omega_SCm · phi_0 · Phi · A_cell',
                '  P_out = E_dd · R_nd · V_active',
                '  COP = P_out / P_in',
                f'  Gamma = {gamma/(2*PI*1e12):.4f} THz',
                f'  P_in = {P_in:.6e} W',
                f'  P_out = {P_out:.6e} W',
                f'  COP = {COP:.4f}',
            ],
        }


# ═══════════════════════════════════════════════════════════════════════════════
# §4  COP PARAMETRIC ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

class COPParametricEngine:
    """Full COP(Gamma) sweep with validation against Kozima/F-P data."""

    def compute(self, dataset: dict) -> dict:
        gamma_sweep_THz = dataset.get('gamma_sweep_THz',
            [0.01, 0.03, 0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50])

        results = []
        peak_cop = 0.0
        peak_gamma = 0.0
        for g_THz in gamma_sweep_THz:
            gamma = 2 * PI * g_THz * 1e12
            ds = dict(dataset)
            ds['gamma'] = gamma
            bal = LENRPowerBalance().compute(ds)
            results.append({
                'gamma_THz': g_THz,
                'P_in_W': bal['P_in_W'],
                'P_out_W': bal['P_out_W'],
                'COP': bal['COP'],
                'R_nd': bal['R_nd'],
            })
            if bal['COP'] > peak_cop:
                peak_cop = bal['COP']
                peak_gamma = g_THz

        # Threshold and ignition
        threshold = MicroPlasmoidThreshold().compute(dataset)

        # Kozima model comparison (predicted COP ~ 1.2-3 for Pd-D)
        kozima_range = (1.2, 3.0)
        # Fleischmann-Pons claimed excess: COP ~ 1.1-10
        fp_range = (1.1, 10.0)

        cop_in_kozima = kozima_range[0] <= peak_cop <= kozima_range[1]
        cop_in_fp = fp_range[0] <= peak_cop <= fp_range[1]

        return {
            'sweep': results,
            'peak_COP': peak_cop,
            'peak_gamma_THz': peak_gamma,
            'threshold': threshold,
            'kozima_compatible': cop_in_kozima,
            'fp_compatible': cop_in_fp,
            'kozima_range': kozima_range,
            'fp_range': fp_range,
            'primary_equations': [
                'COP Parametric Engine:',
                f'  Peak COP = {peak_cop:.4f} at Gamma = {peak_gamma:.4f} THz',
                '',
                'Gamma sweep:',
            ] + [
                f'  Gamma={r["gamma_THz"]:.2f} THz: COP={r["COP"]:.4f}, '
                f'P_in={r["P_in_W"]:.2e}, P_out={r["P_out_W"]:.2e}'
                for r in results
            ] + [
                '',
                f'  Kozima model (COP 1.2-3): {"COMPATIBLE" if cop_in_kozima else "OUTSIDE RANGE"}',
                f'  F-P excess (COP 1.1-10): {"COMPATIBLE" if cop_in_fp else "OUTSIDE RANGE"}',
                '',
                f'  Ignition window: {threshold["gamma_min_THz"]:.4f}-{threshold["gamma_max_THz"]:.4f} THz',
            ],
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    print("=" * 72)
    print("scm_lenr_cop_gamma.py — Self-Tests")
    print("=" * 72)

    ok = True
    passed = 0

    # Test 1: R_nd positive at optimal gamma
    rnd = NeutronDropRate().compute({'gamma': GAMMA_0})
    if rnd['R_nd'] > 0 and math.isfinite(rnd['R_nd']):
        print(f"  [PASS] Test 1: R_nd(Gamma_0) = {rnd['R_nd']:.4e} reactions/m^3/s")
        passed += 1
    else:
        print(f"  [FAIL] Test 1: R_nd = {rnd['R_nd']}"); ok = False

    # Test 2: R_nd peaks at Gamma_0
    rnd_low = NeutronDropRate().compute({'gamma': 0.01 * 2 * PI * 1e12})
    rnd_peak = NeutronDropRate().compute({'gamma': GAMMA_0})
    if rnd_peak['R_nd'] > rnd_low['R_nd']:
        print(f"  [PASS] Test 2: R_nd peaks at Gamma_0 ({rnd_peak['R_nd']:.2e} > {rnd_low['R_nd']:.2e})")
        passed += 1
    else:
        print(f"  [FAIL] Test 2: Wrong peak location"); ok = False

    # Test 3: Micro-plasmoid ignition possible
    thr = MicroPlasmoidThreshold().compute({})
    if thr['ignition_possible']:
        print(f"  [PASS] Test 3: Ignition possible, width = {thr['ignition_width_THz']:.4f} THz")
        passed += 1
    else:
        print(f"  [FAIL] Test 3: No ignition"); ok = False

    # Test 4: Power balance finite
    bal = LENRPowerBalance().compute({'gamma': GAMMA_0})
    if math.isfinite(bal['P_in_W']) and math.isfinite(bal['P_out_W']) and bal['P_in_W'] > 0:
        print(f"  [PASS] Test 4: P_in = {bal['P_in_W']:.4e} W, P_out = {bal['P_out_W']:.4e} W")
        passed += 1
    else:
        print(f"  [FAIL] Test 4: Non-finite power"); ok = False

    # Test 5: COP > 0 at optimal gamma
    if bal['COP'] > 0:
        print(f"  [PASS] Test 5: COP = {bal['COP']:.4f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 5: COP <= 0"); ok = False

    # Test 6: Full parametric sweep has 10 points
    cop = COPParametricEngine().compute({})
    if len(cop['sweep']) == 10:
        print(f"  [PASS] Test 6: Sweep has {len(cop['sweep'])} points")
        passed += 1
    else:
        print(f"  [FAIL] Test 6: Expected 10 points, got {len(cop['sweep'])}"); ok = False

    # Test 7: Peak COP at Gamma ≈ 0.10 THz
    if abs(cop['peak_gamma_THz'] - 0.10) < 0.05:
        print(f"  [PASS] Test 7: Peak at Gamma = {cop['peak_gamma_THz']:.4f} THz")
        passed += 1
    else:
        print(f"  [FAIL] Test 7: Peak at {cop['peak_gamma_THz']} THz (expected ~0.10)"); ok = False

    # Test 8: COP monotonically decreasing away from peak
    sweep = cop['sweep']
    peak_idx = next(i for i, r in enumerate(sweep) if r['gamma_THz'] == cop['peak_gamma_THz'])
    # Check right side decreases
    right_side = [s['COP'] for s in sweep[peak_idx:]]
    right_decr = all(right_side[i] >= right_side[i+1] for i in range(len(right_side)-1))
    if right_decr:
        print(f"  [PASS] Test 8: COP decreases away from peak (right side)")
        passed += 1
    else:
        print(f"  [FAIL] Test 8: COP not monotone right of peak"); ok = False

    # Test 9: Primary equations in all outputs
    all_have = all(
        'primary_equations' in r
        for r in [rnd, thr, bal, cop]
    )
    if all_have:
        print(f"  [PASS] Test 9: All outputs contain primary_equations")
        passed += 1
    else:
        print(f"  [FAIL] Test 9: Missing primary_equations"); ok = False

    # Test 10: COP at off-peak is less than peak COP
    cop_far = LENRPowerBalance().compute({'gamma': 0.50 * 2 * PI * 1e12})
    if cop_far['COP'] < bal['COP']:
        print(f"  [PASS] Test 10: COP(0.50 THz)={cop_far['COP']:.4e} < COP(0.10 THz)={bal['COP']:.4f}")
        passed += 1
    else:
        print(f"  [FAIL] Test 10: Off-peak COP not less"); ok = False

    print(f"\n  scm_lenr_cop_gamma.py: {passed}/10 tests passed")
    return ok


if __name__ == "__main__":
    success = _run_tests()
    raise SystemExit(0 if success else 1)
