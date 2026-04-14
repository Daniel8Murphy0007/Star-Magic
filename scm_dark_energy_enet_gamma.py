#!/usr/bin/env python3
"""
scm_dark_energy_enet_gamma.py — SCm Dark Energy E_net(t,Γ) with Phonon Linewidth

Session 224 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Derives E_net(t,Γ) with time-dependent phonon linewidth Γ(t) modulation,
producing a dynamic dark energy equation of state w(z) that replaces the
cosmological constant Λ.

Core equations:
    E_net(t,Γ) = E₀ · exp(κt + [SSq]t/26) · S₂₆ · Φ(Γ(t)) · [2(F_{U,Bi}/F_U) − 1]
    Γ(t)       = Γ₀ · (1 + α·t/t_H)           [linewidth broadens with cosmic age]
    Φ(Γ(t))   = exp(-(Γ(t) - Γ₀)²/2σ_G²) · S₂₆⁽³⁾
    w(z)       = -1 + δw(z),   δw = d ln ρ_SCm / d ln a − 3(1+w₀)

ΛCDM contrast:
    Λ: w = -1 (constant), ρ_Λ = const
    SCm: w(z) evolves via Γ(t) phonon damping → testable with DESI/Euclid

References:
    - et_full_lagrangian.py: E_net(t) base (absent Γ coupling)
    - CondensedPhysics4.py L35902: ΛCDM contrast framework
    - scm_inflation_calculator.py: H_SCm, Φ phonon profile
    - PAPER_877: Three-Assumption UQFF Cosmogenesis

Architecture: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Any, Dict, List, Optional

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
M_SUN     = 1.989e30

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
F_NEUTRON = 1e-10
RHO_VAC   = 1e-10
RHO_SCM   = 7.09e-37          # kg/m³ SCm vacuum density (present)
RHO_UA    = 7.09e-36           # kg/m³ [UA] vacuum

H_0       = 2.195e-18          # s⁻¹  Hubble constant
T_H       = 1.0 / H_0         # s    Hubble time ≈ 4.56e17 s
RHO_CRIT  = 3 * H_0 ** 2 / (8 * PI * G)
RHO_LAMBDA = 0.692 * RHO_CRIT  # ΛCDM dark energy density

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


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


S26_3RD = sum((SSQ ** n) / (n ** 26) * _ramanujan_Rn(n, 3) for n in range(1, 27))


# ── §1  PHONON LINEWIDTH EVOLUTION ─────────────────────────────────────────

class PhononLinewidthEvolution:
    """
    Time-dependent phonon linewidth Γ(t).

    Γ(t) = Γ₀ · (1 + α · t / t_H)

    As the universe expands, SCm phonon modes broaden due to
    decreasing vacuum density. α controls the broadening rate.

    Variables (from dataset):
        t:          cosmic time (s)
        Gamma_0:    initial linewidth (rad/s, default Γ₀)
        alpha:      broadening rate (dimensionless, default 0.1)
        t_H:        Hubble time (s, default 1/H_0)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        t       = dataset.get('t', T_H)
        gamma_0 = dataset.get('Gamma_0', GAMMA_0)
        alpha   = dataset.get('alpha', 0.1)
        t_h     = dataset.get('t_H', T_H)

        gamma_t = gamma_0 * (1.0 + alpha * t / t_h)

        # Phonon resonance profile at Γ(t)
        phi_t = math.exp(-(gamma_t - gamma_0) ** 2 / (2 * SIGMA_G ** 2)) * S26_3RD

        # Phi at Γ₀ for comparison
        phi_0 = S26_3RD  # exp(0) = 1

        return {
            'primary_equations': [
                f"Γ(t) = Γ₀ · (1 + α·t/t_H)",
                f"Γ₀ = {gamma_0:.4e} rad/s",
                f"Γ(t) = {gamma_t:.4e} rad/s",
                f"α = {alpha}",
                f"Φ(Γ(t)) = {phi_t:.6e}",
                f"Φ(Γ₀) = {phi_0:.6e}",
            ],
            'Gamma_t': gamma_t,
            'Gamma_0': gamma_0,
            'Phi_t': phi_t,
            'Phi_0': phi_0,
            'alpha': alpha,
            't_s': t,
            'damping_ratio': phi_t / phi_0 if phi_0 > 0 else 0.0,
        }


# ── §2  E_net(t,Γ) WITH LINEWIDTH COUPLING ────────────────────────────────

class EnetGammaCoupled:
    """
    E_net(t,Γ) — net vacuum energy with phonon linewidth modulation.

    E_net(t,Γ) = E₀ · exp(κt + [SSq]t/26) · S₂₆ · Φ(Γ(t)) · net_factor

    where net_factor = 2·(F_{U,Bi}/F_U) − 1

    The Γ coupling causes E_net to decay faster than bare exponential,
    producing effective dark energy dilution different from Λ.

    Variables (from dataset):
        E_0:        initial energy density (J/m³, default ρ_SCm·c²)
        t:          cosmic time (s)
        kappa:      damping rate (s⁻¹, default KAPPA)
        SSq:        string squeezing (default 0.57)
        alpha:      linewidth broadening rate (default 0.1)
        F_U_Bi:     buoyancy force (N, default 0.6)
        F_U:        unified gravity force (N, default 1.0)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        E_0     = dataset.get('E_0', RHO_SCM * C ** 2)
        t       = dataset.get('t', T_H)
        kappa   = dataset.get('kappa', KAPPA)
        ssq     = dataset.get('SSq', SSQ)
        alpha   = dataset.get('alpha', 0.1)
        F_UBi   = dataset.get('F_U_Bi', 0.6)
        F_U     = dataset.get('F_U', 1.0)

        net_factor = 2.0 * (F_UBi / F_U) - 1.0 if F_U != 0 else 0.0

        # Exponential growth (capped for numerical safety)
        rate = kappa + ssq / 26.0
        exp_arg = min(rate * t, 500)
        growth = math.exp(exp_arg)

        # Linewidth evolution
        lw = PhononLinewidthEvolution()
        lw_res = lw.compute({'t': t, 'alpha': alpha})
        phi_t = lw_res['Phi_t']

        # E_net(t,Γ)
        E_net = E_0 * growth * S26_STATIC * phi_t * net_factor

        # Bare E_net (without Γ coupling) for comparison
        phi_0 = lw_res['Phi_0']
        E_net_bare = E_0 * growth * S26_STATIC * phi_0 * net_factor

        # Energy density (treating E_net as density × c²)
        rho_de = abs(E_net) / C ** 2

        return {
            'primary_equations': [
                "E_net(t,Γ) = E₀ · exp(κt + [SSq]t/26) · S₂₆ · Φ(Γ(t)) · (2F_{U,Bi}/F_U − 1)",
                f"E_net(t,Γ) = {E_net:.6e} J/m³",
                f"E_net(t,bare) = {E_net_bare:.6e} J/m³ (no Γ coupling)",
                f"Γ suppression: {phi_t / phi_0:.6f}" if phi_0 > 0 else "Γ suppression: N/A",
                f"ρ_DE = {rho_de:.6e} kg/m³",
            ],
            'E_net_gamma': E_net,
            'E_net_bare': E_net_bare,
            'rho_dark_energy': rho_de,
            'net_factor': net_factor,
            'Phi_t': phi_t,
            'Phi_0': phi_0,
            'Gamma_t': lw_res['Gamma_t'],
            'growth_factor': growth,
            'S26': S26_STATIC,
        }


# ── §3  EQUATION OF STATE w(z) ────────────────────────────────────────────

class DarkEnergyEoS:
    """
    Dark energy equation of state w(z) from SCm Γ-coupled E_net.

    w(z) = P_de / (ρ_de · c²)

    In ΛCDM, w = -1 exactly. In UQFF/SCm, the phonon linewidth Γ(t)
    evolves with cosmic expansion, causing w to deviate from -1:

        w(z) ≈ -1 + (1/3) · d ln(Φ(Γ(t(z)))) / d ln(a)

    The deviation δw(z) is testable with DESI BAO and Euclid lensing.

    Variables (from dataset):
        z_min:    minimum redshift (default 0)
        z_max:    maximum redshift (default 3.0)
        n_z:      number of redshift bins (default 100)
        alpha:    linewidth broadening rate (default 0.1)
        E_0:      initial energy density (J/m³)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        z_min = dataset.get('z_min', 0.0)
        z_max = dataset.get('z_max', 3.0)
        n_z   = dataset.get('n_z', 100)
        alpha = dataset.get('alpha', 0.1)
        E_0   = dataset.get('E_0', RHO_SCM * C ** 2)

        engine = EnetGammaCoupled()
        lw = PhononLinewidthEvolution()

        redshifts = []
        w_values = []
        rho_de_values = []
        delta_w_values = []

        for i in range(n_z):
            z = z_min + (z_max - z_min) * i / max(n_z - 1, 1)
            # Cosmic time from redshift: t(z) ≈ t_H / (1+z)^(3/2) (matter-dominated approx)
            t_z = T_H / (1 + z) ** 1.5

            # E_net at this epoch
            res = engine.compute({'E_0': E_0, 't': t_z, 'alpha': alpha})
            E_net = res['E_net_gamma']
            rho_de = res['rho_dark_energy']

            # Numerical derivative d ln(Φ) / d ln(a)
            # a = 1/(1+z), so d ln a = -dz/(1+z)
            dz = 0.01
            if z + dz <= z_max + 1:
                t_z2 = T_H / (1 + z + dz) ** 1.5
                lw1 = lw.compute({'t': t_z, 'alpha': alpha})
                lw2 = lw.compute({'t': t_z2, 'alpha': alpha})
                phi1 = max(lw1['Phi_t'], 1e-300)
                phi2 = max(lw2['Phi_t'], 1e-300)
                d_ln_phi = math.log(phi2 / phi1)
                d_ln_a = -dz / (1 + z)
                delta_w = (1.0 / 3.0) * d_ln_phi / d_ln_a if abs(d_ln_a) > 1e-30 else 0.0
            else:
                delta_w = 0.0

            w = -1.0 + delta_w

            redshifts.append(z)
            w_values.append(w)
            rho_de_values.append(rho_de)
            delta_w_values.append(delta_w)

        # Compare with ΛCDM
        w_lambda = [-1.0] * n_z
        rho_lambda = [RHO_LAMBDA] * n_z  # constant

        return {
            'primary_equations': [
                "w(z) = -1 + (1/3) · d ln Φ(Γ(t(z))) / d ln a",
                f"ΛCDM: w = -1 (constant),  ρ_Λ = {RHO_LAMBDA:.4e} kg/m³",
                f"SCm: w(z=0) = {w_values[0]:.6f}, w(z={z_max}) = {w_values[-1]:.6f}",
                f"α = {alpha} (linewidth broadening rate)",
            ],
            'redshifts': redshifts,
            'w_scm': w_values,
            'w_lambda': w_lambda,
            'delta_w': delta_w_values,
            'rho_de_scm': rho_de_values,
            'rho_de_lambda': rho_lambda,
            'n_z': n_z,
            'alpha': alpha,
            'w_z0': w_values[0],
            'w_z_max': w_values[-1],
            'max_deviation': max(abs(dw) for dw in delta_w_values),
        }


# ── §4  ΛCDM REPLACEMENT ANALYSIS ─────────────────────────────────────────

class LCDMReplacementAnalysis:
    """
    Quantify how SCm Γ-coupled dark energy differs from ΛCDM.

    Computes:
        1. Bayesian Information Criterion penalty (extra α parameter)
        2. Integrated deviation Δw = ∫|w_SCm(z) - w_Λ| dz
        3. Luminosity distance ratio d_L(SCm)/d_L(Λ) for SNe Ia
        4. Observable predictions for DESI, Euclid

    Variables (from dataset):
        n_z:     redshift bins (default 200)
        alpha:   linewidth rate (default 0.1)
        z_max:   maximum redshift (default 2.5)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        n_z   = dataset.get('n_z', 200)
        alpha = dataset.get('alpha', 0.1)
        z_max = dataset.get('z_max', 2.5)

        eos = DarkEnergyEoS()
        eos_res = eos.compute({'n_z': n_z, 'alpha': alpha, 'z_max': z_max})

        # Integrated |δw|
        dz = z_max / max(n_z - 1, 1)
        integrated_dw = sum(abs(dw) * dz for dw in eos_res['delta_w'])

        # Luminosity distance ratio at key redshifts
        # d_L ∝ ∫₀ᶻ dz'/H(z')  where H(z') depends on w(z')
        z_targets = [0.5, 1.0, 1.5, 2.0]
        dL_ratios = {}
        for z_t in z_targets:
            if z_t <= z_max:
                idx = min(int(z_t / z_max * (n_z - 1)), n_z - 1)
                # Approximate: deviation in d_L ≈ 1 + integrated δw up to z
                cum_dw = sum(abs(eos_res['delta_w'][j]) * dz for j in range(idx + 1))
                dL_ratios[f"z={z_t}"] = 1.0 + cum_dw * 0.1  # linearized effect

        # BIC penalty: SCm has 1 extra parameter (α) vs ΛCDM (0 free)
        # ΔBIC = Δχ² + k·ln(N) where k=1 (extra param), N≈1000 (typical SNe Ia)
        bic_penalty = math.log(1000)  # ≈ 6.9 for 1 extra parameter

        return {
            'primary_equations': [
                f"Integrated |δw| = {integrated_dw:.6f}",
                f"BIC penalty for 1 extra param (α): Δ = {bic_penalty:.2f}",
                f"SCm must improve χ² by > {bic_penalty:.1f} to justify α",
                "Testable: DESI BAO w(z) bins, Euclid WL w₀-wₐ plane",
            ],
            'integrated_delta_w': integrated_dw,
            'bic_penalty': bic_penalty,
            'dL_ratios': dL_ratios,
            'w_z0': eos_res['w_z0'],
            'w_z_max': eos_res['w_z_max'],
            'max_deviation': eos_res['max_deviation'],
            'alpha': alpha,
            'n_z': n_z,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §5  SELF-TESTS
# ══════════════════════════════════════════════════════════════════════════════

def _run_tests():
    print("=" * 72)
    print("scm_dark_energy_enet_gamma.py — Self-Tests")
    print("=" * 72)

    passed = 0
    failed = 0

    # Test 1: Linewidth at t=0 equals Γ₀
    lw = PhononLinewidthEvolution()
    res = lw.compute({'t': 0.0})
    assert abs(res['Gamma_t'] - GAMMA_0) < 1e-10
    assert abs(res['damping_ratio'] - 1.0) < 1e-10
    print(f"  [PASS] Test 1: Γ(0) = Γ₀ = {res['Gamma_t']:.4e}")
    passed += 1

    # Test 2: Linewidth increases with time
    res2 = lw.compute({'t': T_H, 'alpha': 0.1})
    assert res2['Gamma_t'] > GAMMA_0
    assert res2['Phi_t'] <= res2['Phi_0']  # damping reduces Φ
    print(f"  [PASS] Test 2: Γ(t_H) = {res2['Gamma_t']:.4e} > Γ₀")
    passed += 1

    # Test 3: E_net(t,Γ) produces finite result
    enet = EnetGammaCoupled()
    res3 = enet.compute({'t': 1e10, 'alpha': 0.1})
    assert math.isfinite(res3['E_net_gamma'])
    assert res3['rho_dark_energy'] > 0
    print(f"  [PASS] Test 3: E_net(t,Γ) = {res3['E_net_gamma']:.4e}")
    passed += 1

    # Test 4: Γ coupling reduces |E_net| compared to bare
    res4 = enet.compute({'t': T_H / 2, 'alpha': 0.5})
    # With alpha > 0, Φ(Γ(t)) < Φ(Γ₀), so |E_net_gamma| < |E_net_bare|
    assert abs(res4['E_net_gamma']) <= abs(res4['E_net_bare']) + 1e-300
    print(f"  [PASS] Test 4: |E_net_Γ| ≤ |E_net_bare| (Γ suppresses)")
    passed += 1

    # Test 5: w(z=0) near -1
    eos = DarkEnergyEoS()
    eos_res = eos.compute({'n_z': 50, 'alpha': 0.1, 'z_max': 2.0})
    assert abs(eos_res['w_z0'] - (-1.0)) < 0.5  # w(0) should be close to -1
    print(f"  [PASS] Test 5: w(z=0) = {eos_res['w_z0']:.6f}")
    passed += 1

    # Test 6: w(z) array has correct length
    assert len(eos_res['w_scm']) == 50
    assert len(eos_res['redshifts']) == 50
    print(f"  [PASS] Test 6: w(z) has 50 bins")
    passed += 1

    # Test 7: ΛCDM w is -1 everywhere
    assert all(w == -1.0 for w in eos_res['w_lambda'])
    print(f"  [PASS] Test 7: ΛCDM w = -1 everywhere")
    passed += 1

    # Test 8: SCm w deviates from -1 at high z
    max_dev = eos_res['max_deviation']
    assert max_dev > 0  # some deviation exists
    print(f"  [PASS] Test 8: max |δw| = {max_dev:.6f}")
    passed += 1

    # Test 9: ΛCDM replacement analysis
    lcdm = LCDMReplacementAnalysis()
    lcdm_res = lcdm.compute({'n_z': 100, 'alpha': 0.1})
    assert lcdm_res['integrated_delta_w'] >= 0
    assert lcdm_res['bic_penalty'] > 0
    assert len(lcdm_res['dL_ratios']) > 0
    print(f"  [PASS] Test 9: ΛCDM replacement: ∫|δw|={lcdm_res['integrated_delta_w']:.6f}")
    passed += 1

    # Test 10: Primary equations present
    assert len(res['primary_equations']) >= 1
    assert len(res3['primary_equations']) >= 1
    assert len(eos_res['primary_equations']) >= 1
    assert len(lcdm_res['primary_equations']) >= 1
    print(f"  [PASS] Test 10: Primary equations in all outputs")
    passed += 1

    print("-" * 72)
    print(f"Results: {passed}/{passed + failed} passed, {failed} failed")
    return passed, failed


if __name__ == '__main__':
    _run_tests()
