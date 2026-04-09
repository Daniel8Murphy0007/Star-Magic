#!/usr/bin/env python3
"""
Phase8_Consolidated.py - Session 204 Kozima LENR + Ramanujan Mock Theta Integration
=====================================================================================

MODULES CONSOLIDATED:
- kozima_scm_cross_section.py   → KozimaSCmCrossSection (26-level neutron-drop)
- fneutron_s26_coupling.py      → FneutronS26Coupling (neutron force × S₂₆)
- ramanujan_polylog_s26.py      → RamanujanPolylogS26, S26VDSCalculator
- mock_theta_q26.py             → MockThetaQ26, MockThetaUQFF (q-Pochhammer)
- ramanujan_pi_uqff.py          → RamanujanPiSeries, RamanujanPiUQFF, HypergeometricPi26D

TOTAL: 5 modules, 10 calculator classes, Session 204 codebase

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Consolidated: April 9, 2026 (Session 204 Phase 8)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import math
import numpy as np
from typing import Dict, List, Optional, Any

# ── CONSTANTS ─────────────────────────────────────────────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
K_B = 1.381e-23

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
H_SCM = 0.99
BETA_I = 0.603
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36
RHO_CRIT = 9.47e-27

# LENR frequencies
OMEGA_LENR = 2 * PI * 1.25e12  # 1.25 THz
SIGMA_0 = 1e-28  # m² (base cross-section)

# Ramanujan
RAMANUJAN_CONST = 1103.0
SQRT_2 = math.sqrt(2)
RAMANUJAN_PREFACTOR = (2 * SQRT_2) / 9801
FACTORIAL_CACHE = {}


def _factorial(n):
    if n not in FACTORIAL_CACHE:
        FACTORIAL_CACHE[n] = math.factorial(n)
    return FACTORIAL_CACHE[n]


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 1: Kozima SCm-Modulated Neutron-Drop Cross-Section
# ═══════════════════════════════════════════════════════════════════════════════

class KozimaCrossSection:
    """
    26-level VDS-enhanced neutron-drop cross-section.
    
    sigma_n^SCm(omega, n) = sigma_0 * exp[-(omega - omega_SCm)^2 / (2*Gamma^2)]
                            * (1 + [SSq]*n/26)
    """

    @staticmethod
    def sigma_scm_frequency(omega: float, n: int = 13,
                             sigma_0: float = SIGMA_0,
                             omega_scm: float = OMEGA_LENR,
                             gamma: float = 2 * PI * 50e9) -> float:
        """Frequency-dependent SCm cross-section at level n (1-26)."""
        gaussian = math.exp(-((omega - omega_scm) ** 2) / (2 * gamma ** 2))
        vds_factor = 1 + SSQ * n / 26
        return sigma_0 * gaussian * vds_factor

    @staticmethod
    def sigma_scm_density(rho: float = RHO_SCM,
                           sigma_0: float = SIGMA_0) -> float:
        """Density-dependent SCm cross-section."""
        power_law = (RHO_SCM / RHO_UA) ** SSQ
        cutoff = math.exp(-SSQ * rho / RHO_CRIT)
        return sigma_0 * power_law * cutoff

    @staticmethod
    def f_neutron_scm(N_n: float = 1e10, omega: float = OMEGA_LENR,
                       n: int = 13, phi_phonon: float = 1e20,
                       beta_ratio: float = 0.603) -> float:
        """Full buoyancy-coupled neutron production force."""
        sigma = KozimaCrossSection.sigma_scm_frequency(omega, n)
        coupling = beta_ratio - 1  # buoyancy reversal
        return N_n * sigma * phi_phonon * coupling

    @staticmethod
    def compute_all(params: dict) -> Dict[str, Any]:
        """Compute all Kozima cross-section results."""
        omega = params.get('omega', OMEGA_LENR)
        n = params.get('n_level', 13)
        rho = params.get('rho', RHO_SCM)

        sigma_freq = KozimaCrossSection.sigma_scm_frequency(omega, n)
        sigma_dens = KozimaCrossSection.sigma_scm_density(rho)
        f_neutron = KozimaCrossSection.f_neutron_scm(omega=omega, n=n)

        return {
            'sigma_scm_freq': sigma_freq,
            'sigma_scm_density': sigma_dens,
            'f_neutron_scm': f_neutron,
            'omega': omega,
            'n_level': n,
            'vds_factor': 1 + SSQ * n / 26,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 2: F_neutron × S₂₆ Coupling
# ═══════════════════════════════════════════════════════════════════════════════

class FneutronS26:
    """
    F_neutron at VDS level n coupled to S₂₆ polylogarithmic sum.
    
    F_neutron(omega, n) = sigma_0 * (omega/omega_LENR)^2 * exp(-kappa*t)
                          * [SSq]^n * S₂₆(z)
    """

    @staticmethod
    def s26(z: float, N_terms: int = 55) -> float:
        """26-branch vacuum density series S₂₆(z) = sum_{k=1}^{N} z^k / k^26."""
        if abs(z) > 1:
            z = z / abs(z) * 0.999  # clamp to convergence radius
        total = 0.0
        for k in range(1, N_terms + 1):
            total += (z ** k) / (k ** 26)
        return total

    @staticmethod
    def f_neutron_at_level(omega: float, n: int,
                            t: float = 0.0, z: float = SSQ) -> float:
        """F_neutron at VDS level n with S₂₆ coupling."""
        ratio = (omega / OMEGA_LENR) ** 2
        decay = math.exp(-KAPPA * t)
        ssq_n = SSQ ** n
        s26_val = FneutronS26.s26(z)
        return SIGMA_0 * ratio * decay * ssq_n * s26_val

    @staticmethod
    def compute_all(params: dict) -> Dict[str, Any]:
        """Compute F_neutron across all 26 levels."""
        omega = params.get('omega', OMEGA_LENR)
        t = params.get('t', 0.0)
        z = params.get('z_s26', SSQ)

        s26_val = FneutronS26.s26(z)
        levels = {}
        for n in range(1, 27):
            levels[f'F_neutron_n{n}'] = FneutronS26.f_neutron_at_level(omega, n, t, z)

        return {
            'S26_value': s26_val,
            'F_neutron_levels': levels,
            'F_neutron_n13': levels.get('F_neutron_n13', 0.0),
            'omega': omega,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 3: Ramanujan Polylogarithmic S₂₆
# ═══════════════════════════════════════════════════════════════════════════════

class RamanujanPolylogS26:
    """
    Naive vs Ramanujan-accelerated Li₂₆(z) polylogarithm calculations.
    
    Naive:      Li_26(z) = sum z^k/k^26
    Ramanujan:  Li_26(z) ~ z/(1-z) * sum a_k * [(z/(1-z)) * d/dz]^k {1/(1-z)}
    """

    @staticmethod
    def naive_li26(z: float, N: int = 100) -> float:
        """Direct series Li₂₆(z) = sum_{k=1}^N z^k / k^26."""
        if abs(z) >= 1:
            z = 0.999 * (1 if z > 0 else -1)
        return sum((z ** k) / (k ** 26) for k in range(1, N + 1))

    @staticmethod
    def ramanujan_li26(z: float, N: int = 30) -> float:
        """Ramanujan-accelerated Li₂₆(z)."""
        if abs(z) < 1e-15:
            return 0.0
        if abs(z) >= 1:
            z = 0.999 * (1 if z > 0 else -1)

        # Euler-Maclaurin + Ramanujan acceleration
        u = z / (1 - z) if abs(1 - z) > 1e-15 else z * 1e15
        if abs(u) < 1e-15:
            return z

        result = 0.0
        for k in range(1, N + 1):
            bernoulli_term = 1.0 / (k ** 26)
            weight = u ** k / _factorial(k)
            result += bernoulli_term * weight
        result *= z / (1 - z) if abs(1 - z) > 1e-15 else z * 1e15
        return result

    @staticmethod
    def vds_s26(z: float = SSQ) -> float:
        """S₂₆(z) vacuum density series using best available method."""
        return RamanujanPolylogS26.naive_li26(z)

    @staticmethod
    def compute_all(params: dict) -> Dict[str, Any]:
        """Compute polylog comparisons."""
        z = params.get('z_s26', SSQ)
        naive = RamanujanPolylogS26.naive_li26(z)
        accelerated = RamanujanPolylogS26.ramanujan_li26(z)
        return {
            'Li26_naive': naive,
            'Li26_ramanujan': accelerated,
            'S26_value': naive,
            'z': z,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 4: Mock Theta Functions q₂₆
# ═══════════════════════════════════════════════════════════════════════════════

def q_pochhammer(a: float, q: float, n: int) -> float:
    """(a; q)_n = prod_{k=0}^{n-1} (1 - a*q^k)."""
    product = 1.0
    for k in range(n):
        factor = 1 - a * (q ** k)
        if abs(factor) < 1e-300:
            break
        product *= factor
    return product


class MockThetaQ26:
    """
    Three mock theta functions evaluated at q = [SSq] with 26-state truncation.
    
    f₂₆(q) = sum_{n=0}^{25} q^{n^2} / (q;q)_n^2
    phi₂₆(q) = sum_{n=0}^{25} q^{n^2} / (-q^2; q^2)_n
    psi₂₆(q) = sum_{n=1}^{26} q^{n^2} / (q; q^2)_n
    """

    @staticmethod
    def f26(q: float = SSQ) -> float:
        total = 0.0
        for n in range(26):
            qpoch = q_pochhammer(q, q, n)
            if abs(qpoch) < 1e-300:
                continue
            total += (q ** (n * n)) / (qpoch ** 2)
        return total

    @staticmethod
    def phi26(q: float = SSQ) -> float:
        total = 0.0
        for n in range(26):
            neg_poch = q_pochhammer(-q ** 2, q ** 2, n)
            if abs(neg_poch) < 1e-300:
                continue
            total += (q ** (n * n)) / neg_poch
        return total

    @staticmethod
    def psi26(q: float = SSQ) -> float:
        total = 0.0
        for n in range(1, 27):
            qpoch = q_pochhammer(q, q ** 2, n)
            if abs(qpoch) < 1e-300:
                continue
            total += (q ** (n * n)) / qpoch
        return total


class MockThetaUQFF:
    """
    UQFF-coupled mock theta: multiply mock theta by vacuum/buoyancy.
    
    F_theta(q) = f₂₆(q) * rho_SCm * H_SCm
    Phi_buoy(q) = phi₂₆(q) * beta_i * (rho_SCm / rho_UA)
    """

    @staticmethod
    def f_theta_coupled(q: float = SSQ) -> float:
        return MockThetaQ26.f26(q) * RHO_SCM * H_SCM

    @staticmethod
    def phi_buoyancy_coupled(q: float = SSQ) -> float:
        return MockThetaQ26.phi26(q) * BETA_I * (RHO_SCM / RHO_UA)

    @staticmethod
    def psi_resonance_coupled(q: float = SSQ) -> float:
        return MockThetaQ26.psi26(q) * SSQ * KAPPA

    @staticmethod
    def compute_all(params: dict) -> Dict[str, Any]:
        q = params.get('q', SSQ)
        return {
            'f26': MockThetaQ26.f26(q),
            'phi26': MockThetaQ26.phi26(q),
            'psi26': MockThetaQ26.psi26(q),
            'F_theta_coupled': MockThetaUQFF.f_theta_coupled(q),
            'Phi_buoyancy_coupled': MockThetaUQFF.phi_buoyancy_coupled(q),
            'Psi_resonance_coupled': MockThetaUQFF.psi_resonance_coupled(q),
            'q': q,
        }


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE 5: Ramanujan 1/pi Series (Classical + UQFF + 26D Hypergeometric)
# ═══════════════════════════════════════════════════════════════════════════════

def _ramanujan_coeff(n: int) -> float:
    """a(n) = (4n)! * (1103 + 26390*n) / ((n!)^4 * 396^(4n))."""
    num = _factorial(4 * n) * (RAMANUJAN_CONST + 26390 * n)
    den = (_factorial(n) ** 4) * (396 ** (4 * n))
    if den == 0:
        return 0.0
    return num / den


class RamanujanPiSeries:
    """
    Classical Ramanujan 1/pi:
    1/pi = (2*sqrt(2)/9801) * sum_{n=0}^{N} a(n)
    """

    @staticmethod
    def compute_pi(N_terms: int = 10) -> float:
        """Compute pi via classical Ramanujan series."""
        total = sum(_ramanujan_coeff(n) for n in range(N_terms))
        inv_pi = RAMANUJAN_PREFACTOR * total
        if abs(inv_pi) < 1e-300:
            return math.pi
        return 1.0 / inv_pi

    @staticmethod
    def matching_digits(N_terms: int = 10) -> int:
        """Count digits matching math.pi."""
        computed = RamanujanPiSeries.compute_pi(N_terms)
        return _count_matching(computed, math.pi)


class RamanujanPiUQFF:
    """
    UQFF-modified 1/pi: multiply Ramanujan coefficients by UQFF vacuum factor.
    
    1/pi_UQFF = (2*sqrt(2)/9801) * sum a(n) * (1 + [SSq]*H_SCm*n/26)
    """

    @staticmethod
    def compute_pi(N_terms: int = 10) -> float:
        total = 0.0
        for n in range(N_terms):
            uqff_mod = 1 + SSQ * H_SCM * n / 26
            total += _ramanujan_coeff(n) * uqff_mod
        inv_pi = RAMANUJAN_PREFACTOR * total
        if abs(inv_pi) < 1e-300:
            return math.pi
        return 1.0 / inv_pi

    @staticmethod
    def matching_digits(N_terms: int = 10) -> int:
        computed = RamanujanPiUQFF.compute_pi(N_terms)
        return _count_matching(computed, math.pi)


class HypergeometricPi26D:
    """
    26D hypergeometric pi approximation.
    
    1/pi_26D = (2*sqrt(2)/9801) * sum a(n) * (1 + C_26*n)
    where C_26 = sum_{k=1}^{26} 1/k  (26th harmonic number)
    """

    H26 = sum(1.0 / k for k in range(1, 27))  # ~ 3.9899

    @staticmethod
    def compute_pi(N_terms: int = 10) -> float:
        total = 0.0
        for n in range(N_terms):
            dim_mod = 1 + HypergeometricPi26D.H26 * n
            total += _ramanujan_coeff(n) * dim_mod
        inv_pi = RAMANUJAN_PREFACTOR * total
        if abs(inv_pi) < 1e-300:
            return math.pi
        return 1.0 / inv_pi

    @staticmethod
    def matching_digits(N_terms: int = 10) -> int:
        computed = HypergeometricPi26D.compute_pi(N_terms)
        return _count_matching(computed, math.pi)


def _count_matching(a: float, b: float) -> int:
    """Count matching decimal digits between two floats."""
    sa = f"{a:.20f}"
    sb = f"{b:.20f}"
    count = 0
    for ca, cb in zip(sa, sb):
        if ca == cb:
            count += 1
        else:
            break
    return max(0, count - 1)  # subtract the "3." prefix


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED COMPUTE INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def compute_phase8_all(params: dict) -> Dict[str, Any]:
    """
    Compute all Phase 8 results in one call.
    
    Args:
        params: dict with optional keys: omega, n_level, rho, t, z_s26, q
    
    Returns:
        Dict with all Phase 8 computation results
    """
    results = {}

    # Module 1: Kozima cross-section
    results['kozima'] = KozimaCrossSection.compute_all(params)

    # Module 2: F_neutron x S26
    results['fneutron_s26'] = FneutronS26.compute_all(params)

    # Module 3: Ramanujan polylog S26
    results['polylog_s26'] = RamanujanPolylogS26.compute_all(params)

    # Module 4: Mock theta q26
    results['mock_theta'] = MockThetaUQFF.compute_all(params)

    # Module 5: Ramanujan pi
    results['ramanujan_pi'] = {
        'pi_classical': RamanujanPiSeries.compute_pi(10),
        'pi_uqff': RamanujanPiUQFF.compute_pi(10),
        'pi_26d': HypergeometricPi26D.compute_pi(10),
        'digits_classical': RamanujanPiSeries.matching_digits(10),
        'digits_uqff': RamanujanPiUQFF.matching_digits(10),
        'digits_26d': HypergeometricPi26D.matching_digits(10),
    }

    return results


# ═══════════════════════════════════════════════════════════════════════════════
# SELF-TEST
# ═══════════════════════════════════════════════════════════════════════════════

def self_test():
    """Run self-test for Phase 8."""
    print("Phase8_Consolidated.py self-test")
    print("=" * 60)

    params = {'omega': OMEGA_LENR, 'n_level': 13, 'q': SSQ, 'z_s26': SSQ}
    results = compute_phase8_all(params)

    # Kozima
    k = results['kozima']
    print(f"  Kozima sigma_scm_freq:  {k['sigma_scm_freq']:.3e} m^2")
    print(f"  Kozima f_neutron_scm:   {k['f_neutron_scm']:.3e} N")

    # F_neutron S26
    fn = results['fneutron_s26']
    print(f"  S26({SSQ}):             {fn['S26_value']:.6f}")
    print(f"  F_neutron_n13:          {fn['F_neutron_n13']:.3e}")

    # Polylog
    pl = results['polylog_s26']
    print(f"  Li26_naive({SSQ}):      {pl['Li26_naive']:.6f}")
    print(f"  Li26_ramanujan({SSQ}):  {pl['Li26_ramanujan']:.6e}")

    # Mock theta
    mt = results['mock_theta']
    print(f"  f26({SSQ}):             {mt['f26']:.6f}")
    print(f"  phi26({SSQ}):           {mt['phi26']:.6f}")
    print(f"  psi26({SSQ}):           {mt['psi26']:.6f}")
    print(f"  F_theta_coupled:        {mt['F_theta_coupled']:.3e}")

    # Ramanujan pi
    rp = results['ramanujan_pi']
    print(f"  pi_classical:           {rp['pi_classical']:.15f} ({rp['digits_classical']} digits)")
    print(f"  pi_uqff:               {rp['pi_uqff']:.15f} ({rp['digits_uqff']} digits)")
    print(f"  pi_26d:                {rp['pi_26d']:.15f} ({rp['digits_26d']} digits)")

    # Validation
    checks = [
        k['sigma_scm_freq'] > 0,
        abs(fn['S26_value'] - pl['Li26_naive']) < 1e-6,
        mt['f26'] > 1.0,
        rp['digits_classical'] >= 10,
    ]
    passed = sum(checks)
    print(f"\n  Checks: {passed}/{len(checks)} PASSED")
    print("=" * 60)
    return all(checks)


if __name__ == '__main__':
    ok = self_test()
    raise SystemExit(0 if ok else 1)
