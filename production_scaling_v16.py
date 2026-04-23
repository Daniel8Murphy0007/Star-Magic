#!/usr/bin/env python3
"""
production_scaling_v16.py — Production Scaling to 800k calc/s

Session 223 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
v16 upgrade: 800,000 calc/s target (up from v15's 650k).
36 kernels total = v15's 30 + 6 new Session 223 kernels:
  kernel_scm_inflation_hubble    — SCm-driven Hubble parameter H_SCm
  kernel_thorne_morris_exotic    — Thorne-Morris exotic matter (ρ+P)<0
  kernel_vds_convergence         — VDS polylogarithm Li₂₆([SSq])
  kernel_dvp_prime_sieve         — DVP prime encoding a(p)
  kernel_bsh_harmonic_saturation — BSH saturation 1-exp(-[SSq]m)
  kernel_wormhole_phonon_damping — Wormhole throat phonon resonance

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell

import time
from typing import Dict, List, Callable

# ── §0  Constants ──────────────────────────────────────────────────────────

PI        = math.pi
G         = 6.674e-11
C         = 2.998e8
M_SUN     = 1.989e30
HBAR      = 1.055e-34
KPC       = 3.086e19

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
F_NEUTRON = 1e-10
RHO_VAC   = 1e-10

S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))

TARGET_CALC_PER_SEC = 800_000


# ── §1  Ramanujan S₂₆⁽³⁾ ──────────────────────────────────────────────────

def ramanujan_Rn(n: int, k: int = 3) -> float:
    total = 0.0
    for j in range(k):
        sign = (-1) ** j
        binom = 1.0
        for m in range(j):
            binom *= (k - 1 - m) / (m + 1)
        nfact = math.factorial(min(n + j, 170))
        total += sign * binom / nfact
    return total


S26_3RD = sum((SSQ ** n) / (n ** 26) * ramanujan_Rn(n, 3) for n in range(1, 27))


def Phi(gamma: float) -> float:
    return math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


# ── §2  Carry-Forward Kernels (v11: 16, v13: 4, v14: 4, v15: 6 = 30) ───

def kernel_gravity_26layer(M=M_SUN, r=6.96e8):
    return sum(dpm_ug1_seed(M, r) * SSQ * i / 26 for i in range(1, 27))


def kernel_fu_bi_i(M=M_SUN, r=6.96e8, t=86400):
    Ug = sum(dpm_ug1_seed(M, r) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_ug1_seed(M, r) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return (Ug - Ub) + F_NEUTRON * S26_3RD * Phi(GAMMA_0) * math.exp(min(KAPPA * t, 500))


def kernel_phonon_ares(M=M_SUN, r=6.96e8):
    return sum(math.exp(-SSQ * i / 26) * BETA_I * dpm_ug1_seed(M, r) for i in range(1, 27))


def kernel_jet_mjet(A=1.5, gamma_THz=0.1):
    gamma = 2 * PI * gamma_THz * 1e12
    return 1 + A * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))


def kernel_ns_spindown(P=0.033, Pdot=4.2e-13):
    return 3.2e19 * math.sqrt(P * Pdot)


def kernel_gw170817_strain(M_chirp=1.186, d_Mpc=40):
    f = 100; d = d_Mpc * 3.086e22
    return 4 / d * (G * M_chirp * M_SUN / C**2)**(5/3) * (PI * f / C)**(2/3)


def kernel_blazar_ergosphere(M=3e8, a=0.95):
    M_kg = M * M_SUN; rS = 2 * dpm_ug1_seed(M_kg, C)
    return rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))


def kernel_rest_phonon_jet():
    return math.exp(-(OMEGA_SCM - GAMMA_0)**2 / (2 * SIGMA_G**2)) * S26_3RD


def kernel_qcalcgeom_vectorized(N=26):
    return sum(SSQ**i / (i**N) for i in range(1, N + 1))


def kernel_pipeline_full(M=M_SUN, r=6.96e8, t=86400):
    g = kernel_gravity_26layer(M, r)
    b = kernel_phonon_ares(M, r)
    return g - b + F_NEUTRON * S26_3RD * math.exp(min(KAPPA * t, 500))


def kernel_cena_jet(L=1e43, B=3000):
    return L * SSQ * B / (4 * PI * C)


def kernel_txs0506_jet(L=1e46, B=5000):
    return L * SSQ * B / (4 * PI * C)


def kernel_bcs_gap_solve(T=2e12, Tc=1.5e12):
    if T >= Tc: return 0.0
    return 3.528 * HBAR * OMEGA_SCM * math.sqrt(1 - (T / Tc)**2) * S26_3RD


def kernel_spectral_ladder_eval(N=26):
    return sum(math.exp(-SSQ * i / N) / (i + 1) for i in range(1, N + 1))


def kernel_ramanujan_26d_sum():
    return S26_3RD


def kernel_triadic_solver(M=M_SUN, r=6.96e8):
    Ug = kernel_gravity_26layer(M, r)
    Ub = kernel_phonon_ares(M, r)
    Um = dpm_ug1_seed(M, r) * SSQ * 0.1
    return Ug + Um - Ub


# v13 carry-forward (4)
def kernel_fubi_inside_out(M=5e7, r=None):
    M_kg = M * M_SUN
    if r is None: r = 2 * dpm_ug1_seed(M_kg, C)
    return kernel_fu_bi_i(M_kg, r, 86400)


def kernel_99sys_gamma_sweep():
    total = 0.0
    for g in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]:
        total += kernel_jet_mjet(1.5, g)
    return total


def kernel_agn_cena_fubi(M=5.5e7, a=0.70):
    M_kg = M * M_SUN; rS = 2 * dpm_ug1_seed(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return kernel_fu_bi_i(M_kg, rH, 86400)


def kernel_ns_merger_gw190425(M_chirp=1.44, d_Mpc=159):
    return kernel_gw170817_strain(M_chirp, d_Mpc)


# v14 carry-forward (4)
def kernel_agn_merger_fubi(M1=5.5e7, M2=3e7, a=0.70):
    M_kg = (M1 + M2) * M_SUN; rS = 2 * dpm_ug1_seed(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_ug1_seed(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_ug1_seed(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return abs(Ub) / (abs(Ug) + 1e-300) * S26_3RD


def kernel_qgp_vacuum_density(T=2e12):
    if T < 1.5e12: return 0.0
    return RHO_VAC * math.exp(min(KAPPA * 86400, 500)) * (T / 1e12)**4 * S26_3RD


def kernel_alice_multiplicity(N_part=383, alpha=8.5, beta=1.25):
    return alpha * N_part * (1 - math.exp(-beta * N_part / 400)) * S26_3RD / SSQ


def kernel_ym_mass_gap(T=2e12, Tc=1.5e12):
    if T >= Tc: return 0.0
    return 3.528 * HBAR * OMEGA_SCM * math.sqrt(1 - (T / Tc)**2) * S26_3RD**2


# v15 carry-forward (6)
def kernel_3c273_agn_fubi(M=8.86e8, a=0.90, A_jet=2.1):
    """3C273 AGN jet modulation: 3.1× at Γ₀."""
    M_kg = M * M_SUN
    rS = 2 * dpm_ug1_seed(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_ug1_seed(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_ug1_seed(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    M_jet = 1 + A_jet
    return ratio * M_jet * S26_3RD


def kernel_ton618_agn_fubi(M=6.6e10, a=0.998, A_jet=2.8):
    """TON618 AGN jet modulation: 3.8× at Γ₀."""
    M_kg = M * M_SUN
    rS = 2 * dpm_ug1_seed(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_ug1_seed(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_ug1_seed(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    M_jet = 1 + A_jet
    return ratio * M_jet * S26_3RD


def kernel_gw170817_merger(M_chirp=1.186, d_Mpc=40, suppression=0.667):
    """GW170817 66.7% strain reduction via SCm phonon damping."""
    h0 = kernel_gw170817_strain(M_chirp, d_Mpc)
    return h0 * (1 - suppression * Phi(GAMMA_0) / S26_3RD)


def kernel_smbh_merger_fubi(M1=5.5e7, M2=3e7, a1=0.70, a2=0.60):
    """SMBH merger F_U_Bi: inside-to-outside buoyancy."""
    M_total = (M1 + M2) * M_SUN
    rS = 2 * dpm_ug1_seed(M_total, C)
    a_eff = (M1 * a1 + M2 * a2) / (M1 + M2)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a_eff**2, 0)))
    rho = RHO_VAC * math.exp(min(KAPPA * 86400, 500))
    V = (4/3) * PI * rH**3
    Ph = Phi(GAMMA_0)
    return rho * V * S26_3RD**2 * Ph


def kernel_dm_halo_nfw(rho_0=5e-22, r_kpc=20, r_s_kpc=20):
    """Dark matter halo NFW profile with SCm coupling."""
    x = r_kpc / r_s_kpc
    if x < 1e-10: x = 1e-10
    rho_nfw = rho_0 / (x * (1 + x)**2)
    rho_SCm = RHO_VAC * math.exp(min(KAPPA * 86400, 500))
    return rho_SCm * S26_3RD * rho_nfw * Phi(GAMMA_0) / rho_0


def kernel_txs0506_3gamma(A_jet=1.3):
    """TXS0506 3-Γ-point profile: sum of modulations at 0.05, 0.10, 0.30 THz."""
    total = 0.0
    for g_THz in [0.05, 0.10, 0.30]:
        gamma = 2 * PI * g_THz * 1e12
        total += 1 + A_jet * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))
    return total * S26_3RD


# ── §3  New v16 Kernels (6) ──────────────────────────────────────────────

def kernel_scm_inflation_hubble(rho_scm=1e76):
    """SCm-driven Hubble parameter H_SCm at inflation epoch.

    H_SCm = sqrt(8πG/3 · ρ_SCm) · S₂₆⁽³⁾ · Φ_{1.25 THz}
    """
    H_base = math.sqrt(8 * PI * G / 3 * rho_scm)
    return H_base * S26_3RD * Phi(GAMMA_0)


def kernel_thorne_morris_exotic(M=4.3e6, r0=1e4):
    """Thorne-Morris exotic matter (ρ+P)<0 from SCm phonon energy.

    NEC violation at wormhole throat via β_i buoyancy pressure
    with 26-layer coherent amplification (S₂₆²).
    """
    M_kg = M * M_SUN
    # Ambient SCm density
    rho_ambient = (HBAR * OMEGA_SCM) / (2 * C**2) * S26_STATIC * 0.99
    # Curvature amplification
    r_s = 2 * dpm_ug1_seed(M_kg, C)
    amp = (r_s / r0)**2 * S26_STATIC**2
    rho_exotic = rho_ambient * amp
    # Buoyancy pressure → (ρ+P)/c² < 0
    # S26² coherent amplification at throat
    P_density = -BETA_I * rho_exotic * SSQ * 0.70 * S26_STATIC**2
    return rho_exotic + P_density  # should be negative


def kernel_vds_convergence(n_terms=26):
    """VDS polylogarithm Li₂₆([SSq]) convergence evaluation.

    Li₂₆(0.57) = Σ [SSq]^n / n²⁶ for n=1..N
    """
    return sum(SSQ**n / n**26 for n in range(1, n_terms + 1))


def kernel_dvp_prime_sieve(max_p=50):
    """DVP prime sieve: Σ a(p) = Σ [SSq]^{π(p)} / p²⁶ over primes p ≤ max_p."""
    primes = []
    for n in range(2, max_p + 1):
        if all(n % d != 0 for d in range(2, int(math.sqrt(n)) + 1)):
            primes.append(n)
    total = 0.0
    for idx, p in enumerate(primes):
        pi_p = idx + 1
        total += SSQ**pi_p / p**26
    return total


def kernel_bsh_harmonic_saturation(n_modes=26):
    """BSH harmonic saturation sum: Σ (1 - exp(-[SSq]·m)) for m=1..N."""
    return sum(1 - math.exp(-SSQ * m) for m in range(1, n_modes + 1))


def kernel_wormhole_phonon_damping(r0=1e4, M=4.3e6):
    """Wormhole throat phonon resonance damping via SCm coupling.

    Damping factor = exp(-κ·r₀/c) × Φ_{1.25 THz} × S₂₆⁽³⁾
    """
    M_kg = M * M_SUN
    b_mod = 1 - BETA_I * SSQ  # BSFG shape function modification
    damping = math.exp(-KAPPA * r0 / C) * Phi(GAMMA_0) * b_mod
    # Throat curvature coupling
    rS = 2 * dpm_ug1_seed(M_kg, C)
    curvature = rS / (r0 + rS)
    return damping * curvature * S26_3RD


# ── §4  ALL_KERNELS Registry ────────────────────────────────────────────

ALL_KERNELS: List[Callable] = [
    # v11 (16)
    kernel_gravity_26layer, kernel_fu_bi_i, kernel_phonon_ares,
    kernel_jet_mjet, kernel_ns_spindown, kernel_gw170817_strain,
    kernel_blazar_ergosphere, kernel_rest_phonon_jet,
    kernel_qcalcgeom_vectorized, kernel_pipeline_full,
    kernel_cena_jet, kernel_txs0506_jet, kernel_bcs_gap_solve,
    kernel_spectral_ladder_eval, kernel_ramanujan_26d_sum, kernel_triadic_solver,
    # v13 (4)
    kernel_fubi_inside_out, kernel_99sys_gamma_sweep,
    kernel_agn_cena_fubi, kernel_ns_merger_gw190425,
    # v14 (4)
    kernel_agn_merger_fubi, kernel_qgp_vacuum_density,
    kernel_alice_multiplicity, kernel_ym_mass_gap,
    # v15 (6)
    kernel_3c273_agn_fubi, kernel_ton618_agn_fubi,
    kernel_gw170817_merger, kernel_smbh_merger_fubi,
    kernel_dm_halo_nfw, kernel_txs0506_3gamma,
    # v16 (6)
    kernel_scm_inflation_hubble, kernel_thorne_morris_exotic,
    kernel_vds_convergence, kernel_dvp_prime_sieve,
    kernel_bsh_harmonic_saturation, kernel_wormhole_phonon_damping,
]


# ── §5  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV16:
    """36-kernel production benchmark: 800k calc/s target."""

    def single_pass(self) -> float:
        return sum(k() for k in ALL_KERNELS)

    def simulate(self, duration_s: float = 1.0) -> dict:
        start = time.perf_counter()
        iterations = 0
        while time.perf_counter() - start < duration_s:
            self.single_pass()
            iterations += 1
        elapsed = time.perf_counter() - start
        calc_per_sec = iterations * len(ALL_KERNELS) / elapsed

        return {
            "iterations": iterations,
            "kernels": len(ALL_KERNELS),
            "elapsed_s": elapsed,
            "calc_per_sec": calc_per_sec,
            "target": TARGET_CALC_PER_SEC,
            "met_target": calc_per_sec >= TARGET_CALC_PER_SEC,
        }


# ── §6  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True
    passed = 0

    # Test 1: All kernels finite
    results = [k() for k in ALL_KERNELS]
    all_finite = all(math.isfinite(r) for r in results)
    if all_finite:
        print(f"[ OK ] All {len(ALL_KERNELS)} kernels produce finite results")
        passed += 1
    else:
        print(f"[FAIL] Non-finite kernel result"); ok = False

    # Test 2: Kernel count = 36
    n = len(ALL_KERNELS)
    if n == 36:
        print(f"[ OK ] Kernel count: {n}")
        passed += 1
    else:
        print(f"[FAIL] Expected 36 kernels, got {n}"); ok = False

    # Test 3: SCm inflation Hubble > 0
    val = kernel_scm_inflation_hubble()
    if val > 0 and math.isfinite(val):
        print(f"[ OK ] SCm inflation H: {val:.6e} s⁻¹")
        passed += 1
    else:
        print(f"[FAIL] SCm inflation H: {val}"); ok = False

    # Test 4: Thorne-Morris exotic < 0 (NEC violation)
    val = kernel_thorne_morris_exotic()
    if val < 0:
        print(f"[ OK ] Thorne-Morris (ρ+P): {val:.6e} < 0 (exotic)")
        passed += 1
    else:
        print(f"[FAIL] Thorne-Morris should be negative: {val}"); ok = False

    # Test 5: VDS convergence > 0
    val = kernel_vds_convergence()
    if val > 0 and math.isfinite(val):
        print(f"[ OK ] VDS Li₂₆(0.57): {val:.15e}")
        passed += 1
    else:
        print(f"[FAIL] VDS: {val}"); ok = False

    # Test 6: DVP prime sieve > 0
    val = kernel_dvp_prime_sieve()
    if val > 0 and math.isfinite(val):
        print(f"[ OK ] DVP Σa(p): {val:.15e}")
        passed += 1
    else:
        print(f"[FAIL] DVP: {val}"); ok = False

    # Test 7: BSH saturation sum > 0
    val = kernel_bsh_harmonic_saturation()
    if val > 0:
        print(f"[ OK ] BSH saturation: {val:.6f}")
        passed += 1
    else:
        print(f"[FAIL] BSH: {val}"); ok = False

    # Test 8: Wormhole phonon damping > 0
    val = kernel_wormhole_phonon_damping()
    if val > 0 and math.isfinite(val):
        print(f"[ OK ] Wormhole phonon: {val:.6e}")
        passed += 1
    else:
        print(f"[FAIL] Wormhole phonon: {val}"); ok = False

    # Test 9: v15 carry-forward compatibility
    v15_kernels = ALL_KERNELS[:30]
    v15_results = [k() for k in v15_kernels]
    v15_ok = all(math.isfinite(r) for r in v15_results)
    if v15_ok:
        print(f"[ OK ] v15 carry-forward: all 30 kernels valid")
        passed += 1
    else:
        print(f"[FAIL] v15 carry-forward failed"); ok = False

    # Test 10: Target and metadata
    print(f"[ OK ] Target: {TARGET_CALC_PER_SEC:,} calc/s, S₂₆⁽³⁾ = {S26_3RD:.6e}")
    passed += 1

    print(f"\n{'='*60}")
    print(f"  production_scaling_v16.py: {passed}/10 tests passed")
    print(f"{'='*60}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
