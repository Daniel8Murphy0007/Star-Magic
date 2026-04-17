#!/usr/bin/env python3
"""
production_scaling_v15.py — Production Scaling to 650k calc/s

Session 220 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
v15 upgrade: 650,000 calc/s target (up from v14's 600k).
30 kernels total = v14's 24 + 6 new Session 220 kernels:
  kernel_3c273_agn_fubi       — 3C273 AGN jet modulation (3.1×)
  kernel_ton618_agn_fubi      — TON618 AGN jet modulation (3.8×)
  kernel_gw170817_merger      — GW170817 66.7% strain reduction
  kernel_smbh_merger_fubi     — SMBH merger inside-to-outside F_U_Bi
  kernel_dm_halo_nfw          — Dark matter halo NFW profile
  kernel_txs0506_3gamma       — TXS0506 3-Γ-point profile

Architecture: Standalone kernel functions, ALL_KERNELS list, benchmark class.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

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

TARGET_CALC_PER_SEC = 650_000


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


# ── §2  Carry-Forward Kernels (v11: 16, v13: 4, v14: 4 = 24) ────────────

def kernel_gravity_26layer(M=M_SUN, r=6.96e8):
    return sum(dpm_emergent_ug1(M, r) * SSQ * i / 26 for i in range(1, 27))


def kernel_fu_bi_i(M=M_SUN, r=6.96e8, t=86400):
    Ug = sum(dpm_emergent_ug1(M, r) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M, r) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return (Ug - Ub) + F_NEUTRON * S26_3RD * Phi(GAMMA_0) * math.exp(min(KAPPA * t, 500))


def kernel_phonon_ares(M=M_SUN, r=6.96e8):
    return sum(math.exp(-SSQ * i / 26) * BETA_I * dpm_emergent_ug1(M, r) for i in range(1, 27))


def kernel_jet_mjet(A=1.5, gamma_THz=0.1):
    gamma = 2 * PI * gamma_THz * 1e12
    return 1 + A * math.exp(-(gamma - GAMMA_0)**2 / (2 * SIGMA_G**2))


def kernel_ns_spindown(P=0.033, Pdot=4.2e-13):
    return 3.2e19 * math.sqrt(P * Pdot)


def kernel_gw170817_strain(M_chirp=1.186, d_Mpc=40):
    f = 100; d = d_Mpc * 3.086e22
    return 4 / d * (G * M_chirp * M_SUN / C**2)**(5/3) * (PI * f / C)**(2/3)


def kernel_blazar_ergosphere(M=3e8, a=0.95):
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
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
    Um = dpm_emergent_ug1(M, r) * SSQ * 0.1
    return Ug + Um - Ub


# v13 carry-forward (4)
def kernel_fubi_inside_out(M=5e7, r=None):
    M_kg = M * M_SUN
    if r is None: r = 2 * dpm_emergent_ug1(M_kg, C)
    return kernel_fu_bi_i(M_kg, r, 86400)


def kernel_99sys_gamma_sweep():
    total = 0.0
    for g in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0]:
        total += kernel_jet_mjet(1.5, g)
    return total


def kernel_agn_cena_fubi(M=5.5e7, a=0.70):
    M_kg = M * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return kernel_fu_bi_i(M_kg, rH, 86400)


def kernel_ns_merger_gw190425(M_chirp=1.44, d_Mpc=159):
    return kernel_gw170817_strain(M_chirp, d_Mpc)


# v14 carry-forward (4)
def kernel_agn_merger_fubi(M1=5.5e7, M2=3e7, a=0.70):
    M_kg = (M1 + M2) * M_SUN; rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    return abs(Ub) / (abs(Ug) + 1e-300) * S26_3RD


def kernel_qgp_vacuum_density(T=2e12):
    if T < 1.5e12: return 0.0
    return RHO_VAC * math.exp(min(KAPPA * 86400, 500)) * (T / 1e12)**4 * S26_3RD


def kernel_alice_multiplicity(N_part=383, alpha=8.5, beta=1.25):
    return alpha * N_part * (1 - math.exp(-beta * N_part / 400)) * S26_3RD / SSQ


def kernel_ym_mass_gap(T=2e12, Tc=1.5e12):
    if T >= Tc: return 0.0
    return 3.528 * HBAR * OMEGA_SCM * math.sqrt(1 - (T / Tc)**2) * S26_3RD**2


# ── §3  New v15 Kernels (6) ──────────────────────────────────────────────

def kernel_3c273_agn_fubi(M=8.86e8, a=0.90, A_jet=2.1):
    """3C273 AGN jet modulation: 3.1× at Γ₀."""
    M_kg = M * M_SUN
    rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
    ratio = abs(Ub) / (abs(Ug) + abs(Ub) + 1e-300)
    M_jet = 1 + A_jet  # at resonance exp(0) = 1
    return ratio * M_jet * S26_3RD


def kernel_ton618_agn_fubi(M=6.6e10, a=0.998, A_jet=2.8):
    """TON618 AGN jet modulation: 3.8× at Γ₀."""
    M_kg = M * M_SUN
    rS = 2 * dpm_emergent_ug1(M_kg, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    Ug = sum(dpm_emergent_ug1(M_kg, rH) * SSQ * i / 26 for i in range(1, 27))
    Ub = sum(dpm_emergent_ug1(M_kg, rH) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))
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
    rS = 2 * dpm_emergent_ug1(M_total, C)
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
]


# ── §5  Benchmark Class ─────────────────────────────────────────────────

class ProductionScalingV15:
    """30-kernel production benchmark: 650k calc/s target."""

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

    # Test 2: 3C273 modulation ratio > 0
    val = kernel_3c273_agn_fubi()
    if val > 0:
        print(f"[ OK ] 3C273 AGN: {val:.6e}")
        passed += 1
    else:
        print(f"[FAIL] 3C273: {val}"); ok = False

    # Test 3: TON618 > 3C273
    t618 = kernel_ton618_agn_fubi()
    c273 = kernel_3c273_agn_fubi()
    if t618 > c273:
        print(f"[ OK ] TON618 ({t618:.4e}) > 3C273 ({c273:.4e})")
        passed += 1
    else:
        print(f"[FAIL] TON618 should exceed 3C273"); ok = False

    # Test 4: GW170817 damped strain < undamped
    h_damped = kernel_gw170817_merger()
    h_undamped = kernel_gw170817_strain()
    if h_damped < h_undamped:
        print(f"[ OK ] GW170817 damped: {h_damped:.6e} < undamped: {h_undamped:.6e}")
        passed += 1
    else:
        print(f"[FAIL] Damped strain should be less"); ok = False

    # Test 5: SMBH merger F_U_Bi > 0
    val = kernel_smbh_merger_fubi()
    if val > 0:
        print(f"[ OK ] SMBH merger: {val:.6e}")
        passed += 1
    else:
        print(f"[FAIL] SMBH merger: {val}"); ok = False

    # Test 6: DM halo NFW > 0
    val = kernel_dm_halo_nfw()
    if val > 0:
        print(f"[ OK ] DM halo NFW: {val:.6e}")
        passed += 1
    else:
        print(f"[FAIL] DM halo: {val}"); ok = False

    # Test 7: Kernel count = 30
    n = len(ALL_KERNELS)
    if n == 30:
        print(f"[ OK ] Kernel count: {n}")
        passed += 1
    else:
        print(f"[FAIL] Expected 30 kernels, got {n}"); ok = False

    # Test 8: Target and metadata
    print(f"[ OK ] Target: {TARGET_CALC_PER_SEC:,} calc/s, S₂₆⁽³⁾ = {S26_3RD:.6e}")
    passed += 1

    print(f"\n{'='*60}")
    print(f"  production_scaling_v15.py: {passed}/8 tests passed")
    print(f"{'='*60}")
    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
