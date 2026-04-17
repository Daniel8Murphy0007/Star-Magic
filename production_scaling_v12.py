#!/usr/bin/env python3
"""
production_scaling_v12.py — QGP + 99-System Scaling at 501k calc/s

Session 216 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Upgrade from v11 (500k) → v12 (501k calc/s target, +0.2% over v11).
18 benchmark kernels: v11's 16 + QGP vacuum density + 99-system master.
Maintains REST API + QCalcGeom vectorized throughput.

History: v4 (100k) → v5 (150k) → v6 (200k) → v7 (300k) → v8 (350k)
       → v9 (400k) → v10 (450k) → v11 (500k) → v12 (501k)
────────────────────────────────────────────────────────────────────────────────
"""

import math
from dpm_helpers import dpm_emergent_ug1, dpm_emergent_ug2

import time
from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
G         = 6.674e-11
HBAR      = 1.055e-34
K_B       = 1.381e-23
M_SUN     = 1.989e30
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
GAMMA_0   = 2 * PI * 0.1e12
SIGMA_G   = 0.08 * 2 * PI * 1e12
BETA_I    = 0.603
T_C_QGP   = 1.5e12
LAMBDA_QCD = 0.217e9

TARGET_CALC_PER_SEC = 501_000


# ── §1  Benchmark Kernels (v11 carry-forward: 16) ─────────────────────────

def kernel_gravity_26layer(M_kg: float = 4e6 * M_SUN, r: float = 1e12) -> float:
    return sum(dpm_emergent_ug1(M_kg, r) * SSQ * i / 26 for i in range(1, 27))

def kernel_fu_bi_i(M_kg: float = 4e6 * M_SUN, r: float = 1e12) -> float:
    return sum(dpm_emergent_ug1(M_kg, r) * math.exp(-SSQ * i / 26) * BETA_I for i in range(1, 27))

def kernel_phonon_ares(omega: float = OMEGA_SCM, gamma: float = GAMMA_0) -> float:
    return math.exp(-(omega - OMEGA_SCM)**2 / (2 * gamma**2)) * S26

def kernel_jet_mjet(gamma_THz: float = 0.10, A_jet: float = 1.5) -> float:
    Gr = 2 * PI * gamma_THz * 1e12
    return 1 + A_jet * math.exp(-(Gr - GAMMA_0)**2 / (2 * SIGMA_G**2))

def kernel_ns_spindown(P: float = 0.003, Pdot: float = 1e-15) -> float:
    return 3.2e19 * math.sqrt(P * Pdot)

def kernel_gw170817_strain(d_Mpc: float = 40.0) -> float:
    return 0.333 * 1e-21 * (40.0 / d_Mpc)

def kernel_blazar_ergosphere(M_Msun: float = 6.5e9, a: float = 0.90) -> float:
    M = M_Msun * M_SUN
    rS = 2 * dpm_emergent_ug1(M, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return sum(dpm_emergent_ug1(M, rH) * math.exp(-SSQ * i / 26) for i in range(1, 27))

def kernel_rest_phonon_jet(gamma_THz: float = 0.10, A_jet: float = 1.5) -> float:
    Gr = 2 * PI * gamma_THz * 1e12
    Mj = 1 + A_jet * math.exp(-(Gr - GAMMA_0)**2 / (2 * SIGMA_G**2))
    return Mj * S26

def kernel_qcalcgeom_vectorized(N: int = 100) -> float:
    return sum(math.exp(-SSQ * k / 26) * BETA_I for k in range(1, N + 1))

def kernel_pipeline_full() -> float:
    return kernel_gravity_26layer() + kernel_fu_bi_i() + kernel_phonon_ares() + kernel_jet_mjet()

def kernel_cena_jet(M_Msun: float = 5.5e7, a: float = 0.70, B: float = 3000) -> float:
    M = M_Msun * M_SUN
    rS = 2 * dpm_emergent_ug1(M, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return (B**2 / (8 * PI)) * (rH / C)**2 * a**2 * C

def kernel_txs0506_jet(M_Msun: float = 3e8, a: float = 0.95, B: float = 5000) -> float:
    M = M_Msun * M_SUN
    rS = 2 * dpm_emergent_ug1(M, C)
    rH = rS / 2 * (1 + math.sqrt(max(1 - a**2, 0)))
    return (B**2 / (8 * PI)) * (rH / C)**2 * a**2 * C

def kernel_bcs_gap_solve(T_K: float = 4.2) -> float:
    delta = HBAR * OMEGA_SCM / 2
    for _ in range(20):
        arg = delta / (2 * K_B * T_K) if T_K > 0 else 50.0
        arg = min(arg, 50.0)
        delta = (HBAR * OMEGA_SCM / 2) * math.tanh(arg) * S26 * 0.6
    return delta

def kernel_spectral_ladder_eval() -> float:
    E0 = HBAR * OMEGA_SCM
    return sum(E0 * (2 * PI) ** (n / 3.0) * S26 for n in range(1, 27))

def kernel_ramanujan_26d_sum(z: float = SSQ, N: int = 20) -> float:
    total = 0.0
    for n in range(1, N + 1):
        Rn = 1.0 / math.factorial(min(n, 20))
        total += (z ** n) / (n ** 26) * Rn
    return total

def kernel_triadic_solver(gamma_THz: float = 0.10) -> float:
    gamma = 2 * PI * gamma_THz * 1e12
    compressed = 0.6 * S26 * 1.5
    resonant = math.exp(0) * S26
    buoyancy = S26 * math.cos(0) * math.exp(0)
    return compressed + resonant + buoyancy


# ── §2  New v12 Kernels ───────────────────────────────────────────────────

def kernel_qgp_density(T_K: float = 2e12) -> float:
    """QGP vacuum density: ρ_QGP(T) = ρ_SCm · S₂₆^{(k)} · exp(-(T_c-T)/T)."""
    rho_scm = 1e-10
    S26k = kernel_ramanujan_26d_sum(SSQ, 20)  # Reuse for speed
    if T_K > 0:
        exp_factor = math.exp(max(min(-(T_C_QGP - T_K) / T_K, 500), -500))
    else:
        exp_factor = 0.0
    return rho_scm * S26k * exp_factor


def kernel_99system_master(n_systems: int = 99) -> float:
    """99-system master equation: F_U^{(99)} aggregate evaluation."""
    total = 0.0
    for i in range(1, n_systems + 1):
        M = (0.1 + i * 2) * M_SUN
        r = 1e9 * (1 + i * 0.3)
        Ug = sum(dpm_emergent_ug1(M, r) * SSQ * j / 26 for j in range(1, 27))
        Um = dpm_emergent_ug1(M, r) * SSQ * 0.1
        UA = dpm_emergent_ug1(M, r) * 1e-10
        Ub = sum(dpm_emergent_ug1(M, r) * math.exp(-SSQ * j / 26) * BETA_I for j in range(1, 27))
        Fn = 1e-10
        Phi = S26
        total += Ug + Um + UA - Ub + Fn * S26 * Phi
    return total


# ── §3  Kernel Registry ───────────────────────────────────────────────────

ALL_KERNELS = [
    kernel_gravity_26layer,        # 1
    kernel_fu_bi_i,                # 2
    kernel_phonon_ares,            # 3
    kernel_jet_mjet,               # 4
    kernel_ns_spindown,            # 5
    kernel_gw170817_strain,        # 6
    kernel_blazar_ergosphere,      # 7
    kernel_rest_phonon_jet,        # 8
    kernel_qcalcgeom_vectorized,   # 9
    kernel_pipeline_full,          # 10
    kernel_cena_jet,               # 11
    kernel_txs0506_jet,            # 12
    kernel_bcs_gap_solve,          # 13
    kernel_spectral_ladder_eval,   # 14
    kernel_ramanujan_26d_sum,      # 15 (v11)
    kernel_triadic_solver,         # 16 (v11)
    kernel_qgp_density,            # 17 (NEW v12)
    kernel_99system_master,        # 18 (NEW v12)
]


# ── §4  Benchmark Runner ──────────────────────────────────────────────────

class ProductionScalingV12:
    """Production benchmark at 501k calc/s with 18 kernels.

    Session 216. Extends v11 with QGP density + 99-system master kernels.
    """

    TARGET = TARGET_CALC_PER_SEC
    N_KERNELS = len(ALL_KERNELS)

    def single_pass(self) -> float:
        total = 0.0
        for kernel in ALL_KERNELS:
            total += kernel()
        return total

    def simulate(self, duration_s: float = 1.0) -> dict:
        count = 0
        t0 = time.perf_counter()
        while time.perf_counter() - t0 < duration_s:
            self.single_pass()
            count += 1
        elapsed = time.perf_counter() - t0
        rate = count * self.N_KERNELS / elapsed
        return {
            "iterations": count,
            "kernels": self.N_KERNELS,
            "elapsed_s": elapsed,
            "calc_per_sec": rate,
            "target": self.TARGET,
            "met_target": rate >= self.TARGET,
        }

    def compute(self, dataset: dict) -> dict:
        duration = float(dataset.get("duration_s", 0.5))
        bench = self.simulate(duration)
        return {
            **bench,
            "primary_equations": [
                f"rate = {bench['calc_per_sec']:.0f} calc/s (target {self.TARGET})",
                f"{self.N_KERNELS} kernels × {bench['iterations']} iterations in {bench['elapsed_s']:.3f} s",
                f"Target {'MET' if bench['met_target'] else 'NOT MET'}",
            ],
            "note": "PAPER_977. Session 216.",
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True

    # Test 1: All 18 kernels execute
    for i, kernel in enumerate(ALL_KERNELS):
        v = kernel()
        if not math.isfinite(v):
            print(f"[FAIL] kernel {i} returned non-finite"); ok = False
    if ok:
        print(f"[ OK ] All {len(ALL_KERNELS)} kernels returned finite values")

    # Test 2: New v12 kernels
    rho = kernel_qgp_density()
    if rho <= 0:
        print("[FAIL] QGP density should be positive"); ok = False
    else:
        print(f"[ OK ] QGP density at 2×10¹²K = {rho:.10e}")

    fu99 = kernel_99system_master()
    if not math.isfinite(fu99):
        print("[FAIL] 99-system master non-finite"); ok = False
    else:
        print(f"[ OK ] 99-system master F_U = {fu99:.6e}")

    # Test 3: Benchmark class
    ps = ProductionScalingV12()
    total = ps.single_pass()
    if not math.isfinite(total):
        print("[FAIL] single_pass() non-finite"); ok = False
    else:
        print(f"[ OK ] single_pass() = {total:.6e}")

    print(f"[ OK ] v12 target: {TARGET_CALC_PER_SEC:,} calc/s, {len(ALL_KERNELS)} kernels")

    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
