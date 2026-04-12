#!/usr/bin/env python3
"""
Production Scaling v8 — 350k calc/s Benchmark Harness with REST API + QCalcGeom Vectorized Kernels

Session 212 | PAPER_EXP_BENCH_V8 | Daniel Murphy
PURPOSE: Benchmark harness measuring calculations/second for the extended UQFF
         computation pipeline with REST API throughput and QCalcGeom vectorized
         phonon batch kernels. Target: 350k calc/s (3.5× v4 target, +17% over v7).

ARCHITECTURE:
  production_scaling_v7.py (300k baseline)  ──┐
  blazar_jet_phonon.py (blazar jets)          ├──→ this module ──→ benchmark report
  ns_phonon_gw170817_wstp.py (GW170817)      │
  uqff_server.js /api/phonon/jet (REST)  ────┘

BENCHMARKS:
  1. 26-Layer Gravity (baseline, from v4/v7)
  2. F_U_Bi_i Assembly (baseline, from v4/v7)
  3. Phonon Resonance a_res (from v7)
  4. Jet Modulation M_jet(Γ) sweep (from v7)
  5. NS Spindown Ω̇ phonon correction (from v7)
  6. GW170817 h_UQFF(t) strain (new, replacing GW190425)
  7. Blazar Ergosphere Phonon Coupling (new)
  8. REST API /api/phonon/jet simulated throughput (new)
  9. QCalcGeom Vectorized Phonon Batch (new, numpy)
  10. Full Pipeline v8 (all kernels combined)

TARGET: 350,000 calc/s (3.5× v4's 100k)
"""

import math
import time
import json
import os
import sys
import statistics
from typing import Dict, List, Tuple, Optional, Any, Callable
from dataclasses import dataclass, field
from datetime import datetime

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

# ── §1  CONSTANTS ───────────────────────────────────────────────────────────

G       = 6.67430e-11
c       = 2.99792e8
hbar    = 1.05457e-34
PI      = math.pi
M_sun   = 1.98892e30
k_B     = 1.38065e-23

KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_UA  = 7.09e-36
RHO_SCM = 7.09e-37
ETA_AETHER = 1e-22
E_REACT = 1e46
N_LEVELS = 26

OMEGA_SCM = 2 * PI * 1.25e12
GAMMA_SCM = 2 * PI * 0.1e12
PHI_0     = 1e20

A_JET_DEFAULT = 1.5
SIGMA_G_RAD = 0.08 * 2 * PI * 1e12

TARGET_CALC_PER_SEC = 350_000


# ── §2  BENCHMARK KERNELS (v7 baseline) ────────────────────────────────────

def kernel_26layer_gravity(M: float, r: float, mu_s: float,
                           omega_s: float, t: float) -> Dict[str, float]:
    """26-layer compressed gravity (from v4/v7)."""
    Ug_total = 0.0
    Ubi_total = 0.0
    for layer in range(1, 27):
        qi = SSQ * layer / N_LEVELS
        Ug_i = G * M / (r**2) * qi * (1 + KAPPA * t)
        Ubi_i = BETA_I * Ug_i * math.exp(-U_UA * layer)
        Ug_total += Ug_i
        Ubi_total += Ubi_i
    return {"Ug_total": Ug_total, "Ubi_total": Ubi_total}


def kernel_fu_bi_i(M: float, r: float, mu_s: float,
                   omega_s: float, t: float) -> Dict[str, float]:
    """F_U_Bi_i assembly (from v4/v7)."""
    Ug1 = mu_0 * mu_s / (4 * PI * r**3)
    Ug2 = E_REACT * ETA_AETHER / (r**2)
    Ug3 = mu_s * omega_s * math.sin(omega_s * t) / (r**2)
    Ug4 = RHO_SCM * G * M / (r * RHO_UA)
    Ug_sum = Ug1 + Ug2 + Ug3 + Ug4
    F_U = Ug_sum * (1 + KAPPA * t) * H_SCM
    F_U_Bi = BETA_I * F_U
    return {"F_U": F_U, "F_U_Bi": F_U_Bi}


def kernel_phonon_ares(omega: float, Gamma: float) -> float:
    """Phonon resonance acceleration a_res (from v7)."""
    delta = omega - OMEGA_SCM
    gaussian = math.exp(-delta**2 / (2 * Gamma**2))
    S26 = 0.57**26 / math.factorial(26) if SSQ < 1 else 1e-20
    return gaussian * S26 * PHI_0


def kernel_jet_modulation(Gamma: float) -> float:
    """Jet modulation M_jet(Γ) (from v7)."""
    delta = Gamma - GAMMA_SCM
    return 1 + A_JET_DEFAULT * math.exp(-delta**2 / (2 * SIGMA_G_RAD**2))


def kernel_ns_spindown(Omega_dot: float, Phi: float, S26: float) -> float:
    """NS spindown phonon correction (from v7)."""
    correction = Phi * S26 * SSQ / N_LEVELS
    return Omega_dot * (1 + correction)


# ── §3  NEW v8 KERNELS ─────────────────────────────────────────────────────

def kernel_gw170817_strain(t: float = 0.0) -> float:
    """GW170817 h_UQFF(t) = h_GR · 0.333 · exp([SSq]·t/26)."""
    h_GR = 5.4176e-22
    D_total = 0.333
    return h_GR * D_total * math.exp(SSQ * t / 26.0)


def kernel_blazar_ergosphere(M: float, a_spin: float,
                              Gamma_bulk: float, theta_obs: float) -> float:
    """Blazar ergosphere phonon coupling energy."""
    r_S = 2 * G * M / c**2
    r_H = r_S / 2 * (1 + math.sqrt(max(1 - a_spin**2, 0)))
    Omega_H = a_spin * c / (2 * r_H)
    beta = math.sqrt(1 - 1 / max(Gamma_bulk**2, 1.001))
    delta_D = 1 / (Gamma_bulk * (1 - beta * math.cos(theta_obs)))
    S26_val = SSQ**26 / math.factorial(26) if SSQ < 1 else 1e-20
    Phi_norm = S26_val  # at resonance gaussian=1
    return (a_spin / 2) * M * c**2 * Phi_norm * S26_val * delta_D


def kernel_rest_api_phonon_jet(M_bh: float, a_spin: float,
                                B: float, Gamma_THz: float) -> Dict[str, float]:
    """Simulated REST API /api/phonon/jet computation kernel."""
    Gamma_rad = Gamma_THz * 2 * PI * 1e12
    r_S = 2 * G * M_bh / c**2
    r_H = r_S / 2 * (1 + math.sqrt(max(1 - a_spin**2, 0)))
    P_BZ = (B**2 / (8 * PI)) * (r_H / c)**2 * a_spin**2 * c
    delta = Gamma_rad - GAMMA_SCM
    M_jet = 1 + A_JET_DEFAULT * math.exp(-delta**2 / (2 * SIGMA_G_RAD**2))
    P_jet = P_BZ * (1 + M_jet)
    return {"P_BZ": P_BZ, "P_jet": P_jet, "M_jet": M_jet, "enhancement": P_jet / max(P_BZ, 1e-50)}


def kernel_qcalcgeom_vectorized(n: int = 1000) -> float:
    """QCalcGeom vectorized phonon batch kernel (numpy)."""
    if not HAS_NUMPY:
        total = 0.0
        for i in range(n):
            omega = OMEGA_SCM + (i - n / 2) * 1e9
            delta = omega - OMEGA_SCM
            gaussian = math.exp(-delta**2 / (2 * GAMMA_SCM**2))
            total += gaussian * PHI_0 * SSQ
        return total / n

    omegas = np.linspace(OMEGA_SCM - 0.5e12, OMEGA_SCM + 0.5e12, n)
    deltas = omegas - OMEGA_SCM
    gaussians = np.exp(-deltas**2 / (2 * GAMMA_SCM**2))
    Phi_array = gaussians * PHI_0 * SSQ
    return float(np.mean(Phi_array))


# ── §4  BENCHMARK RUNNER ───────────────────────────────────────────────────

@dataclass
class BenchmarkResult:
    kernel_name: str
    rate_per_sec: float
    n_iterations: int
    elapsed_s: float
    meets_target: bool


def run_benchmark(name: str, func: Callable, n: int,
                  target: float = TARGET_CALC_PER_SEC) -> BenchmarkResult:
    """Run a single kernel benchmark."""
    t0 = time.perf_counter()
    for _ in range(n):
        func()
    elapsed = time.perf_counter() - t0
    rate = n / elapsed if elapsed > 0 else 0
    return BenchmarkResult(
        kernel_name=name,
        rate_per_sec=rate,
        n_iterations=n,
        elapsed_s=elapsed,
        meets_target=rate >= target,
    )


class ProductionScalingV8:
    """
    v8 production benchmark: 10 kernels targeting 350k calc/s.
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        n = dataset.get('n_iterations', 10000)
        target = dataset.get('target', TARGET_CALC_PER_SEC)

        M = 4e6 * M_sun
        r = 1e12
        mu_s = 1e25
        omega_s = 1.0
        t = 0.0

        benchmarks = []

        # 1. 26-layer gravity
        b = run_benchmark("26-Layer Gravity", lambda: kernel_26layer_gravity(M, r, mu_s, omega_s, t), n, target)
        benchmarks.append(b)

        # 2. F_U_Bi_i assembly
        b = run_benchmark("F_U_Bi_i Assembly", lambda: kernel_fu_bi_i(M, r, mu_s, omega_s, t), n, target)
        benchmarks.append(b)

        # 3. Phonon a_res
        b = run_benchmark("Phonon a_res", lambda: kernel_phonon_ares(OMEGA_SCM, GAMMA_SCM), n, target)
        benchmarks.append(b)

        # 4. Jet M_jet(Γ) sweep
        b = run_benchmark("Jet M_jet(Γ)", lambda: kernel_jet_modulation(GAMMA_SCM), n, target)
        benchmarks.append(b)

        # 5. NS spindown
        b = run_benchmark("NS Spindown", lambda: kernel_ns_spindown(-4.2e-15, PHI_0, 1e-20), n, target)
        benchmarks.append(b)

        # 6. GW170817 strain
        b = run_benchmark("GW170817 h_UQFF", lambda: kernel_gw170817_strain(0.0), n, target)
        benchmarks.append(b)

        # 7. Blazar ergosphere
        b = run_benchmark("Blazar Ergosphere", lambda: kernel_blazar_ergosphere(1e9 * M_sun, 0.95, 10.0, 0.1), n, target)
        benchmarks.append(b)

        # 8. REST API phonon/jet
        b = run_benchmark("REST /api/phonon/jet", lambda: kernel_rest_api_phonon_jet(6.5e9 * M_sun, 0.9, 50.0, 0.1), n, target)
        benchmarks.append(b)

        # 9. QCalcGeom vectorized
        b = run_benchmark("QCalcGeom Vectorized", lambda: kernel_qcalcgeom_vectorized(100), n // 10, target // 10)
        benchmarks.append(b)

        # 10. Full pipeline
        def full_pipeline():
            kernel_26layer_gravity(M, r, mu_s, omega_s, t)
            kernel_fu_bi_i(M, r, mu_s, omega_s, t)
            kernel_phonon_ares(OMEGA_SCM, GAMMA_SCM)
            kernel_jet_modulation(GAMMA_SCM)
            kernel_ns_spindown(-4.2e-15, PHI_0, 1e-20)
            kernel_gw170817_strain(0.0)
            kernel_blazar_ergosphere(1e9 * M_sun, 0.95, 10.0, 0.1)
            kernel_rest_api_phonon_jet(6.5e9 * M_sun, 0.9, 50.0, 0.1)

        b = run_benchmark("Full Pipeline v8", full_pipeline, n, target // 8)
        benchmarks.append(b)

        rates = [bm.rate_per_sec for bm in benchmarks[:8]]
        avg_rate = statistics.mean(rates) if rates else 0
        meets_target = avg_rate >= target

        return {
            "benchmarks": [
                {
                    "kernel": bm.kernel_name,
                    "rate_per_s": bm.rate_per_sec,
                    "n_iterations": bm.n_iterations,
                    "elapsed_s": bm.elapsed_s,
                    "meets_target": bm.meets_target,
                }
                for bm in benchmarks
            ],
            "avg_rate_per_s": avg_rate,
            "target_per_s": target,
            "meets_target": meets_target,
            "n_kernels": len(benchmarks),
            "primary_equations": [
                f"{bm.kernel_name}: {bm.rate_per_sec:.0f} calc/s" for bm in benchmarks
            ] + [
                f"Average (first 8): {avg_rate:.0f} calc/s",
                f"Target: {target} calc/s  {'MET' if meets_target else 'NOT MET'}",
            ],
            "note": "Production scaling v8 — 350k target. Session 212.",
        }

    def simulate(self, sweep=None, **kw):
        results = []
        for n in (sweep or [1000, 5000, 10000, 50000]):
            res = self.compute({"n_iterations": n})
            res["sweep_val"] = n
            results.append(res)
        return results


# ── §5  MAIN / DEMO ────────────────────────────────────────────────────────

def main():
    """Run v8 production scaling benchmark."""
    print("=" * 72)
    print(f"Production Scaling v8 — Target: {TARGET_CALC_PER_SEC:,} calc/s")
    print("=" * 72)

    v8 = ProductionScalingV8()
    result = v8.compute({"n_iterations": 10000})

    for bm in result["benchmarks"]:
        status = "✓" if bm["meets_target"] else "✗"
        print(f"  {status} {bm['kernel']:<25s}: {bm['rate_per_s']:>12,.0f} calc/s")

    print(f"\n  Average (8 core kernels): {result['avg_rate_per_s']:>12,.0f} calc/s")
    print(f"  Target:                   {result['target_per_s']:>12,} calc/s")
    print(f"  Status:                   {'MET' if result['meets_target'] else 'NOT MET'}")

    print(f"\n{'=' * 72}")
    print("PRODUCTION SCALING v8 COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
