#!/usr/bin/env python3
"""
Production Scaling v7 — 300k calc/s Benchmark Harness with Phonon + Jet + NS Batch Kernels

Session 211 | PAPER_EXP_BENCH_V7 | Daniel Murphy
PURPOSE: Benchmark harness measuring calculations/second for the extended UQFF
         computation pipeline including phonon resonance, jet modulation, and
         NS spindown kernels. Target: 300k calc/s (3× v4 target).

ARCHITECTURE:
  production_scaling_v4.py (100k baseline) ──┐
  scm_phonon_resonance.py (phonon)           ├──→ this module ──→ benchmark report
  quasar_jet_phonon.py (jets)                │
  ns_phonon_gw190425_wstp.py (NS/GW)  ──────┘

BENCHMARKS:
  1. 26-Layer Gravity (baseline, from v4)
  2. F_U_Bi_i Assembly (baseline, from v4)
  3. Phonon Resonance a_res (new)
  4. Jet Modulation M_jet(Γ) sweep (new)
  5. NS Spindown Ω̇ phonon correction (new)
  6. GW190425 h_UQFF(t) strain (new)
  7. Full Pipeline v7 (all kernels combined)
  8. Vectorized Phonon Batch (numpy, new)

TARGET: 300,000 calc/s (3× v4's 100k)
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

OMEGA_SCM = 2 * PI * 1.25e12  # 1.25 THz
GAMMA_SCM = 2 * PI * 0.1e12   # 0.1 THz linewidth
PHI_0     = 1e20               # phonon fluence

A_JET_DEFAULT = 1.5
SIGMA_G_RAD = 0.08 * 2 * PI * 1e12


# ── §2  BENCHMARK KERNELS (v4 baseline) ────────────────────────────────────

def kernel_26layer_gravity(M: float, r: float, mu_s: float,
                           omega_s: float, t: float) -> Dict[str, float]:
    """26-layer compressed gravity (from v4)."""
    Ug_total = 0.0
    Ubi_total = 0.0
    for layer in range(1, 27):
        Q_i = 1.0 / (1.0 + layer * 0.01)
        theta_i = 2 * PI * layer / 26.0
        proj = math.cos(theta_i)
        r_safe = max(r, 1e-10)

        Ug1_i = G * mu_s**2 / r_safe**4 * Q_i * proj
        Ug2_i = G * mu_s * E_REACT / (r_safe**3 * c**2) * Q_i * proj
        Ug3_i = G * M / r_safe**2 * math.sin(omega_s * t) * Q_i * proj
        Ug4_i = G * M / r_safe**2 * SSQ * (RHO_SCM / max(RHO_UA, 1e-50)) * Q_i * proj

        Ug_i = Ug1_i + Ug2_i + Ug3_i + Ug4_i
        Ug_total += Ug_i

        Ubi_i = -BETA_I * Ug_i * 7.27e-12 * M / max(r_safe, 1e-10) * U_UA * Q_i * math.cos(PI * 0.5)
        Ubi_total += Ubi_i

    return {"Ug_sum": Ug_total, "Ubi": Ubi_total}


def kernel_fubii_assembly(M: float, r: float, mu_s: float,
                          omega_s: float, t: float) -> Dict[str, float]:
    """F_U_Bi_i force assembly (simplified from v4)."""
    grav = kernel_26layer_gravity(M, r, mu_s, omega_s, t)
    Ug_sum = grav["Ug_sum"]
    Ubi = grav["Ubi"]

    # Additional force terms
    F_LENR = 1e-15 * math.sin(1.25e12 * t + KAPPA * r)
    F_mag = mu_s * 1e-7 / max(r, 1e-10)**3
    F_aether = ETA_AETHER * c * max(r, 1e-10)

    F_total = abs(Ug_sum) + abs(Ubi) + abs(F_LENR) + abs(F_mag) + abs(F_aether)
    return {"F_total": F_total, "Ug_sum": Ug_sum, "Ubi": Ubi}


# ── §3  NEW PHONON + JET + NS KERNELS ──────────────────────────────────────

def kernel_phonon_a_res(F_UBi: float, F_U: float, omega: float,
                        Gamma: float) -> Dict[str, float]:
    """Phonon resonance acceleration a_res (from scm_phonon_resonance.py)."""
    ssq = SSQ
    # S₂₆ approximation (Li₂₆ evaluation)
    S26 = 0.57  # PolyLog[26, 0.57] ≈ 0.57 for fast eval

    delta_omega = omega - OMEGA_SCM
    gaussian = math.exp(-delta_omega**2 / (2.0 * Gamma**2))
    Phi_norm = gaussian * S26

    ratio = F_UBi / max(F_U, 1e-50)
    a_res = ratio * Phi_norm * S26

    Q = omega / (2.0 * Gamma) if Gamma > 0 else 0.0
    return {"a_res": a_res, "Phi": Phi_norm, "Q": Q}


def kernel_jet_modulation(Gamma: float, A_jet: float = A_JET_DEFAULT) -> Dict[str, float]:
    """Jet modulation M_jet(Γ) (from quasar_jet_phonon.py)."""
    delta_G = Gamma - GAMMA_SCM
    gauss = math.exp(-delta_G**2 / (2.0 * SIGMA_G_RAD**2))
    M_jet = 1.0 + A_jet * gauss

    # BZ power (generic 10⁹ M☉ BH)
    M_bh = 1e9 * M_sun
    r_S = 2.0 * G * M_bh / c**2
    r_H = r_S / 2.0 * (1.0 + math.sqrt(max(1.0 - 0.9**2, 0.0)))
    B = 1e4
    P_BZ = (B**2 / (8.0 * PI)) * (r_H / c)**2 * 0.9**2 * c
    P_jet = P_BZ * (1.0 + M_jet)

    return {"M_jet": M_jet, "P_BZ": P_BZ, "P_jet": P_jet}


def kernel_ns_spindown(M_ns: float, R_ns: float, B_surf: float,
                       P_spin: float) -> Dict[str, float]:
    """NS spindown with phonon correction (from ns_phonon_gw190425_wstp.py)."""
    S26 = SSQ  # fast approx
    Phi_norm = S26  # at resonance

    correction = Phi_norm * S26 * SSQ / N_LEVELS
    I_ns = 0.4 * M_ns * R_ns**2
    Omega = 2.0 * PI / P_spin

    Omega_dot = -B_surf**2 * R_ns**6 * Omega**3 / (6.0 * I_ns * c**3)
    Omega_dot_phonon = Omega_dot * (1.0 + correction)

    return {"Omega_dot": Omega_dot, "Omega_dot_phonon": Omega_dot_phonon}


def kernel_gw190425_strain(t: float) -> Dict[str, float]:
    """GW190425 strain h_UQFF(t) (from ns_phonon_gw190425_wstp.py)."""
    h_GR = 3.0e-22
    suppression = 0.5297
    exp_factor = math.exp(SSQ * t / 26.0) if abs(t) < 1e10 else 1.0
    h_UQFF = h_GR * suppression * exp_factor

    return {"h_GR": h_GR, "h_UQFF": h_UQFF, "exp_factor": exp_factor}


def kernel_full_pipeline_v7(params: dict) -> Dict[str, float]:
    """Full v7 pipeline: gravity + F_UBi + phonon + jet + NS + GW strain."""
    M = params.get("M", 4.3e6 * M_sun)
    r = params.get("r", 2.44e20)
    mu_s = params.get("mu_s", 1e25)
    omega_s = params.get("omega_s", 1e-7)
    t = params.get("t", 0.0)

    # Gravity
    grav = kernel_26layer_gravity(M, r, mu_s, omega_s, t)

    # Force assembly
    forces = kernel_fubii_assembly(M, r, mu_s, omega_s, t)

    # Phonon
    phonon = kernel_phonon_a_res(0.6, 1.0, OMEGA_SCM, GAMMA_SCM)

    # Jet
    jet = kernel_jet_modulation(GAMMA_SCM)

    # NS spindown
    ns = kernel_ns_spindown(1.6 * M_sun, 12e3, 1e8, 0.01)

    # GW strain
    gw = kernel_gw190425_strain(t)

    return {
        "Ug_sum": grav["Ug_sum"],
        "F_total": forces["F_total"],
        "a_res": phonon["a_res"],
        "M_jet": jet["M_jet"],
        "Omega_dot_phonon": ns["Omega_dot_phonon"],
        "h_UQFF": gw["h_UQFF"],
    }


def kernel_phonon_batch_vectorized(N: int) -> Dict[str, Any]:
    """Vectorized phonon a_res evaluation (numpy)."""
    if not HAS_NUMPY:
        return {"error": "numpy not available"}

    omegas = np.full(N, OMEGA_SCM)
    Gammas = np.random.uniform(GAMMA_SCM * 0.5, GAMMA_SCM * 2.0, N)
    F_UBis = np.random.uniform(0.3, 0.8, N)

    delta = omegas - OMEGA_SCM
    gaussian = np.exp(-delta**2 / (2.0 * Gammas**2))
    Phi = gaussian * SSQ
    a_res = (F_UBis / 1.0) * Phi * SSQ

    return {"mean_a_res": float(np.mean(a_res)), "N": N}


# ── §4  BENCHMARK FRAMEWORK ────────────────────────────────────────────────

@dataclass
class BenchmarkResult:
    name: str
    n_iterations: int
    total_seconds: float
    calc_per_second: float
    avg_us: float
    min_us: float
    max_us: float
    std_us: float
    target_met: bool

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "n_iterations": self.n_iterations,
            "total_seconds": round(self.total_seconds, 6),
            "calc_per_second": round(self.calc_per_second, 1),
            "avg_us": round(self.avg_us, 3),
            "min_us": round(self.min_us, 3),
            "max_us": round(self.max_us, 3),
            "std_us": round(self.std_us, 3),
            "target_met": self.target_met,
        }


class ProductionBenchmarkV7:
    TARGET_CALC_PER_SEC = 300_000  # 3× v4 target

    def __init__(self, n_warmup: int = 200, n_iterations: int = 50_000):
        self.n_warmup = n_warmup
        self.n_iterations = n_iterations
        self.results: List[BenchmarkResult] = []

    def _bench_kernel(self, name: str, fn: Callable, n: int = None) -> BenchmarkResult:
        n = n or self.n_iterations

        # Warmup
        for _ in range(self.n_warmup):
            fn()

        # Timed iterations
        timings = []
        t0_all = time.perf_counter()
        for _ in range(n):
            t0 = time.perf_counter()
            fn()
            elapsed = time.perf_counter() - t0
            timings.append(elapsed)
        total = time.perf_counter() - t0_all

        timings_us = [t * 1e6 for t in timings]
        calc_per_sec = n / total if total > 0 else 0.0
        avg_us = statistics.mean(timings_us)
        min_us = min(timings_us)
        max_us = max(timings_us)
        std_us = statistics.stdev(timings_us) if len(timings_us) > 1 else 0.0

        result = BenchmarkResult(
            name=name, n_iterations=n,
            total_seconds=total,
            calc_per_second=calc_per_sec,
            avg_us=avg_us, min_us=min_us, max_us=max_us, std_us=std_us,
            target_met=calc_per_sec >= self.TARGET_CALC_PER_SEC,
        )
        self.results.append(result)
        return result

    def run_all(self) -> Dict[str, Any]:
        self.results.clear()
        t0_all = time.perf_counter()

        M = 4.3e6 * M_sun
        r = 2.44e20
        mu_s = 1e25
        omega_s = 1e-7
        t = 0.0

        # 1. 26-Layer Gravity
        self._bench_kernel(
            "26-Layer Gravity",
            lambda: kernel_26layer_gravity(M, r, mu_s, omega_s, t),
        )

        # 2. F_U_Bi_i Assembly
        self._bench_kernel(
            "F_U_Bi_i Assembly",
            lambda: kernel_fubii_assembly(M, r, mu_s, omega_s, t),
        )

        # 3. Phonon a_res
        self._bench_kernel(
            "Phonon a_res",
            lambda: kernel_phonon_a_res(0.6, 1.0, OMEGA_SCM, GAMMA_SCM),
        )

        # 4. Jet M_jet(Γ)
        self._bench_kernel(
            "Jet M_jet(Γ)",
            lambda: kernel_jet_modulation(GAMMA_SCM),
        )

        # 5. NS Spindown Ω̇
        self._bench_kernel(
            "NS Spindown Ω̇_phonon",
            lambda: kernel_ns_spindown(1.6 * M_sun, 12e3, 1e8, 0.01),
        )

        # 6. GW190425 h_UQFF
        self._bench_kernel(
            "GW190425 h_UQFF(t)",
            lambda: kernel_gw190425_strain(0.0),
        )

        # 7. Full Pipeline v7
        params = {"M": M, "r": r, "mu_s": mu_s, "omega_s": omega_s, "t": t}
        self._bench_kernel(
            "Full Pipeline v7",
            lambda: kernel_full_pipeline_v7(params),
        )

        # 8. Vectorized Phonon Batch
        if HAS_NUMPY:
            N_batch = 1000
            self._bench_kernel(
                f"Vectorized Phonon ({N_batch} batch)",
                lambda: kernel_phonon_batch_vectorized(N_batch),
                n=5000,
            )

        total_elapsed = time.perf_counter() - t0_all

        all_met = all(r.target_met for r in self.results)
        fastest = max(self.results, key=lambda r: r.calc_per_second)
        slowest = min(self.results, key=lambda r: r.calc_per_second)

        return {
            "version": "v7",
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "target_calc_per_sec": self.TARGET_CALC_PER_SEC,
            "n_warmup": self.n_warmup,
            "n_iterations": self.n_iterations,
            "total_elapsed_seconds": round(total_elapsed, 3),
            "all_targets_met": all_met,
            "fastest_kernel": fastest.name,
            "fastest_calc_per_sec": round(fastest.calc_per_second, 1),
            "slowest_kernel": slowest.name,
            "slowest_calc_per_sec": round(slowest.calc_per_second, 1),
            "has_numpy": HAS_NUMPY,
            "python_version": sys.version.split()[0],
            "benchmarks": [r.to_dict() for r in self.results],
        }

    def export_report(self, report: Dict[str, Any],
                      filepath: str = "production_benchmark_v7_report.json") -> str:
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        return filepath


# ── §5  MAIN ────────────────────────────────────────────────────────────────

def main():
    print("=" * 72)
    print("Production Scaling v7 — UQFF Benchmark (Target: 300k calc/s)")
    print("=" * 72)
    print(f"  Python:  {sys.version.split()[0]}")
    print(f"  NumPy:   {'yes' if HAS_NUMPY else 'no'}")
    print(f"  Target:  {ProductionBenchmarkV7.TARGET_CALC_PER_SEC:,} calc/s")

    bench = ProductionBenchmarkV7(n_warmup=200, n_iterations=50_000)

    print(f"\nRunning 7+1 benchmarks "
          f"({bench.n_iterations:,} iterations each, {bench.n_warmup} warmup)...\n")

    report = bench.run_all()

    print(f"{'Kernel':<40s} {'calc/s':>12s} {'avg μs':>10s} {'min μs':>10s} {'target':>8s}")
    print("-" * 82)
    for r in report["benchmarks"]:
        status = "  PASS" if r["target_met"] else " *FAIL"
        print(f"{r['name']:<40s} {r['calc_per_second']:>12,.1f} {r['avg_us']:>10.3f} "
              f"{r['min_us']:>10.3f} {status:>8s}")

    print("-" * 82)
    print(f"{'All targets met:':<40s} {'YES' if report['all_targets_met'] else 'NO':>12s}")
    print(f"{'Fastest:':<40s} {report['fastest_kernel']:>30s}  ({report['fastest_calc_per_sec']:,.1f} calc/s)")
    print(f"{'Slowest:':<40s} {report['slowest_kernel']:>30s}  ({report['slowest_calc_per_sec']:,.1f} calc/s)")
    print(f"{'Total elapsed:':<40s} {report['total_elapsed_seconds']:>12.3f} s")

    out_path = bench.export_report(report)
    print(f"\n[OK] Report exported: {out_path}")

    print(f"\n{'=' * 72}")
    print(f"BENCHMARK V7 COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
