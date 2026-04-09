"""
Production Scaling v4 — 100k calc/s Benchmark Harness with Timing Instrumentation

Session 204 | PAPER_EXP_BENCH | Daniel Murphy
PURPOSE: Benchmark harness measuring calculations/second for vectorized UQFF
         computation. Profiles 26-layer gravity, F_U_Bi_i assembly, BSFG
         geodesic tracing, and VDS/BSH series. Target: 100k calc/s.

ARCHITECTURE:
  QCalc_Performance.py (vectorization) ──┐
  bsfg_wormhole_geodesic.py (geodesics)  ├──→ this module ──→ benchmark report
  vds_dvp_bsh_lenr_synthesis.py (VDS)  ──┘

BENCHMARKS:
  1. 26-Layer Gravity (Ug1–Ug4 × 26 dimensions)
  2. F_U_Bi_i Assembly (11 force terms from 9-sector E-L)
  3. BSFG Geodesic Step (Christoffel + RK4 integration)
  4. VDS/DVP/BSH Series (Li_26, BSH harmonics, DVP primes)
  5. Full Pipeline (1→4 combined single-system evaluation)
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

KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_UA  = 7.09e-36
RHO_SCM = 7.09e-37
ETA_AETHER = 1e-22
E_REACT = 1e46


# ── §2  BENCHMARK KERNELS ──────────────────────────────────────────────────

def kernel_26layer_gravity(M: float, r: float, mu_s: float,
                           omega_s: float, t: float,
                           ssq: float, rho_scm: float, rho_ua: float,
                           beta_i: float, omega_g: float, d_g: float,
                           ua: float, t_n: float) -> Dict[str, float]:
    """26-layer compressed gravity: g(r,t) = Σ_{i=1}^{26} [Ug1_i + Ug2_i + Ug3_i + Ug4_i].

    Each layer has quantum state factors Q_i, [UA]_i, [SCm]_i.
    """
    Ug1_total = 0.0
    Ug2_total = 0.0
    Ug3_total = 0.0
    Ug4_total = 0.0
    Ubi_total = 0.0
    Um_total  = 0.0

    for layer in range(1, 27):
        # Layer-dependent quantum factors
        Q_i   = 1.0 / (1.0 + layer * 0.01)
        UA_i  = ua * Q_i
        SCm_i = H_SCM * Q_i

        # Layer projection
        theta_i = 2 * PI * layer / 26.0
        proj    = math.cos(theta_i)

        # 4 gravity channels per layer
        Ug1_i = G * mu_s ** 2 / max(r, 1e-10) ** 4 * Q_i * proj
        Ug2_i = G * mu_s * E_REACT / (max(r, 1e-10) ** 3 * c ** 2) * Q_i * proj
        Ug3_i = G * M / max(r, 1e-10) ** 2 * math.sin(omega_s * t) * Q_i * proj
        Ug4_i = G * M / max(r, 1e-10) ** 2 * ssq * (rho_scm / max(rho_ua, 1e-50)) * Q_i * proj

        Ug1_total += Ug1_i
        Ug2_total += Ug2_i
        Ug3_total += Ug3_i
        Ug4_total += Ug4_i

        # Buoyancy per layer
        Ug_i_sum = Ug1_i + Ug2_i + Ug3_i + Ug4_i
        Ubi_i = -beta_i * Ug_i_sum * omega_g * M / max(d_g, 1e-10) * UA_i * math.cos(PI * t_n)
        Ubi_total += Ubi_i

        # Magnetic torque per layer
        Um_i = G * M / max(r, 1e-10) ** 2 * SCm_i * proj
        Um_total += Um_i

    return {
        "Ug1": Ug1_total, "Ug2": Ug2_total,
        "Ug3": Ug3_total, "Ug4": Ug4_total,
        "Ubi": Ubi_total, "Um": Um_total,
        "Ug_sum": Ug1_total + Ug2_total + Ug3_total + Ug4_total,
    }


def kernel_fubii_assembly(M: float, r: float, mu_s: float,
                           omega_s: float, t: float,
                           ssq: float, beta_i: float,
                           omega_g: float, d_g: float,
                           ua: float, t_n: float, scm: float,
                           eta: float, ts00: float,
                           chi: float, k_lenr: float,
                           omega_lenr: float, omega_act: float,
                           lambda_act: float, sigma_n: float,
                           alpha_s: float) -> float:
    """F_U_Bi_i = Σ Ug + Σ Ubi + Um + Tr(A) + F_external (11 terms)."""
    r_safe = max(r, 1e-10)

    # 4 gravity channels
    Ug1 = G * mu_s ** 2 / r_safe ** 4
    Ug2 = G * mu_s * E_REACT / (r_safe ** 3 * c ** 2)
    Ug3 = G * M / r_safe ** 2 * math.sin(omega_s * t)
    Ug4 = G * M / r_safe ** 2 * ssq * (RHO_SCM / max(RHO_UA, 1e-50))

    # 4 buoyancy
    Ug_sum = Ug1 + Ug2 + Ug3 + Ug4
    cos_tn = math.cos(PI * t_n)
    Ubi = -beta_i * Ug_sum * omega_g * M / max(d_g, 1e-10) * ua * cos_tn

    # Magnetic torque
    Um = G * M / r_safe ** 2 * scm

    # Aether trace
    A_trace = G * M / r_safe ** 2 * eta * ts00 * cos_tn

    # LENR
    F_lenr = k_lenr * omega_lenr ** 2 * chi
    F_act  = lambda_act * math.cos(omega_act * t)
    F_res  = sigma_n * chi

    # SM
    F_quark = alpha_s / r_safe ** 2

    return Ug_sum + Ubi + Um + A_trace + F_lenr + F_act + F_res + F_quark


def kernel_bsfg_geodesic_step(r: float, M: float, dr: float,
                               b0: float, r0: float,
                               a_trace: float) -> Tuple[float, float, float]:
    """Single RK4 geodesic step in BSFG metric.

    ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + r²(dθ² + sin²θ dφ²)
    Φ(r) = -GM/(rc²)(1 + η A_trace)
    b(r) = b₀(r₀/r)(1 - β_i [SSq])
    """
    r_safe = max(r, 1e-10)

    # Redshift function
    Phi = -G * M / (r_safe * c ** 2) * (1 + ETA_AETHER * a_trace)

    # Shape function
    b_r = b0 * (r0 / r_safe) * (1 - BETA_I * SSQ)

    # Metric component g_rr
    grr = 1.0 / max(1.0 - b_r / r_safe, 1e-30)

    # Effective potential derivative (radial geodesic equation)
    dPhi_dr = G * M / (r_safe ** 2 * c ** 2) * (1 + ETA_AETHER * a_trace)
    db_dr   = -b0 * r0 / r_safe ** 2 * (1 - BETA_I * SSQ)

    # RK4 step for dr/dlambda (simplified null geodesic)
    k1 = dr
    r_mid = r_safe + 0.5 * k1 * 1e-3
    Phi_mid = -G * M / (max(r_mid, 1e-10) * c ** 2) * (1 + ETA_AETHER * a_trace)
    k2 = dr + 0.5 * dPhi_dr * 1e-3
    k3 = dr + 0.5 * dPhi_dr * 1e-3
    k4 = dr + dPhi_dr * 1e-3

    r_new  = r_safe + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0 * 1e-3
    dr_new = dr + dPhi_dr * 1e-3

    return r_new, dr_new, Phi


def kernel_vds_bsh(ssq: float, rho_ratio: float, beta_i: float,
                   m_kg: float, r_m: float) -> Dict[str, float]:
    """VDS Li_26 + BSH(26) + DVP identity."""
    # Li_26([SSq])
    li26 = 0.0
    for n in range(27):
        gamma_n = 1.0 / (1.0 + n * 0.1)
        li26 += (ssq ** n / math.factorial(n)) * rho_ratio * gamma_n

    # BSH(26)
    f_Ub = beta_i * G * m_kg / max(r_m, 1e-10) ** 2
    bsh = 0.0
    sat = 1.0 - math.exp(-ssq * m_kg / M_sun) if m_kg >= M_sun else 1.0 - math.exp(-ssq * m_kg * 1e27)
    for j in range(1, 27):
        bsh += (1.0 / j) * f_Ub * sat * math.cos(2 * PI * j / 26.0)

    # DVP 26! mod 113
    dvp = 1
    for i in range(1, 27):
        dvp = (dvp * i) % 113

    return {"li_26": li26, "bsh_26": bsh, "dvp_26mod113": dvp, "f_Ub": f_Ub}


def kernel_full_pipeline(params: Dict[str, float]) -> Dict[str, Any]:
    """Full pipeline: 26-layer gravity + F_U_Bi_i + BSFG step + VDS/BSH."""
    M     = params.get("M", 4.3e6 * M_sun)
    r     = params.get("r", 2.44e20)
    mu_s  = params.get("mu_s", 1e25)
    omega_s = params.get("omega_s", 1e-7)
    t     = params.get("t", 0.0)
    ssq   = params.get("ssq", SSQ)
    beta_i = params.get("beta_i", BETA_I)
    omega_g = params.get("omega_g", 7.27e-12)
    d_g   = params.get("d_g", 2.44e20)
    ua    = params.get("ua", U_UA)
    t_n   = params.get("t_n", 0.5)

    # 1. 26-layer gravity
    grav = kernel_26layer_gravity(M, r, mu_s, omega_s, t, ssq, RHO_SCM, RHO_UA,
                                   beta_i, omega_g, d_g, ua, t_n)

    # 2. F_U_Bi_i
    fubi = kernel_fubii_assembly(M, r, mu_s, omega_s, t, ssq, beta_i,
                                  omega_g, d_g, ua, t_n, H_SCM,
                                  ETA_AETHER, 1.27e3,
                                  1e-15, 1.0, 1.25e12, 300 * 2 * PI,
                                  1e-20, 1e-28, 0.118)

    # 3. BSFG geodesic step
    r_new, dr_new, Phi = kernel_bsfg_geodesic_step(r, M, 1e5, 1e4, 1e4, 4.04)

    # 4. VDS/BSH
    vds = kernel_vds_bsh(ssq, RHO_SCM / max(RHO_UA, 1e-50), beta_i, M, r)

    return {
        "Ug_sum": grav["Ug_sum"],
        "F_U_Bi_i": fubi,
        "Phi_BSFG": Phi,
        "li_26": vds["li_26"],
        "bsh_26": vds["bsh_26"],
    }


# ── §3  VECTORIZED KERNELS (numpy) ────────────────────────────────────────

def kernel_26layer_gravity_vectorized(
        M: np.ndarray, r: np.ndarray, mu_s: float,
        omega_s: float, t: float, ssq: float) -> np.ndarray:
    """Vectorized 26-layer gravity over N systems. Returns Ug_sum array."""
    r_safe = np.maximum(r, 1e-10)
    Ug_sum = np.zeros_like(M, dtype=np.float64)

    for layer in range(1, 27):
        Q_i = 1.0 / (1.0 + layer * 0.01)
        theta_i = 2 * PI * layer / 26.0
        proj = math.cos(theta_i)

        Ug1 = G * mu_s ** 2 / r_safe ** 4 * Q_i * proj
        Ug2 = G * mu_s * E_REACT / (r_safe ** 3 * c ** 2) * Q_i * proj
        Ug3 = G * M / r_safe ** 2 * math.sin(omega_s * t) * Q_i * proj
        Ug4 = G * M / r_safe ** 2 * ssq * (RHO_SCM / max(RHO_UA, 1e-50)) * Q_i * proj

        Ug_sum += Ug1 + Ug2 + Ug3 + Ug4

    return Ug_sum


# ── §4  BENCHMARK HARNESS ─────────────────────────────────────────────────

@dataclass
class BenchmarkResult:
    """Result of a single benchmark run."""
    name: str
    n_iterations: int
    total_seconds: float
    calc_per_second: float
    avg_us: float           # microseconds per call
    min_us: float
    max_us: float
    std_us: float
    target_met: bool        # ≥ 100k calc/s ?

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "n_iterations": self.n_iterations,
            "total_seconds": round(self.total_seconds, 4),
            "calc_per_second": round(self.calc_per_second, 1),
            "avg_us": round(self.avg_us, 3),
            "min_us": round(self.min_us, 3),
            "max_us": round(self.max_us, 3),
            "std_us": round(self.std_us, 3),
            "target_met": self.target_met,
        }


class ProductionBenchmark:
    """Production benchmark harness for UQFF computation kernels.

    Measures wall-clock throughput for each kernel, reports calc/s,
    and checks against 100k calc/s target.
    """

    TARGET_CALC_PER_SEC = 100_000

    def __init__(self, n_warmup: int = 100, n_iterations: int = 10_000):
        self.n_warmup = n_warmup
        self.n_iterations = n_iterations
        self.results: List[BenchmarkResult] = []

    def _bench_kernel(self, name: str, kernel_fn: Callable, n: Optional[int] = None) -> BenchmarkResult:
        """Benchmark a single kernel function."""
        N = n or self.n_iterations

        # Warmup
        for _ in range(self.n_warmup):
            kernel_fn()

        # Timed runs — collect per-iteration times
        times_us: List[float] = []
        t0 = time.perf_counter()
        for _ in range(N):
            t_start = time.perf_counter_ns()
            kernel_fn()
            t_end = time.perf_counter_ns()
            times_us.append((t_end - t_start) / 1000.0)
        total = time.perf_counter() - t0

        calc_per_sec = N / max(total, 1e-15)
        avg_us = statistics.mean(times_us)
        min_us = min(times_us)
        max_us = max(times_us)
        std_us = statistics.stdev(times_us) if len(times_us) > 1 else 0.0

        result = BenchmarkResult(
            name=name,
            n_iterations=N,
            total_seconds=total,
            calc_per_second=calc_per_sec,
            avg_us=avg_us,
            min_us=min_us,
            max_us=max_us,
            std_us=std_us,
            target_met=calc_per_sec >= self.TARGET_CALC_PER_SEC,
        )
        self.results.append(result)
        return result

    def run_all(self) -> Dict[str, Any]:
        """Run all benchmark kernels and produce report."""
        self.results.clear()
        t0_all = time.perf_counter()

        # Default parameters for Sgr A*
        M = 4.3e6 * M_sun
        r = 2.44e20
        mu_s = 1e25
        omega_s = 1e-7
        t = 0.0
        beta_i = BETA_I
        omega_g = 7.27e-12
        d_g = 2.44e20
        ua = U_UA
        t_n = 0.5

        # ── Benchmark 1: 26-Layer Gravity ──
        self._bench_kernel(
            "26-Layer Gravity",
            lambda: kernel_26layer_gravity(
                M, r, mu_s, omega_s, t, SSQ, RHO_SCM, RHO_UA,
                beta_i, omega_g, d_g, ua, t_n),
        )

        # ── Benchmark 2: F_U_Bi_i Assembly ──
        self._bench_kernel(
            "F_U_Bi_i Assembly",
            lambda: kernel_fubii_assembly(
                M, r, mu_s, omega_s, t, SSQ, beta_i,
                omega_g, d_g, ua, t_n, H_SCM,
                ETA_AETHER, 1.27e3,
                1e-15, 1.0, 1.25e12, 300 * 2 * PI,
                1e-20, 1e-28, 0.118),
        )

        # ── Benchmark 3: BSFG Geodesic Step ──
        self._bench_kernel(
            "BSFG Geodesic Step",
            lambda: kernel_bsfg_geodesic_step(r, M, 1e5, 1e4, 1e4, 4.04),
        )

        # ── Benchmark 4: VDS/DVP/BSH Series ──
        self._bench_kernel(
            "VDS/DVP/BSH Series",
            lambda: kernel_vds_bsh(SSQ, RHO_SCM / max(RHO_UA, 1e-50), beta_i, M, r),
        )

        # ── Benchmark 5: Full Pipeline ──
        params = {"M": M, "r": r, "mu_s": mu_s, "omega_s": omega_s,
                  "t": t, "ssq": SSQ, "beta_i": beta_i,
                  "omega_g": omega_g, "d_g": d_g, "ua": ua, "t_n": t_n}
        self._bench_kernel(
            "Full Pipeline",
            lambda: kernel_full_pipeline(params),
        )

        # ── Benchmark 6: Vectorized 26-Layer (numpy) ──
        if HAS_NUMPY:
            N_batch = 1000
            M_arr = np.full(N_batch, M)
            r_arr = np.full(N_batch, r)
            self._bench_kernel(
                f"Vectorized 26-Layer ({N_batch} systems)",
                lambda: kernel_26layer_gravity_vectorized(
                    M_arr, r_arr, mu_s, omega_s, t, SSQ),
                n=1000,  # fewer iterations since each call does 1000 systems
            )

        total_elapsed = time.perf_counter() - t0_all

        # Aggregate
        all_met = all(r.target_met for r in self.results)
        fastest = max(self.results, key=lambda r: r.calc_per_second)
        slowest = min(self.results, key=lambda r: r.calc_per_second)

        report = {
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

        return report

    def export_report(self, report: Dict[str, Any],
                      filepath: str = "production_benchmark_v4_report.json") -> str:
        """Export benchmark report to JSON."""
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(report, f, indent=2)
        return filepath


# ── §5  MAIN ────────────────────────────────────────────────────────────────

def main():
    """Run production benchmark suite."""
    print("=" * 72)
    print("Production Scaling v4 — UQFF Benchmark Harness (Target: 100k calc/s)")
    print("=" * 72)
    print(f"  Python:  {sys.version.split()[0]}")
    print(f"  NumPy:   {'yes' if HAS_NUMPY else 'no'}")
    print(f"  Target:  {ProductionBenchmark.TARGET_CALC_PER_SEC:,} calc/s")

    bench = ProductionBenchmark(n_warmup=100, n_iterations=10_000)

    print(f"\nRunning {len(['26L', 'FUBi', 'BSFG', 'VDS', 'Full']) + (1 if HAS_NUMPY else 0)} benchmarks "
          f"({bench.n_iterations:,} iterations each, {bench.n_warmup} warmup)...\n")

    report = bench.run_all()

    # Print results table
    print(f"{'Kernel':<35s} {'calc/s':>12s} {'avg μs':>10s} {'min μs':>10s} {'target':>8s}")
    print("-" * 77)
    for r in report["benchmarks"]:
        status = "  PASS" if r["target_met"] else " *FAIL"
        print(f"{r['name']:<35s} {r['calc_per_second']:>12,.1f} {r['avg_us']:>10.3f} "
              f"{r['min_us']:>10.3f} {status:>8s}")

    print("-" * 77)
    print(f"{'All targets met:':<35s} {'YES' if report['all_targets_met'] else 'NO':>12s}")
    print(f"{'Fastest:':<35s} {report['fastest_kernel']:>35s}  ({report['fastest_calc_per_sec']:,.1f} calc/s)")
    print(f"{'Slowest:':<35s} {report['slowest_kernel']:>35s}  ({report['slowest_calc_per_sec']:,.1f} calc/s)")
    print(f"{'Total elapsed:':<35s} {report['total_elapsed_seconds']:>12.3f} s")

    # Export
    out_path = bench.export_report(report)
    print(f"\n[OK] Report exported: {out_path}")

    print(f"\n{'=' * 72}")
    print(f"BENCHMARK COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
