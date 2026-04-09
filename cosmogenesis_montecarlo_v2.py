"""
Cosmogenesis Monte Carlo v2 — VDS/DVP/BSH Sampling Constrained by GW Data

Session 204 | PAPER_EXP_MC | Daniel Murphy
PURPOSE: Monte Carlo integration of Vacuum Density Series (VDS), Dipole Vortex
         Primes (DVP), and Buoyancy Shell Harmonics (BSH) parameters constrained
         by gravitational-wave observational data from PAPER_001–009.

ARCHITECTURE:
  PAPER_001–009 (GW damping data)  ──┐
  vds_dvp_bsh_lenr_synthesis.py ─────┼──→ this module ──→ posterior distributions
  CondensedPhysics2.py MC wrapper ───┘                    on ρ_SCm/ρ_UA ratio

GW CONSTRAINTS (PAPER_001–009):
  D_TRZ      = 0.900   (10% reduction, topological resonance zone)
  D_String   = 0.370   (63% reduction, dominant string channel)
  F_combined = 0.333   (66.7% strain reduction, product of all D channels)
  GW170817:  h_UQFF/h_GR = 0.333,  h_GR = 5.4176e-22
  GW150914:  h_UQFF/h_GR = 0.333,  h_GR = 1.2499e-21
  GW190425:  D_total = 0.530  (mass-gap BNS, weaker damping)
"""

import math
import json
import os
import sys
import time
import hashlib
from typing import Dict, List, Tuple, Optional, Any
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

# UQFF calibrated
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_UA  = 7.09e-36
RHO_SCM = 7.09e-37
ETA_AETHER = 1e-22

# VDS/DVP/BSH (from vds_dvp_bsh_lenr_synthesis.py)
RHO_VAC_SCM = 9.47e-27
RHO_VAC_UA  = 5e-27
BASE_RATIO  = RHO_VAC_SCM / RHO_VAC_UA  # ≈ 1.894

# GW damping (PAPER_001–009) — OBSERVATIONAL CONSTRAINTS
D_AETHER    = 1.000
D_SCM_GW    = 1.000
D_TRZ       = 0.900
D_STRING    = 0.370
F_COMBINED  = D_AETHER * D_SCM_GW * D_TRZ * D_STRING  # 0.333

# DVP primes
PRIMES_30 = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
             53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]
MAGIC_NUMBERS = {2, 8, 20, 28, 50, 82, 126}


# ── §2  GW EVENT DATA ──────────────────────────────────────────────────────

@dataclass
class GWEvent:
    """Gravitational-wave event observational constraint."""
    name: str
    m1_Msun: float
    m2_Msun: float
    DL_Mpc: float
    D_total: float
    h_GR: float
    h_UQFF: float
    event_type: str  # "BNS" or "BBH"

GW_EVENTS = [
    GWEvent("GW170817", 1.46, 1.27, 40.0, 0.333, 5.4176e-22, 1.804e-22, "BNS"),
    GWEvent("GW150914", 36.0, 29.0, 410.0, 0.810, 1.2499e-21, 4.1622e-22, "BBH"),
    GWEvent("GW190425", 1.7, 1.5, 159.0, 0.530, 3.0e-22, 1.59e-22, "BNS"),
]


# ── §3  VDS / DVP / BSH PHYSICS FUNCTIONS ──────────────────────────────────

def li_26_ssq(ssq: float, base_ratio: float) -> float:
    """Li_26([SSq]) — 26-term vacuum density polylogarithm series."""
    total = 0.0
    for n in range(27):
        gamma_n = 1.0 / (1.0 + n * 0.1)
        total += (ssq ** n / math.factorial(n)) * base_ratio * gamma_n
    return total


def dvp_prime(Z: int) -> int:
    """DVP prime for atomic number Z."""
    idx = min(Z - 1, len(PRIMES_30) - 1)
    return PRIMES_30[max(idx, 0)]


def nuclear_identity(mod_base: int = 113) -> int:
    """26! mod 113 = 12  (DVP proto-shell identity)."""
    result = 1
    for i in range(1, 27):
        result = (result * i) % mod_base
    return result


def bsh_harmonic(j: int, m_kg: float, f_Ub: float, ssq: float = SSQ) -> float:
    """Single BSH term at harmonic j."""
    if m_kg >= M_sun:
        saturation = 1.0 - math.exp(-ssq * m_kg / M_sun)
    else:
        saturation = 1.0 - math.exp(-ssq * m_kg * 1e27)
    layer_proj = math.cos(2 * PI * j / 26.0)
    return (1.0 / j) * f_Ub * saturation * layer_proj


def bsh_sum(k: int, m_kg: float, f_Ub: float, ssq: float = SSQ) -> float:
    """Cumulative BSH up to harmonic k."""
    return sum(bsh_harmonic(j, m_kg, f_Ub, ssq) for j in range(1, k + 1))


def compute_d_string(ssq: float, rho_ratio: float) -> float:
    """D_String damping factor from VDS ratio and [SSq].

    Model: D_String = SSq × (1 - rho_SCm/rho_total) maps VDS vacuum
    structure to GW string-channel damping.
    """
    rho_total = 1.0 + rho_ratio  # normalized: rho_UA=1
    return ssq * (1.0 - rho_ratio / rho_total)


def compute_d_trz(kappa: float, t_observation: float) -> float:
    """D_TRZ damping from time-reversal zone decay.

    Model: D_TRZ = exp(-κ × t_eff)  where t_eff captures
    topological resonance zone lifetime.
    """
    return math.exp(-kappa * t_observation)


def compute_d_total(d_aether: float, d_scm: float, d_trz: float, d_string: float) -> float:
    """Product damping: D_total = D_Aether × D_SCm × D_TRZ × D_String."""
    return d_aether * d_scm * d_trz * d_string


def h_uqff_from_gr(h_gr: float, d_total: float) -> float:
    """h_UQFF = h_GR × D_total."""
    return h_gr * d_total


def chirp_mass(m1_kg: float, m2_kg: float) -> float:
    """Chirp mass M_c = (m1 m2)^(3/5) / (m1 + m2)^(1/5)."""
    return (m1_kg * m2_kg) ** 0.6 / (m1_kg + m2_kg) ** 0.2


# ── §4  MONTE CARLO ENGINE ────────────────────────────────────────────────

@dataclass
class MCConfig:
    """Monte Carlo configuration."""
    n_samples: int = 100_000
    seed: int = 42
    ssq_range: Tuple[float, float] = (0.40, 0.75)
    rho_ratio_range: Tuple[float, float] = (1.0, 3.0)
    kappa_range: Tuple[float, float] = (1e-10, 1e-7)
    beta_i_range: Tuple[float, float] = (0.3, 0.9)
    # GW constraint tolerances
    f_combined_target: float = 0.333
    f_combined_tol: float = 0.05      # ±5% acceptance
    d_trz_target: float = 0.900
    d_trz_tol: float = 0.05
    d_string_target: float = 0.370
    d_string_tol: float = 0.05


class CosmogenesisMonteCarloV2:
    """GW-constrained Monte Carlo over VDS/DVP/BSH parameter space.

    Samples [SSq], ρ_SCm/ρ_UA ratio, κ, β_i uniformly within prior ranges,
    then rejects samples that produce D_total outside the observational
    acceptance window around F_combined = 0.333 (PAPER_001 GW170817).
    """

    def __init__(self, config: Optional[MCConfig] = None):
        self.config = config or MCConfig()
        if HAS_NUMPY:
            self.rng = np.random.default_rng(self.config.seed)
        else:
            import random
            self._pyrng = random.Random(self.config.seed)

    # ── sampling ──

    def _uniform(self, low: float, high: float) -> float:
        if HAS_NUMPY:
            return float(self.rng.uniform(low, high))
        return self._pyrng.uniform(low, high)

    def _uniform_batch(self, low: float, high: float, n: int):
        if HAS_NUMPY:
            return self.rng.uniform(low, high, size=n)
        return [self._pyrng.uniform(low, high) for _ in range(n)]

    # ── single sample ──

    def sample_and_evaluate(self) -> Optional[Dict[str, float]]:
        """Draw one parameter sample, compute damping, apply GW constraint.

        Returns dict if accepted, None if rejected.
        """
        cfg = self.config

        ssq       = self._uniform(*cfg.ssq_range)
        rho_ratio = self._uniform(*cfg.rho_ratio_range)
        kappa     = self._uniform(*cfg.kappa_range)
        beta_i    = self._uniform(*cfg.beta_i_range)

        # t_observation calibrated so D_TRZ(kappa_nominal) ≈ 0.900
        # ln(0.900) / (-kappa) ≈ 1.82e7 for kappa=5.787e-9
        t_obs = -math.log(cfg.d_trz_target) / KAPPA  # reference timescale

        d_trz    = compute_d_trz(kappa, t_obs)
        d_string = compute_d_string(ssq, rho_ratio)
        d_total  = compute_d_total(D_AETHER, D_SCM_GW, d_trz, d_string)

        # Acceptance: D_total within tolerance of F_combined target
        if abs(d_total - cfg.f_combined_target) > cfg.f_combined_tol:
            return None

        # VDS computation
        li26 = li_26_ssq(ssq, rho_ratio)

        # BSH for representative stellar mass (1 M_sun, 1 AU)
        f_Ub = beta_i * G * M_sun / (1.496e11) ** 2
        bsh_26_val = bsh_sum(26, M_sun, f_Ub, ssq)

        # DVP identity check
        dvp_12 = nuclear_identity()

        return {
            "ssq": ssq,
            "rho_ratio": rho_ratio,
            "kappa": kappa,
            "beta_i": beta_i,
            "d_trz": d_trz,
            "d_string": d_string,
            "d_total": d_total,
            "li_26": li26,
            "bsh_26": bsh_26_val,
            "f_Ub": f_Ub,
            "dvp_26mod113": dvp_12,
        }

    # ── full ensemble ──

    def run(self, n_samples: Optional[int] = None) -> Dict[str, Any]:
        """Run full MC ensemble with GW-constrained rejection sampling.

        Returns statistics, accepted samples, and convergence diagnostics.
        """
        N = n_samples or self.config.n_samples
        cfg = self.config

        t0 = time.perf_counter()

        accepted: List[Dict[str, float]] = []
        rejected = 0

        if HAS_NUMPY:
            # Vectorized batch sampling + serial constraint check
            ssq_all       = self._uniform_batch(*cfg.ssq_range, N)
            rho_ratio_all = self._uniform_batch(*cfg.rho_ratio_range, N)
            kappa_all     = self._uniform_batch(*cfg.kappa_range, N)
            beta_i_all    = self._uniform_batch(*cfg.beta_i_range, N)

            t_obs = -math.log(cfg.d_trz_target) / KAPPA

            # Vectorized damping computation
            d_trz_all    = np.exp(-kappa_all * t_obs)
            rho_total    = 1.0 + rho_ratio_all
            d_string_all = ssq_all * (1.0 - rho_ratio_all / rho_total)
            d_total_all  = D_AETHER * D_SCM_GW * d_trz_all * d_string_all

            # Acceptance mask
            mask = np.abs(d_total_all - cfg.f_combined_target) <= cfg.f_combined_tol
            accepted_idx = np.where(mask)[0]
            rejected = int(N - len(accepted_idx))

            for i in accepted_idx:
                ssq_i = float(ssq_all[i])
                rr_i  = float(rho_ratio_all[i])
                ki    = float(kappa_all[i])
                bi    = float(beta_i_all[i])

                li26 = li_26_ssq(ssq_i, rr_i)
                f_Ub = bi * G * M_sun / (1.496e11) ** 2
                bsh_val = bsh_sum(26, M_sun, f_Ub, ssq_i)

                accepted.append({
                    "ssq": ssq_i,
                    "rho_ratio": rr_i,
                    "kappa": ki,
                    "beta_i": bi,
                    "d_trz": float(d_trz_all[i]),
                    "d_string": float(d_string_all[i]),
                    "d_total": float(d_total_all[i]),
                    "li_26": li26,
                    "bsh_26": bsh_val,
                    "f_Ub": f_Ub,
                    "dvp_26mod113": 12,
                })
        else:
            for _ in range(N):
                result = self.sample_and_evaluate()
                if result is not None:
                    accepted.append(result)
                else:
                    rejected += 1

        elapsed = time.perf_counter() - t0

        # Compute posterior statistics
        stats = self._compute_statistics(accepted)

        # Per-event validation
        event_validation = self._validate_gw_events(stats)

        return {
            "config": {
                "n_samples": N,
                "seed": self.config.seed,
                "ssq_range": list(self.config.ssq_range),
                "rho_ratio_range": list(self.config.rho_ratio_range),
                "kappa_range": list(self.config.kappa_range),
                "beta_i_range": list(self.config.beta_i_range),
                "f_combined_target": cfg.f_combined_target,
                "f_combined_tol": cfg.f_combined_tol,
            },
            "n_accepted": len(accepted),
            "n_rejected": rejected,
            "acceptance_rate": len(accepted) / max(N, 1),
            "elapsed_seconds": elapsed,
            "samples_per_second": N / max(elapsed, 1e-9),
            "statistics": stats,
            "gw_event_validation": event_validation,
            "convergence": self._convergence_diagnostic(accepted),
            "dvp_identity": {"26! mod 113": 12, "verified": nuclear_identity() == 12},
        }

    # ── statistics ──

    def _compute_statistics(self, samples: List[Dict[str, float]]) -> Dict[str, Any]:
        """Posterior statistics for all sampled parameters."""
        if not samples:
            return {"error": "no accepted samples"}

        keys = ["ssq", "rho_ratio", "kappa", "beta_i",
                "d_trz", "d_string", "d_total", "li_26", "bsh_26"]
        stats = {}

        if HAS_NUMPY:
            for key in keys:
                values = np.array([s[key] for s in samples])
                stats[key] = {
                    "mean": float(np.mean(values)),
                    "std": float(np.std(values)),
                    "median": float(np.median(values)),
                    "ci_lower_95": float(np.percentile(values, 2.5)),
                    "ci_upper_95": float(np.percentile(values, 97.5)),
                    "p05": float(np.percentile(values, 5)),
                    "p95": float(np.percentile(values, 95)),
                    "min": float(np.min(values)),
                    "max": float(np.max(values)),
                }
        else:
            for key in keys:
                values = sorted([s[key] for s in samples])
                n = len(values)
                mean = sum(values) / n
                var = sum((v - mean) ** 2 for v in values) / max(n - 1, 1)
                std = var ** 0.5
                stats[key] = {
                    "mean": mean,
                    "std": std,
                    "median": values[n // 2],
                    "ci_lower_95": values[max(0, int(0.025 * n))],
                    "ci_upper_95": values[min(n - 1, int(0.975 * n))],
                    "min": values[0],
                    "max": values[-1],
                }

        return stats

    def _validate_gw_events(self, stats: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Validate posterior against each GW event."""
        if "error" in stats:
            return [{"error": "no statistics available"}]

        results = []
        d_mean = stats.get("d_total", {}).get("mean", 0.333)

        for ev in GW_EVENTS:
            h_predicted = ev.h_GR * d_mean
            residual = abs(h_predicted - ev.h_UQFF) / max(ev.h_UQFF, 1e-30)
            results.append({
                "event": ev.name,
                "type": ev.event_type,
                "D_total_observed": ev.D_total,
                "D_total_posterior_mean": d_mean,
                "h_GR": ev.h_GR,
                "h_UQFF_observed": ev.h_UQFF,
                "h_UQFF_predicted": h_predicted,
                "residual_fraction": residual,
                "consistent": residual < 0.15,  # 15% tolerance
            })

        return results

    def _convergence_diagnostic(self, samples: List[Dict[str, float]]) -> Dict[str, Any]:
        """Assess MC convergence by comparing first/second half statistics."""
        n = len(samples)
        if n < 20:
            return {"converged": False, "reason": f"too few samples ({n})"}

        half = n // 2
        first_half = samples[:half]
        second_half = samples[half:]

        def _mean(lst, key):
            return sum(s[key] for s in lst) / len(lst)

        key = "d_total"
        mean_1 = _mean(first_half, key)
        mean_2 = _mean(second_half, key)
        drift = abs(mean_1 - mean_2) / max(abs(mean_1), 1e-30)

        return {
            "converged": drift < 0.02,
            "d_total_mean_first_half": mean_1,
            "d_total_mean_second_half": mean_2,
            "drift_fraction": drift,
            "n_accepted": n,
        }

    # ── export ──

    def export_results(self, results: Dict[str, Any],
                       filepath: str = "cosmogenesis_mc_v2_results.json") -> str:
        """Export MC results to JSON."""
        # Remove non-serializable values
        clean = json.loads(json.dumps(results, default=str))
        with open(filepath, "w", encoding="utf-8") as f:
            json.dump(clean, f, indent=2)
        return filepath


# ── §5  MAIN ────────────────────────────────────────────────────────────────

def main():
    """Run GW-constrained Monte Carlo on VDS/DVP/BSH parameters."""
    print("=" * 72)
    print("Cosmogenesis Monte Carlo v2 — VDS × GW Constrained Sampling")
    print("=" * 72)

    # Configure
    config = MCConfig(
        n_samples=100_000,
        seed=42,
        ssq_range=(0.40, 0.75),
        rho_ratio_range=(1.0, 3.0),
        kappa_range=(1e-10, 1e-7),
        beta_i_range=(0.3, 0.9),
        f_combined_target=0.333,
        f_combined_tol=0.05,
    )

    engine = CosmogenesisMonteCarloV2(config)

    print(f"\nSampling {config.n_samples:,} parameter sets...")
    print(f"  [SSq]     ∈ [{config.ssq_range[0]}, {config.ssq_range[1]}]")
    print(f"  ρ_ratio   ∈ [{config.rho_ratio_range[0]}, {config.rho_ratio_range[1]}]")
    print(f"  κ         ∈ [{config.kappa_range[0]:.1e}, {config.kappa_range[1]:.1e}]")
    print(f"  β_i       ∈ [{config.beta_i_range[0]}, {config.beta_i_range[1]}]")
    print(f"  F_target  = {config.f_combined_target} ± {config.f_combined_tol}")

    results = engine.run()

    # Report
    print(f"\n── Results ──")
    print(f"  Accepted:        {results['n_accepted']:,} / {config.n_samples:,}")
    print(f"  Acceptance rate: {results['acceptance_rate']:.4f}")
    print(f"  Elapsed:         {results['elapsed_seconds']:.2f} s")
    print(f"  Throughput:      {results['samples_per_second']:,.0f} samples/s")

    stats = results["statistics"]
    if "error" not in stats:
        print(f"\n── Posterior Distributions ──")
        for key in ["ssq", "rho_ratio", "d_total", "d_trz", "d_string", "li_26"]:
            s = stats[key]
            print(f"  {key:12s}:  μ={s['mean']:.6f}  σ={s['std']:.6f}  "
                  f"95%CI=[{s['ci_lower_95']:.6f}, {s['ci_upper_95']:.6f}]")

        print(f"\n── GW Event Validation ──")
        for ev in results["gw_event_validation"]:
            status = "PASS" if ev.get("consistent", False) else "FAIL"
            print(f"  {ev['event']:10s} ({ev['type']}):  "
                  f"h_pred={ev.get('h_UQFF_predicted', 0):.4e}  "
                  f"h_obs={ev.get('h_UQFF_observed', 0):.4e}  "
                  f"resid={ev.get('residual_fraction', 0):.4f}  [{status}]")

        conv = results["convergence"]
        print(f"\n── Convergence ──")
        print(f"  Converged: {conv['converged']}")
        print(f"  Drift:     {conv.get('drift_fraction', 0):.6f}")

    # DVP identity
    dvp = results["dvp_identity"]
    print(f"\n── DVP Identity ──")
    print(f"  26! mod 113 = {dvp['26! mod 113']}  (verified: {dvp['verified']})")

    # Export
    out_path = engine.export_results(results)
    print(f"\n[OK] Results exported: {out_path}")

    print(f"\n{'=' * 72}")
    print(f"MONTE CARLO COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
