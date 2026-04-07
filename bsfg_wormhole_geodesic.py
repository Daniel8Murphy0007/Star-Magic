#!/usr/bin/env python3
"""
bsfg_wormhole_geodesic.py — BSFG Metric Wormhole Geodesic Tracer
═════════════════════════════════════════════════════════════════

PURPOSE: Trace null and timelike geodesics through traversable wormholes using
the updated BSFG (Buoyancy-SCm-Fluid-Geometry) metric, which extends the
Morris-Thorne metric with UQFF vacuum buoyancy and aether tensor terms.

BSFG METRIC (updated):
  ds² = -e^{2Φ(r)} dt² + dr²/(1 - b(r)/r) + (b² + r²)(dθ² + sin²θ dφ²)

  where:
    Φ(r) = -GM/(rc²) × (1 + A_trace)          [redshift function, aether-corrected]
    b(r) = b₀ × (r₀/r) × (1 - β_i SSq)       [shape function, buoyancy-modified]

  Traversability condition (PAPER_877):
    ρ + P = -1.75 × 10⁵ kg/m³  (exotic matter from SCm vacuum)

REFERENCES:
  - CondensedPhysics4.py L1806: MorrisThorneWormholeNullGeodesicsCalculator
  - CondensedPhysics4.py L10769-11223: BSFG 5-calculator cascade (#149-153)
  - source4.cpp L1560: compute_a_wormhole(r, b, f_worm, Evac_neb)
  - PAPER_373: Morris-Thorne null geodesics
  - PAPER_554-558: BSFG complete system
  - PAPER_877: Production wormhole metrics (ρ+P = -1.75e5)

SESSION: 203 | April 7, 2026 | PAPER_859-877 integration
"""

import math
import json
import sys
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass

# ═══════════════════════════════════════════════════════════════════════════════
# §1  CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

G       = 6.67430e-11
c       = 2.99792e8
hbar    = 1.05457e-34
PI      = math.pi
M_sun   = 1.98892e30

# UQFF
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_UA  = 7.09e-36
RHO_SCM = 7.09e-37
ETA_AETHER = 1e-22

# Wormhole
RHO_PLUS_P_EXOTIC = -1.75e5   # kg/m³ — exotic matter density (PAPER_877)


# ═══════════════════════════════════════════════════════════════════════════════
# §2  BSFG METRIC FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class BSFGMetric:
    """
    BSFG wormhole metric parameters.

    The BSFG metric extends Morris-Thorne with:
    - B: Buoyancy-modified shape function via β_i and [SSq]
    - S: SCm redshift correction via H_SCm
    - F: Fluid density profile from ρ_SCm/ρ_UA ratio
    - G: Geometric SO(3)×U(1)²³ symmetry from 26D compactification
    """
    M_kg: float = 4.3e6 * M_sun        # central mass
    r0_m: float = 1e4                    # throat radius
    b0_m: float = 1e4                    # shape function scale
    A_trace: float = 4.04               # aether tensor trace


def redshift_function(r: float, metric: BSFGMetric) -> float:
    """
    Redshift function Φ(r) with aether correction:
      Φ(r) = -GM/(rc²) × (1 + η·A_trace)
    """
    if r == 0:
        return -1e10
    return -G * metric.M_kg / (r * c**2) * (1 + ETA_AETHER * metric.A_trace)


def shape_function(r: float, metric: BSFGMetric) -> float:
    """
    Buoyancy-modified shape function:
      b(r) = b₀ × (r₀/r) × (1 - β_i × [SSq])

    Must satisfy b(r) < r for traversability (outside throat).
    """
    if r == 0:
        return metric.b0_m
    buoyancy_mod = 1 - BETA_I * SSQ
    return metric.b0_m * (metric.r0_m / r) * buoyancy_mod


def metric_g_tt(r: float, metric: BSFGMetric) -> float:
    """g_tt = -exp(2Φ(r))"""
    Phi = redshift_function(r, metric)
    return -math.exp(2 * Phi)


def metric_g_rr(r: float, metric: BSFGMetric) -> float:
    """g_rr = 1 / (1 - b(r)/r)"""
    b = shape_function(r, metric)
    denom = 1 - b / max(r, 1e-10)
    if abs(denom) < 1e-30:
        return 1e30  # near throat
    return 1.0 / denom


def metric_g_angular(r: float, metric: BSFGMetric) -> float:
    """g_θθ = b² + r²"""
    b = shape_function(r, metric)
    return b**2 + r**2


def proper_distance(r: float, metric: BSFGMetric) -> float:
    """
    Proper radial distance: dl = dr / √(1 - b(r)/r)
    Integrand at radius r.
    """
    g_rr = metric_g_rr(r, metric)
    return math.sqrt(abs(g_rr))


# ═══════════════════════════════════════════════════════════════════════════════
# §3  EXOTIC MATTER & TRAVERSABILITY
# ═══════════════════════════════════════════════════════════════════════════════

def exotic_matter_density(r: float, metric: BSFGMetric) -> float:
    """
    Exotic matter density required for traversability:
      ρ_exotic = -(1/8πG) × b'(r)/r²

    PAPER_877: ρ + P = -1.75 × 10⁵ kg/m³ achievable via SCm vacuum.
    """
    dr = max(r * 1e-6, 1e-10)
    b_plus = shape_function(r + dr, metric)
    b_minus = shape_function(r - dr, metric)
    b_prime = (b_plus - b_minus) / (2 * dr)
    return -(1.0 / (8 * PI * G)) * b_prime / max(r**2, 1e-20)


def nec_violation(r: float, metric: BSFGMetric) -> float:
    """
    Null Energy Condition violation: ρ + P.
    Negative values = NEC violated = wormhole possible.
    """
    rho = exotic_matter_density(r, metric)
    # Pressure from SCm vacuum: P = -ρ_SCm c² + buoyancy correction
    P = -RHO_SCM * c**2 + BETA_I * SSQ * rho
    return rho + P / c**2


def traversability_check(metric: BSFGMetric) -> Dict:
    """Check BSFG wormhole traversability conditions."""
    r_throat = metric.r0_m
    r_test = [r_throat, 2 * r_throat, 5 * r_throat, 10 * r_throat]

    checks = []
    for r in r_test:
        b = shape_function(r, metric)
        Phi = redshift_function(r, metric)
        rho_exotic = exotic_matter_density(r, metric)
        nec = nec_violation(r, metric)

        checks.append({
            "r_m": r,
            "r/r0": r / r_throat,
            "b(r)": b,
            "b(r)/r": b / r,
            "Phi(r)": Phi,
            "rho_exotic": rho_exotic,
            "NEC_violation": nec,
            "traversable": b < r and nec < 0,
        })

    return {
        "throat_radius_m": r_throat,
        "mass_kg": metric.M_kg,
        "buoyancy_mod": 1 - BETA_I * SSQ,
        "checks": checks,
        "all_traversable": all(ch["traversable"] for ch in checks),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# §4  GEODESIC INTEGRATOR
# ═══════════════════════════════════════════════════════════════════════════════

class GeodesicTracer:
    """
    Traces null (photon) and timelike (massive particle) geodesics
    through the BSFG wormhole.

    Uses Euler integration of the geodesic equation:
      d²x^μ/dλ² + Γ^μ_{αβ} (dx^α/dλ)(dx^β/dλ) = 0
    """

    def __init__(self, metric: BSFGMetric = None):
        self.metric = metric or BSFGMetric()

    def christoffel_r_tt(self, r: float) -> float:
        """Γ^r_tt = -g^rr × (∂g_tt/∂r) / 2"""
        dr = max(r * 1e-6, 1e-10)
        dg_tt = (metric_g_tt(r + dr, self.metric) - metric_g_tt(r - dr, self.metric)) / (2 * dr)
        g_rr = metric_g_rr(r, self.metric)
        return -dg_tt / (2 * g_rr)

    def christoffel_r_rr(self, r: float) -> float:
        """Γ^r_rr = g^rr × (∂g_rr/∂r) / 2"""
        dr = max(r * 1e-6, 1e-10)
        dg_rr = (metric_g_rr(r + dr, self.metric) - metric_g_rr(r - dr, self.metric)) / (2 * dr)
        g_rr = metric_g_rr(r, self.metric)
        return dg_rr / (2 * g_rr)

    def effective_potential(self, r: float, L: float, mu: float = 0) -> float:
        """
        Effective potential for radial motion:
          V_eff(r) = (1 - b(r)/r) × [L²/(b²+r²) + μ²c²]

        mu=0 for null geodesics, mu=1 for timelike.
        """
        b = shape_function(r, self.metric)
        grr_inv = 1 - b / max(r, 1e-10)
        g_ang = b**2 + r**2
        return grr_inv * (L**2 / max(g_ang, 1e-10) + mu**2 * c**2)

    def trace_null(self, r_start: float, dr_dlambda_0: float,
                   L: float, n_steps: int = 1000,
                   d_lambda: float = 1e-3) -> List[Dict]:
        """
        Trace a null geodesic (photon path) through the wormhole.

        Uses simplified radial equation:
          (dr/dλ)² = E² - V_eff(r, L, mu=0)
        """
        trajectory = []
        r = r_start
        dr = dr_dlambda_0
        lam = 0.0

        for step in range(n_steps):
            b = shape_function(r, self.metric)
            Phi = redshift_function(r, self.metric)
            V_eff = self.effective_potential(r, L, mu=0)

            trajectory.append({
                "step": step,
                "lambda": lam,
                "r_m": r,
                "r/r0": r / self.metric.r0_m,
                "b(r)": b,
                "Phi(r)": Phi,
                "V_eff": V_eff,
                "dr_dlambda": dr,
            })

            # Geodesic acceleration
            Gamma_r = self.christoffel_r_tt(r) + self.christoffel_r_rr(r) * dr**2
            dr += -Gamma_r * d_lambda
            r += dr * d_lambda
            lam += d_lambda

            # Prevent negative radius or runaway
            if r < self.metric.r0_m * 0.1:
                r = self.metric.r0_m * 0.1
                dr = abs(dr)  # bounce at throat
            if r > self.metric.r0_m * 1e6:
                break

        return trajectory

    def trace_timelike(self, r_start: float, v_initial: float,
                       L: float, n_steps: int = 1000,
                       d_tau: float = 1e-3) -> List[Dict]:
        """
        Trace a timelike geodesic (massive particle) through the wormhole.
        """
        trajectory = []
        r = r_start
        dr = v_initial / c  # normalized
        tau = 0.0

        for step in range(n_steps):
            b = shape_function(r, self.metric)
            V_eff = self.effective_potential(r, L, mu=1)

            trajectory.append({
                "step": step,
                "tau": tau,
                "r_m": r,
                "r/r0": r / self.metric.r0_m,
                "b(r)": b,
                "V_eff": V_eff,
                "dr_dtau": dr,
            })

            Gamma_r = self.christoffel_r_tt(r) + self.christoffel_r_rr(r) * dr**2
            dr += -Gamma_r * d_tau
            r += dr * c * d_tau
            tau += d_tau

            if r < self.metric.r0_m * 0.1:
                r = self.metric.r0_m * 0.1
                dr = abs(dr)
            if r > self.metric.r0_m * 1e6:
                break

        return trajectory

    def embedding_diagram(self, r_min: float = None, r_max: float = None,
                          n_points: int = 200) -> List[Dict]:
        """
        Compute z(r) embedding diagram for the wormhole equatorial slice.
          dz/dr = ±1/√(r/b(r) - 1)
        """
        r_min = r_min or self.metric.r0_m * 1.01
        r_max = r_max or self.metric.r0_m * 20

        dr = (r_max - r_min) / n_points
        points = []
        z_val = 0.0

        for i in range(n_points):
            r = r_min + i * dr
            b = shape_function(r, self.metric)
            ratio = r / max(b, 1e-10) - 1
            if ratio > 0:
                dz_dr = 1.0 / math.sqrt(ratio)
                z_val += dz_dr * dr
            points.append({"r_m": r, "z_m": z_val, "b(r)": b})

        return points

    def print_report(self):
        """Print wormhole geodesic trace report."""
        print("=" * 78)
        print("BSFG WORMHOLE GEODESIC TRACER")
        print("=" * 78)
        print(f"  M       = {self.metric.M_kg:.3e} kg ({self.metric.M_kg/M_sun:.1f} M_Sun)")
        print(f"  r₀      = {self.metric.r0_m:.3e} m (throat)")
        print(f"  β_i     = {BETA_I}")
        print(f"  [SSq]   = {SSQ}")
        print(f"  b_mod   = 1 - β_i·[SSq] = {1 - BETA_I * SSQ:.4f}")
        print(f"  ρ+P     = {RHO_PLUS_P_EXOTIC:.2e} kg/m³ (exotic matter)")

        # Traversability
        trav = traversability_check(self.metric)
        print(f"\n▶ Traversability Check:")
        for ch in trav["checks"]:
            status = "PASS" if ch["traversable"] else "FAIL"
            print(f"    r/r₀={ch['r/r0']:.1f}: b/r={ch['b(r)/r']:.4f} NEC={ch['NEC_violation']:.4e} [{status}]")
        print(f"    Overall: {'TRAVERSABLE' if trav['all_traversable'] else 'NOT TRAVERSABLE'}")

        # Null geodesic
        print(f"\n▶ Null Geodesic Trace (first/last 5 steps):")
        null_traj = self.trace_null(
            r_start=self.metric.r0_m * 5,
            dr_dlambda_0=-1e3,
            L=1e10,
            n_steps=200
        )
        for t in null_traj[:5]:
            print(f"    λ={t['lambda']:.3f} r/r₀={t['r/r0']:.4f} dr/dλ={t['dr_dlambda']:.4e}")
        if len(null_traj) > 10:
            print(f"    ...")
            for t in null_traj[-5:]:
                print(f"    λ={t['lambda']:.3f} r/r₀={t['r/r0']:.4f} dr/dλ={t['dr_dlambda']:.4e}")
        print(f"    Total steps: {len(null_traj)}")

        # Embedding
        emb = self.embedding_diagram(n_points=50)
        print(f"\n▶ Embedding Diagram (first/last 5 points):")
        for p in emb[:5]:
            print(f"    r={p['r_m']:.3e} z={p['z_m']:.3e}")
        if len(emb) > 10:
            print(f"    ...")
            for p in emb[-5:]:
                print(f"    r={p['r_m']:.3e} z={p['z_m']:.3e}")

        print("=" * 78)


# ═══════════════════════════════════════════════════════════════════════════════
# §5  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    metric = BSFGMetric()
    tracer = GeodesicTracer(metric)

    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        result = {
            "metric": {
                "M_kg": metric.M_kg, "r0_m": metric.r0_m,
                "b0_m": metric.b0_m, "A_trace": metric.A_trace,
            },
            "traversability": traversability_check(metric),
            "null_geodesic": tracer.trace_null(metric.r0_m * 5, -1e3, 1e10, 200),
            "embedding": tracer.embedding_diagram(n_points=50),
        }
        outfile = sys.argv[2] if len(sys.argv) > 2 else "bsfg_wormhole_results.json"
        with open(outfile, "w") as f:
            json.dump(result, f, indent=2, default=str)
        print(f"Exported to {outfile}")
    else:
        tracer.print_report()


if __name__ == "__main__":
    main()
