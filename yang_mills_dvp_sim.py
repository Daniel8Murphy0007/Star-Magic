#!/usr/bin/env python3
"""
yang_mills_dvp_sim.py — Yang-Mills Dynamical Gap Simulations with DVP Primes
═════════════════════════════════════════════════════════════════════════════

PURPOSE: Simulate Yang-Mills mass gap dynamics using DVP (Dipole Vortex Prime)
encoding as the discretization lattice. Maps SU(2) gauge sectors to UQFF Ug3
string rotation nodes, computes dynamical gap m_gap² from SCm parameters.

KEY EQUATION (PAPER_183 §3.2):
  m_gap² = 2γ × H_SCm(0) / v_SCm²

  where γ = confinement string tension, H_SCm = SCm metric (0.99),
  v_SCm = SCm velocity.

DVP PRIME LATTICE:
  Each DVP prime p_n defines a discrete gauge connection spacing.
  The Yang-Mills action S_YM = -(1/4) ∫ F^a_μν F_a^μν is discretized
  on a lattice with spacing a_n = l_Planck × p_n / 113 (normalized to H).

REFERENCES:
  - CondensedPhysics4.py L3445: YangMillsMassGapVacuumDensityEvolutionCalculator
  - CondensedPhysics4.py L8651: DipoleVortexPrimeEncodingCalculator
  - PAPER_183: Yang-Mills Hamiltonian SCm/UA Framework
  - PAPER_388: Yang-Mills VDS Evolution
  - PAPER_841 L165: Millennium Prize mass gap

SESSION: 203 | April 7, 2026 | PAPER_859-877 integration
"""

import math
import json
import sys
from typing import Dict, List, Tuple

# ═══════════════════════════════════════════════════════════════════════════════
# §1  CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

G       = 6.67430e-11
c       = 2.99792e8
hbar    = 1.05457e-34
k_B     = 1.38065e-23
PI      = math.pi
M_sun   = 1.98892e30

# UQFF
KAPPA   = 5.787e-9
SSQ     = 0.57
H_SCM   = 0.99
BETA_I  = 0.603
U_UA    = 1e-4
RHO_SCM = 7.09e-37
RHO_UA  = 7.09e-36

# Yang-Mills
LAMBDA_QCD = 0.2           # GeV — QCD confinement scale
SIGMA_STRING = 0.18        # GeV² — string tension (lattice QCD value)
G_YM = 1.0                 # Yang-Mills coupling (normalized)
L_PLANCK = 1.616e-35       # m

# DVP primes
DVP_PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113]


# ═══════════════════════════════════════════════════════════════════════════════
# §2  DVP LATTICE
# ═══════════════════════════════════════════════════════════════════════════════

class DVPLattice:
    """
    Discretizes Yang-Mills gauge field on a DVP prime lattice.
    Each site n has spacing a_n = l_Planck × p_n / 113.
    """

    def __init__(self, N_sites: int = 26):
        self.N = min(N_sites, len(DVP_PRIMES))
        self.primes = DVP_PRIMES[:self.N]

    def lattice_spacing(self, n: int) -> float:
        """Lattice spacing at site n (meters)."""
        return L_PLANCK * self.primes[n] / 113.0

    def plaquette_action(self, n: int, B_field: float) -> float:
        """
        Wilson plaquette action at site n:
          S_n = (g²/2) × a_n⁴ × B²
        """
        a = self.lattice_spacing(n)
        return (G_YM**2 / 2.0) * a**4 * B_field**2

    def total_action(self, B_field: float) -> float:
        """Total lattice action summed over all sites."""
        return sum(self.plaquette_action(n, B_field) for n in range(self.N))

    def effective_coupling(self, n: int, energy_scale_GeV: float = 1.0) -> float:
        """
        Running coupling at site n:
          g²(μ) = g²(Λ) / (1 + (11/16π²) g²(Λ) ln(μ/Λ))
        Asymptotic freedom for SU(2).
        """
        mu = energy_scale_GeV
        Lambda = LAMBDA_QCD
        if mu <= Lambda:
            return G_YM**2  # Confinement regime
        beta_0 = 11.0 / (16 * PI**2)
        return G_YM**2 / (1 + beta_0 * G_YM**2 * math.log(mu / Lambda))


# ═══════════════════════════════════════════════════════════════════════════════
# §3  YANG-MILLS MASS GAP ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

class YangMillsDVPGapSimulator:
    """
    Computes Yang-Mills mass gap using DVP prime lattice discretization.

    Gap equation (PAPER_183):
      m_gap² = 2σ × H_SCm / v_SCm²

    where σ = string tension, H_SCm = SCm metric, v_SCm = SCm velocity.
    """

    def __init__(self, N_sites: int = 26):
        self.lattice = DVPLattice(N_sites)

    def mass_gap_scm(self, v_scm: float = None) -> Dict:
        """
        Compute mass gap from SCm parameters.

        m_gap² = 2σ × H_SCm / v_SCm²
        """
        if v_scm is None:
            v_scm = U_UA * c  # SCm velocity from aether fraction
        m_gap_sq = 2 * SIGMA_STRING * H_SCM / (v_scm / (1.602e-10 / c**2))**2
        # Convert to natural units: σ in GeV², v in natural
        sigma_nat = SIGMA_STRING  # GeV²
        v_nat = v_scm / c  # dimensionless
        m_gap_sq_nat = 2 * sigma_nat * H_SCM / max(v_nat**2, 1e-30)
        m_gap_GeV = math.sqrt(abs(m_gap_sq_nat))

        return {
            "m_gap_squared_GeV2": m_gap_sq_nat,
            "m_gap_GeV": m_gap_GeV,
            "sigma_GeV2": sigma_nat,
            "H_SCm": H_SCM,
            "v_SCm_ms": v_scm,
            "v_SCm_natural": v_nat,
            "equation": "m_gap^2 = 2*sigma*H_SCm / v_SCm^2",
        }

    def dvp_gap_spectrum(self) -> List[Dict]:
        """
        Compute mass gap at each DVP prime lattice site.
        Each site has a different effective coupling → different gap.
        """
        spectrum = []
        for n in range(self.lattice.N):
            p = self.lattice.primes[n]
            a = self.lattice.lattice_spacing(n)
            g_eff = self.lattice.effective_coupling(n, energy_scale_GeV=1.0 / (a * 1e15))

            # Gap at this scale: m_n = σ^{1/2} × (p_n/113) × H_SCm
            m_n = math.sqrt(SIGMA_STRING) * (p / 113.0) * H_SCM

            spectrum.append({
                "site": n,
                "dvp_prime": p,
                "lattice_spacing_m": a,
                "g_effective": g_eff,
                "m_gap_GeV": m_n,
                "resonant": p > 26,
            })
        return spectrum

    def vacuum_density_evolution(self, t_steps: int = 100, t_max: float = 1e15) -> List[Dict]:
        """
        VDS evolution of Yang-Mills vacuum through DVP lattice.
        Tracks ρ_vac(t) = ρ_SCm × exp(-κt) × (1 + Σ_n S_n(B))
        """
        dt = t_max / t_steps
        evolution = []
        B = 1e-5  # background magnetic field

        for step in range(t_steps + 1):
            t = step * dt
            rho_base = RHO_SCM * math.exp(-KAPPA * t)
            S_total = self.lattice.total_action(B * math.exp(-KAPPA * t))
            rho_vac = rho_base * (1 + S_total)
            ssq_eff = rho_vac / max(RHO_UA, 1e-50)

            evolution.append({
                "step": step,
                "t_s": t,
                "rho_vac": rho_vac,
                "S_lattice": S_total,
                "ssq_effective": ssq_eff,
            })
        return evolution

    def print_report(self):
        """Print mass gap simulation report."""
        print("=" * 78)
        print("YANG-MILLS DYNAMICAL GAP — DVP PRIME LATTICE SIMULATION")
        print("=" * 78)

        gap = self.mass_gap_scm()
        print(f"\n▶ SCm Mass Gap (PAPER_183 §3.2):")
        print(f"    {gap['equation']}")
        print(f"    σ (string tension) = {gap['sigma_GeV2']:.3f} GeV²")
        print(f"    H_SCm              = {gap['H_SCm']}")
        print(f"    v_SCm              = {gap['v_SCm_ms']:.2e} m/s")
        print(f"    m_gap              = {gap['m_gap_GeV']:.4f} GeV")
        print(f"    m_gap²             = {gap['m_gap_squared_GeV2']:.4e} GeV²")

        print(f"\n{'─' * 78}")
        print("DVP PRIME SPECTRUM (gap at each lattice site)")
        print(f"{'─' * 78}")
        print(f"  {'Site':>4} {'Prime':>6} {'Spacing(m)':>12} {'g_eff':>10} {'m_gap(GeV)':>12} {'Resonant':>9}")
        for s in self.dvp_gap_spectrum():
            print(f"  {s['site']:>4} {s['dvp_prime']:>6} {s['lattice_spacing_m']:>12.3e} "
                  f"{s['g_effective']:>10.4f} {s['m_gap_GeV']:>12.6f} {'YES' if s['resonant'] else 'no':>9}")

        print(f"\n{'─' * 78}")
        print("VDS VACUUM EVOLUTION (first/last 5 steps)")
        print(f"{'─' * 78}")
        evo = self.vacuum_density_evolution(t_steps=50)
        for e in evo[:5] + [{"step": "...", "t_s": "...", "rho_vac": "...", "ssq_effective": "..."}] + evo[-5:]:
            if isinstance(e["step"], int):
                print(f"  step={e['step']:>3} t={e['t_s']:.3e}s ρ_vac={e['rho_vac']:.4e} [SSq]_eff={e['ssq_effective']:.4e}")
            else:
                print(f"  ...")

        print("=" * 78)


# ═══════════════════════════════════════════════════════════════════════════════
# §4  CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    sim = YangMillsDVPGapSimulator(N_sites=26)

    if len(sys.argv) > 1 and sys.argv[1] == "--json":
        result = {
            "mass_gap": sim.mass_gap_scm(),
            "dvp_spectrum": sim.dvp_gap_spectrum(),
            "vds_evolution": sim.vacuum_density_evolution(t_steps=50),
        }
        outfile = sys.argv[2] if len(sys.argv) > 2 else "yang_mills_dvp_results.json"
        with open(outfile, "w") as f:
            json.dump(result, f, indent=2, default=str)
        print(f"Exported to {outfile}")
    else:
        sim.print_report()


if __name__ == "__main__":
    main()
