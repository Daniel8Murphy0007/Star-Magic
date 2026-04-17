#!/usr/bin/env python3
"""
vds_dvp_bsh_expanded.py — Unified VDS/DVP/BSH Number Systems with Applications

Session 226 | Daniel Murphy
────────────────────────────────────────────────────────────────────────────────
Standalone module for the three UQFF number systems:
  - VDS (Vacuum Density Series)
  - DVP (Dipole Vortex Primes)
  - BSH (Buoyancy Harmonics)

with applications to AGN, NS (neutron stars), QGP, DM (dark matter), and CMB.

Gap closed:
  - Unified VDS/DVP/BSH engine
  - Cross-system correlation matrix
  - 5 astrophysical application domains

Physics:
  VDS: a_n^{VDS} = [SSq]^n / n^26        (vacuum spectral weight per mode n)
  DVP: a(p)^{DVP} = [SSq]^π(p) / p^26    (prime-indexed vortex amplitudes)
  BSH: h_k^{BSH} = β_i^k · H_SCm^{26-k}  (buoyancy harmonic modes)

ARCHITECTURE: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
SOLAR_MASS = 1.989e30
eV        = 1.602e-19

SSQ       = 0.57
BETA_I    = 0.603
H_SCM     = 0.99
F_UBI_RATIO = 0.6
OMEGA_SCM = 2 * PI * 1.25e12


def _is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def _prime_counting(p: int) -> int:
    """π(p): number of primes ≤ p."""
    return sum(1 for k in range(2, p + 1) if _is_prime(k))


# ── §1  VACUUM DENSITY SERIES (VDS) ───────────────────────────────────────

class VacuumDensitySeries:
    """VDS: a_n = [SSq]^n / n^26.

    The VDS encodes the spectral weight of the vacuum condensate at each
    harmonic mode n. The series converges super-exponentially due to
    [SSq] < 1 and the 26th-power suppression.
    """

    def __init__(self, ssq: float = SSQ, n_max: int = 50):
        self.ssq = ssq
        self.n_max = n_max

    def amplitude(self, n: int) -> float:
        """a_n^{VDS} = [SSq]^n / n^26."""
        return self.ssq ** n / n ** 26

    def partial_sum(self, n_terms: int = None) -> float:
        """S_VDS(N) = Σ_{n=1}^{N} a_n."""
        N = n_terms or self.n_max
        return sum(self.amplitude(n) for n in range(1, N + 1))

    def convergence_rate(self) -> float:
        """Ratio a_{n+1}/a_n → [SSq]·(n/(n+1))^26 ≈ [SSq] for large n."""
        return self.ssq

    def spectrum(self, n_max: int = None) -> list:
        """Full spectrum of VDS amplitudes."""
        N = n_max or self.n_max
        return [{"n": n, "a_n": self.amplitude(n), "log10_a_n": math.log10(self.amplitude(n)) if self.amplitude(n) > 0 else -999}
                for n in range(1, N + 1)]

    def compute(self, dataset: dict) -> dict:
        ssq = float(dataset.get("SSq", self.ssq))
        self.ssq = ssq
        n_max = int(dataset.get("n_max", self.n_max))
        S = self.partial_sum(n_max)
        a1 = self.amplitude(1)
        spectrum = self.spectrum(min(n_max, 30))
        return {
            "S_VDS": S,
            "a_1": a1,
            "n_max": n_max,
            "convergence_rate": self.convergence_rate(),
            "spectrum_len": len(spectrum),
            "primary_equations": [
                f"a_n^{{VDS}} = [SSq]^n / n^26",
                f"a_1 = {a1:.6e}",
                f"S_VDS({n_max}) = {S:.6e}",
                f"Convergence ratio = {self.convergence_rate():.4f}",
            ],
        }


# ── §2  DIPOLE VORTEX PRIMES (DVP) ────────────────────────────────────────

class DipoleVortexPrimes:
    """DVP: a(p) = [SSq]^π(p) / p^26.

    Prime-indexed vortex amplitudes. The π(p) exponent means higher
    primes contribute less, giving super-exponential convergence.
    Only primes p > 26 contribute to the vorticity bound.
    """

    def __init__(self, ssq: float = SSQ, p_max: int = 200):
        self.ssq = ssq
        self.p_max = p_max
        self.primes = [p for p in range(2, p_max + 1) if _is_prime(p)]

    def amplitude(self, p: int) -> float:
        """a(p)^{DVP} = [SSq]^π(p) / p^26."""
        pi_p = _prime_counting(p)
        return self.ssq ** pi_p / p ** 26

    def primes_above_26(self) -> list:
        """Primes contributing to vorticity bound (p > 26)."""
        return [p for p in self.primes if p > 26]

    def partial_sum(self, min_p: int = 2) -> float:
        """Σ_{p prime, p ≥ min_p} a(p)."""
        return sum(self.amplitude(p) for p in self.primes if p >= min_p)

    def vorticity_bound_sum(self) -> float:
        """Sum over primes > 26 only (relevant for enstrophy bound)."""
        return sum(self.amplitude(p) for p in self.primes_above_26())

    def compute(self, dataset: dict) -> dict:
        ssq = float(dataset.get("SSq", self.ssq))
        self.ssq = ssq
        S_all = self.partial_sum(2)
        S_gt26 = self.vorticity_bound_sum()
        n_primes = len(self.primes)
        n_gt26 = len(self.primes_above_26())
        return {
            "S_DVP_all": S_all,
            "S_DVP_gt26": S_gt26,
            "n_primes": n_primes,
            "n_primes_gt26": n_gt26,
            "primary_equations": [
                f"a(p)^{{DVP}} = [SSq]^π(p) / p^26",
                f"Σ a(p) [all primes ≤ {self.p_max}] = {S_all:.6e}",
                f"Σ a(p) [primes > 26] = {S_gt26:.6e}",
                f"N primes > 26: {n_gt26}",
            ],
        }


# ── §3  BUOYANCY HARMONICS (BSH) ──────────────────────────────────────────

class BuoyancyHarmonics:
    """BSH: h_k = β_i^k · H_SCm^{26-k}.

    Buoyancy harmonic modes interpolating between the buoyancy index β_i
    and the SCm coherence factor H_SCm across 26 dimensions.
    """

    def __init__(self, beta_i: float = BETA_I, h_scm: float = H_SCM):
        self.beta_i = beta_i
        self.h_scm = h_scm

    def harmonic(self, k: int) -> float:
        """h_k = β_i^k · H_SCm^{26-k}."""
        return self.beta_i ** k * self.h_scm ** (26 - k)

    def spectrum(self, k_max: int = 26) -> list:
        """Full harmonic spectrum k = 0..26."""
        return [{"k": k, "h_k": self.harmonic(k)} for k in range(k_max + 1)]

    def total_weight(self, k_max: int = 26) -> float:
        """W_BSH = Σ_{k=0}^{26} h_k."""
        return sum(self.harmonic(k) for k in range(k_max + 1))

    def dominant_mode(self) -> int:
        """Mode k* that maximises h_k.

        d/dk [k·ln(β_i) + (26-k)·ln(H_SCm)] = 0
        k* = (ln(H_SCm) - ln(β_i)) · 26 / (ln(H_SCm) - ln(β_i))
        Since β_i < H_SCm, h_k is monotonically decreasing → k* = 0.
        """
        if self.beta_i >= self.h_scm:
            return 26
        return 0  # h_0 = H_SCm^26 is the maximum

    def compute(self, dataset: dict) -> dict:
        beta = float(dataset.get("beta_i", self.beta_i))
        hscm = float(dataset.get("H_SCm", self.h_scm))
        self.beta_i = beta
        self.h_scm = hscm
        W = self.total_weight()
        h0 = self.harmonic(0)
        h26 = self.harmonic(26)
        return {
            "W_BSH": W,
            "h_0": h0,
            "h_26": h26,
            "h_0_over_h_26": h0 / h26 if h26 > 0 else float('inf'),
            "dominant_mode": self.dominant_mode(),
            "primary_equations": [
                f"h_k = β_i^k · H_SCm^(26-k)",
                f"h_0 = H_SCm^26 = {h0:.6e}",
                f"h_26 = β_i^26 = {h26:.6e}",
                f"W_BSH = Σ h_k = {W:.6e}",
            ],
        }


# ── §4  CROSS-SYSTEM CORRELATION ──────────────────────────────────────────

class VDS_DVP_BSH_Correlator:
    """Unified correlator across VDS, DVP, and BSH number systems.

    Computes cross-system products and ratios useful for astrophysical
    applications (AGN, NS, QGP, DM, CMB).
    """

    def __init__(self):
        self.vds = VacuumDensitySeries()
        self.dvp = DipoleVortexPrimes()
        self.bsh = BuoyancyHarmonics()

    def cross_product(self) -> float:
        """S_VDS · S_DVP · W_BSH — combined spectral weight."""
        return self.vds.partial_sum() * self.dvp.partial_sum() * self.bsh.total_weight()

    def agn_application(self, dataset: dict) -> dict:
        """AGN feedback: uses VDS for vacuum energy + BSH for jet structure."""
        M_bh = float(dataset.get("M_BH", 1e8 * SOLAR_MASS))
        L_edd = 4 * PI * G * M_bh * M_P * C / (6.652e-29)  # σ_T
        S_vds = self.vds.partial_sum()
        W_bsh = self.bsh.total_weight()
        f_feedback = S_vds * W_bsh * F_UBI_RATIO
        return {
            "domain": "AGN",
            "L_Eddington_W": L_edd,
            "f_feedback": f_feedback,
            "equation": f"f_AGN = S_VDS·W_BSH·(F_UBi/F_U) = {f_feedback:.6e}",
        }

    def ns_application(self, dataset: dict) -> dict:
        """Neutron star: DVP vortex modes for superfluid interior."""
        P_ms = float(dataset.get("P_ms", 33.0))
        omega_ns = 2 * PI / (P_ms * 1e-3)
        S_dvp = self.dvp.partial_sum()
        n_vortices = S_dvp * omega_ns / PI  # rough estimate
        return {
            "domain": "NS",
            "omega_NS_rad_s": omega_ns,
            "S_DVP": S_dvp,
            "n_vortex_estimate": n_vortices,
            "equation": f"n_vortex ∝ S_DVP·Ω_NS/π = {n_vortices:.6e}",
        }

    def qgp_application(self, dataset: dict) -> dict:
        """QGP: BSH modes as thermal harmonics."""
        T_MeV = float(dataset.get("T_MeV", 300))
        W = self.bsh.total_weight()
        eta_s = W / (4 * PI)
        return {
            "domain": "QGP",
            "T_MeV": T_MeV,
            "W_BSH": W,
            "eta_over_s": eta_s,
            "equation": f"η/s = W_BSH/(4π) = {eta_s:.6e}",
        }

    def dm_application(self, dataset: dict) -> dict:
        """Dark matter: VDS vacuum energy contribution to ρ_DM."""
        S_vds = self.vds.partial_sum()
        rho_vac = HBAR * OMEGA_SCM * S_vds / C ** 2
        return {
            "domain": "DM",
            "rho_vac_equivalence": rho_vac,
            "S_VDS": S_vds,
            "equation": f"ρ_DM_UQFF = ℏω_SCm·S_VDS/c² = {rho_vac:.6e} kg/m³",
        }

    def cmb_application(self, dataset: dict) -> dict:
        """CMB: BSH + VDS spectral imprint on angular power spectrum."""
        l_max = int(dataset.get("l_max", 2500))
        S_vds = self.vds.partial_sum()
        W_bsh = self.bsh.total_weight()
        delta_Cl = S_vds * W_bsh * 1e-10  # fractional perturbation
        return {
            "domain": "CMB",
            "l_max": l_max,
            "delta_Cl_fractional": delta_Cl,
            "equation": f"δC_l/C_l = S_VDS·W_BSH·10⁻¹⁰ = {delta_Cl:.6e}",
        }

    def compute(self, dataset: dict) -> dict:
        """Full unified VDS/DVP/BSH computation with astrophysical applications."""
        vds_r = self.vds.compute(dataset)
        dvp_r = self.dvp.compute(dataset)
        bsh_r = self.bsh.compute(dataset)

        apps = {
            "AGN": self.agn_application(dataset),
            "NS": self.ns_application(dataset),
            "QGP": self.qgp_application(dataset),
            "DM": self.dm_application(dataset),
            "CMB": self.cmb_application(dataset),
        }

        cross = self.cross_product()

        eqs = vds_r["primary_equations"] + dvp_r["primary_equations"] + bsh_r["primary_equations"]
        eqs.append(f"S_VDS · S_DVP · W_BSH = {cross:.6e}")
        for name, app in apps.items():
            eqs.append(f"  [{name}] {app['equation']}")

        return {
            "VDS": vds_r,
            "DVP": dvp_r,
            "BSH": bsh_r,
            "cross_product": cross,
            "applications": apps,
            "primary_equations": eqs,
        }


# ── §5  RUNNER ─────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("=" * 72)
    print("Unified VDS / DVP / BSH Number Systems with Applications")
    print("=" * 72)

    correlator = VDS_DVP_BSH_Correlator()
    result = correlator.compute({"T_MeV": 300, "P_ms": 33.0})

    for eq in result["primary_equations"]:
        print(f"  {eq}")

    print(f"\n5 astrophysical domains applied")
    print("\n✓ VDS/DVP/BSH expanded module complete")
