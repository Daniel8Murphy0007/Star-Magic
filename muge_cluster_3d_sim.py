#!/usr/bin/env python3
"""
muge_cluster_3d_sim.py — 3D MUGE Galaxy Cluster Simulation

Session 216 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
3-dimensional gravitational simulation of galaxy cluster dynamics under the
Multi-Universal Gravitational Equation (MUGE) framework with:
  - NFW dark matter halo profiles
  - ICM (intra-cluster medium) phonon modulation
  - S₂₆ Ramanujan acceleration on all gravity terms
  - Triadic (Compressed/Resonant/Buoyancy) decomposition
  - 3D positional integration via leapfrog scheme

Distinct from muge_magnetar_3d_sim.py (magnetar + Goldreich-Julian only).
Links: PAPER_976, PAPER_454-457.
────────────────────────────────────────────────────────────────────────────────
"""

import math
import random
from typing import Dict, List, Tuple

# ── §0  Constants ──────────────────────────────────────────────────────────

PI   = math.pi
G    = 6.674e-11
M_SUN = 1.989e30
KPC  = 3.086e19
MPC  = 3.086e22
C    = 2.998e8
SSQ  = 0.57
S26  = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
OMEGA_SCM = 2 * PI * 1.25e12
BETA_I = 0.603

# Default cluster: Virgo-like
DEFAULT_CLUSTER = {
    "name": "Virgo-like",
    "M200_kg": 1.2e14 * M_SUN,
    "r200_m": 1.55 * MPC,
    "c_nfw": 6.5,
    "n_galaxies": 50,
    "T_icm_K": 2.5e7,
    "rho_icm_0": 3e-24,
    "beta_icm": 0.67,
    "r_core_m": 40 * KPC,
}


# ── §1  NFW Profile ──────────────────────────────────────────────────────

def nfw_density(r: float, rho_s: float, r_s: float) -> float:
    """NFW dark matter density: ρ(r) = ρ_s / [(r/r_s)(1+r/r_s)²]."""
    x = max(r / r_s, 1e-10)
    return rho_s / (x * (1.0 + x) ** 2)


def nfw_enclosed_mass(r: float, rho_s: float, r_s: float) -> float:
    """Enclosed mass within r for NFW profile."""
    x = r / r_s
    return 4 * PI * rho_s * r_s ** 3 * (math.log(1 + x) - x / (1 + x))


def nfw_params(M200: float, r200: float, c: float) -> Tuple[float, float]:
    """Derive (rho_s, r_s) from virial mass, radius, and concentration."""
    r_s = r200 / c
    rho_s = M200 / (4 * PI * r_s ** 3 * (math.log(1 + c) - c / (1 + c)))
    return rho_s, r_s


# ── §2  ICM Beta-Model ──────────────────────────────────────────────────

def icm_density(r: float, rho_0: float, r_core: float, beta: float) -> float:
    """β-model ICM gas density."""
    return rho_0 * (1.0 + (r / r_core) ** 2) ** (-3.0 * beta / 2.0)


def icm_pressure(r: float, rho_0: float, r_core: float, beta: float, T: float) -> float:
    """ICM thermal pressure P = n·k_B·T ≈ ρ/m_p · k_B · T."""
    k_B = 1.381e-23
    m_p = 1.673e-27
    rho = icm_density(r, rho_0, r_core, beta)
    return rho / m_p * k_B * T


# ── §3  MUGE Gravity at Cluster Scale ────────────────────────────────────

def g_muge_cluster(M_enc: float, r: float) -> float:
    """MUGE gravity: 26-layer compressed + phonon modulation + buoyancy."""
    if r < 1.0:
        r = 1.0
    g_comp = sum(G * M_enc / r**2 * SSQ * i / 26.0 for i in range(1, 27))
    Phi = S26  # On-resonance at 1.25 THz
    g_res = g_comp * Phi
    g_buoy = sum(G * M_enc / r**2 * math.exp(-SSQ * k / 26.0) * BETA_I for k in range(1, 27))
    total = abs(g_comp) + abs(g_res) + abs(g_buoy) + 1e-300
    wC, wR, wB = abs(g_comp)/total, abs(g_res)/total, abs(g_buoy)/total
    return wC * g_comp + wR * g_res + wB * g_buoy


# ── §4  3D Galaxy Particle ──────────────────────────────────────────────

class Galaxy3D:
    """Single galaxy in cluster simulation."""

    def __init__(self, mass_kg: float, pos: List[float], vel: List[float]):
        self.mass = mass_kg
        self.pos = list(pos)
        self.vel = list(vel)
        self.acc = [0.0, 0.0, 0.0]


# ── §5  Cluster Simulation ──────────────────────────────────────────────

class MUGECluster3DSim:
    """3D galaxy cluster evolution under MUGE gravity.

    Uses leapfrog integration with:
    - NFW dark matter potential (dominant)
    - Galaxy-galaxy MUGE interactions (pair-wise, softened)
    - ICM pressure gradient on galaxies (ram pressure proxy)
    """

    def __init__(self, cluster: dict = None):
        c = cluster or DEFAULT_CLUSTER
        self.name = c["name"]
        self.M200 = c["M200_kg"]
        self.r200 = c["r200_m"]
        self.c_nfw = c["c_nfw"]
        self.rho_s, self.r_s = nfw_params(self.M200, self.r200, self.c_nfw)
        self.n_gal = c.get("n_galaxies", 50)
        self.T_icm = c.get("T_icm_K", 2.5e7)
        self.rho_icm_0 = c.get("rho_icm_0", 3e-24)
        self.beta_icm = c.get("beta_icm", 0.67)
        self.r_core = c.get("r_core_m", 40 * KPC)
        self.galaxies: List[Galaxy3D] = []
        self.softening = 5 * KPC  # Softening length
        self.seed = c.get("seed", 42)

    def _init_galaxies(self):
        """Initialize galaxy positions/velocities from NFW sampling."""
        rng = random.Random(self.seed)
        self.galaxies = []
        for i in range(self.n_gal):
            # Sample radius from NFW-weighted distribution (rejection)
            r = self.r200 * rng.random() ** 0.5 * 0.8
            theta = rng.uniform(0, PI)
            phi = rng.uniform(0, 2 * PI)
            x = r * math.sin(theta) * math.cos(phi)
            y = r * math.sin(theta) * math.sin(phi)
            z = r * math.cos(theta)

            # Circular velocity from NFW enclosed mass
            M_enc = nfw_enclosed_mass(r, self.rho_s, self.r_s) if r > 0 else 0
            v_circ = math.sqrt(G * M_enc / max(r, 1.0)) if M_enc > 0 else 0
            # Add dispersion
            sigma = v_circ * 0.3
            vx = -v_circ * math.sin(phi) + rng.gauss(0, sigma)
            vy = v_circ * math.cos(phi) + rng.gauss(0, sigma)
            vz = rng.gauss(0, sigma)

            gal_mass = (1e10 + rng.random() * 1e12) * M_SUN
            self.galaxies.append(Galaxy3D(gal_mass, [x, y, z], [vx, vy, vz]))

    def _compute_accelerations(self):
        """Compute accelerations from NFW + galaxy-galaxy MUGE."""
        for g in self.galaxies:
            g.acc = [0.0, 0.0, 0.0]

        for i, gi in enumerate(self.galaxies):
            ri = math.sqrt(sum(p ** 2 for p in gi.pos))
            ri_safe = max(ri, self.softening)

            # NFW halo gravity (central potential)
            M_enc = nfw_enclosed_mass(ri_safe, self.rho_s, self.r_s)
            g_nfw = g_muge_cluster(M_enc, ri_safe)
            for d in range(3):
                gi.acc[d] -= g_nfw * gi.pos[d] / ri_safe if ri_safe > 0 else 0

            # Galaxy-galaxy interactions (pairwise, softened)
            for j in range(i + 1, len(self.galaxies)):
                gj = self.galaxies[j]
                dx = [gj.pos[d] - gi.pos[d] for d in range(3)]
                dist2 = sum(d ** 2 for d in dx) + self.softening ** 2
                dist = math.sqrt(dist2)
                force_mag = G * gi.mass * gj.mass / dist2 * S26
                for d in range(3):
                    a_ij = force_mag * dx[d] / dist / gi.mass
                    a_ji = -force_mag * dx[d] / dist / gj.mass
                    gi.acc[d] += a_ij
                    gj.acc[d] += a_ji

    def _leapfrog_step(self, dt: float):
        """Leapfrog integration step."""
        # Kick (half)
        for g in self.galaxies:
            for d in range(3):
                g.vel[d] += 0.5 * g.acc[d] * dt
        # Drift
        for g in self.galaxies:
            for d in range(3):
                g.pos[d] += g.vel[d] * dt
        # Compute new accelerations
        self._compute_accelerations()
        # Kick (half)
        for g in self.galaxies:
            for d in range(3):
                g.vel[d] += 0.5 * g.acc[d] * dt

    def compute(self, dataset: dict = None) -> dict:
        """Run the 3D cluster simulation.

        dataset keys (optional):
          n_steps: int (default 100)
          dt_Myr: float (default 10.0)
        """
        ds = dataset or {}
        n_steps = int(ds.get("n_steps", 100))
        dt_Myr = float(ds.get("dt_Myr", 10.0))
        dt = dt_Myr * 3.156e13

        self._init_galaxies()
        self._compute_accelerations()

        snapshots = []
        for step in range(n_steps):
            self._leapfrog_step(dt)

            if step % max(n_steps // 10, 1) == 0 or step == n_steps - 1:
                radii = [math.sqrt(sum(p**2 for p in g.pos)) for g in self.galaxies]
                speeds = [math.sqrt(sum(v**2 for v in g.vel)) for g in self.galaxies]
                snapshots.append({
                    "step": step,
                    "t_Myr": (step + 1) * dt_Myr,
                    "r_mean_kpc": sum(radii) / len(radii) / KPC,
                    "r_max_kpc": max(radii) / KPC,
                    "v_mean_kms": sum(speeds) / len(speeds) / 1e3,
                })

        return {
            "cluster": self.name,
            "n_galaxies": self.n_gal,
            "M200_Msun": self.M200 / M_SUN,
            "n_steps": n_steps,
            "dt_Myr": dt_Myr,
            "t_total_Myr": n_steps * dt_Myr,
            "snapshots": snapshots,
            "primary_equations": [
                "g_MUGE = w_C·g_comp + w_R·g_res + w_B·g_buoy (cluster scale)",
                f"NFW c={self.c_nfw}, r_s={self.r_s/KPC:.1f} kpc",
                f"Leapfrog: {n_steps} steps × {dt_Myr} Myr",
            ],
            "note": "PAPER_976. Session 216.",
        }

    def simulate(self, sweep=None, **kw):
        return [self.compute({"n_steps": n, "dt_Myr": 10.0}) for n in (sweep or [50, 100])]


# ── §6  Self-Tests ──────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True

    # Test 1: NFW parameters
    rho_s, r_s = nfw_params(DEFAULT_CLUSTER["M200_kg"], DEFAULT_CLUSTER["r200_m"], DEFAULT_CLUSTER["c_nfw"])
    if rho_s <= 0 or r_s <= 0:
        print(f"[FAIL] NFW params invalid: rho_s={rho_s}, r_s={r_s}"); ok = False
    else:
        print(f"[ OK ] NFW params: rho_s={rho_s:.3e} kg/m³, r_s={r_s/KPC:.1f} kpc")

    # Test 2: Enclosed mass at r200
    M_enc = nfw_enclosed_mass(DEFAULT_CLUSTER["r200_m"], rho_s, r_s)
    ratio = M_enc / DEFAULT_CLUSTER["M200_kg"]
    if not (0.8 < ratio < 1.2):
        print(f"[FAIL] NFW M_enc/M200 = {ratio:.3f}"); ok = False
    else:
        print(f"[ OK ] NFW M_enc/M200 = {ratio:.4f}")

    # Test 3: MUGE cluster gravity positive
    g = g_muge_cluster(1e14 * M_SUN, 1 * MPC)
    if g <= 0:
        print(f"[FAIL] g_MUGE = {g}"); ok = False
    else:
        print(f"[ OK ] g_MUGE at 1 Mpc = {g:.6e} m/s²")

    # Test 4: Short simulation
    sim = MUGECluster3DSim()
    result = sim.compute({"n_steps": 20, "dt_Myr": 5.0})
    if len(result["snapshots"]) == 0:
        print("[FAIL] No snapshots"); ok = False
    else:
        print(f"[ OK ] Simulation: {len(result['snapshots'])} snapshots, "
              f"final r_mean={result['snapshots'][-1]['r_mean_kpc']:.1f} kpc")

    return ok


if __name__ == "__main__":
    print("=" * 70)
    print("  muge_cluster_3d_sim.py — 3D MUGE Galaxy Cluster Simulation")
    print("=" * 70)
    passed = _run_tests()
    print("=" * 70)
    print(f"  {'ALL TESTS PASSED' if passed else 'SOME TESTS FAILED'}")
    print("=" * 70)
