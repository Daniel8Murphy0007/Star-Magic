#!/usr/bin/env python3
"""
muge_magnetar_3d_sim.py — 3D MUGE Magnetar Simulation Generator

Session 215 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Generates 3D MUGE magnetar simulation data: SCm core, magnetic vortex tubes,
and phonon resonance shells.

Layers:
  Core:   SCm superconducting interior with Δ(T) BCS gap
  Vortex: Abrikosov-type magnetic flux tubes at B > B_crit
  Shells: 26-state phonon resonance shells (HRes spectral ladder)

Output: 3D coordinate arrays for visualization.
────────────────────────────────────────────────────────────────────────────────
"""

import math
from typing import Dict, List, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
HBAR      = 1.055e-34
K_B       = 1.381e-23
C         = 2.998e8
G         = 6.674e-11
M_SUN     = 1.989e30
OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
S26       = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))
B_CRIT    = 4.4e13   # T, critical magnetic field for magnetars
R_NS      = 1.2e4    # m, typical neutron star radius


# ── §1  SCm Core Model ────────────────────────────────────────────────────

class SCmCoreModel:
    """SCm superconducting core with radial BCS gap profile.

    Δ(r) = Δ_0 · (1 - (r/R_core)²) for r < R_core, else 0
    Session 215.
    """

    def __init__(self, R_core: float = 0.5 * R_NS, Delta_0_eV: float = 1.0e-3):
        self.R_core = R_core
        self.Delta_0 = Delta_0_eV * 1.602e-19  # J

    def gap_profile(self, r: float) -> float:
        if r >= self.R_core:
            return 0.0
        return self.Delta_0 * (1 - (r / self.R_core) ** 2)

    def generate_shell(self, N_theta: int = 20, N_phi: int = 40) -> List[Dict]:
        """Generate core surface coordinates."""
        points = []
        for i in range(N_theta + 1):
            theta = PI * i / N_theta
            for j in range(N_phi):
                phi = 2 * PI * j / N_phi
                x = self.R_core * math.sin(theta) * math.cos(phi)
                y = self.R_core * math.sin(theta) * math.sin(phi)
                z = self.R_core * math.cos(theta)
                points.append({"x": x, "y": y, "z": z, "Delta_J": self.Delta_0})
        return points


# ── §2  Magnetic Vortex Tube Model ─────────────────────────────────────────

class MagneticVortexModel:
    """Abrikosov-type magnetic flux tubes for B > B_crit.

    Flux quantization: Φ_0 = h/(2e) ≈ 2.068e-15 Wb
    Vortex density: n_v = B / Φ_0
    Session 215.
    """

    PHI_0 = 2.068e-15  # Wb, flux quantum

    def __init__(self, B_surface: float = 1e14, R_star: float = R_NS):
        self.B = B_surface
        self.R_star = R_star

    def vortex_density(self) -> float:
        return self.B / self.PHI_0

    def generate_tubes(self, N_tubes: int = 12, N_z: int = 20) -> List[Dict]:
        """Generate vortex tube center-line coordinates along magnetic axis."""
        tubes = []
        for i in range(N_tubes):
            angle = 2 * PI * i / N_tubes
            r_tube = 0.3 * self.R_star
            for j in range(N_z + 1):
                z_frac = -1.0 + 2.0 * j / N_z
                x = r_tube * math.cos(angle)
                y = r_tube * math.sin(angle)
                z = z_frac * self.R_star
                B_local = self.B * math.exp(-abs(z_frac))
                tubes.append({
                    "tube": i, "x": x, "y": y, "z": z,
                    "B_T": B_local, "supercritical": B_local > B_CRIT,
                })
        return tubes


# ── §3  Phonon Resonance Shells ────────────────────────────────────────────

class PhononResonanceShells:
    """26-state HRes phonon resonance shells around the magnetar.

    Shell n at radius R_n = R_NS · (1 + 0.05n), with energy E_n from spectral ladder.
    Session 215.
    """

    E0 = HBAR * OMEGA_SCM

    def shell_radius(self, n: int, R_base: float = R_NS) -> float:
        return R_base * (1 + 0.05 * n)

    def shell_energy(self, n: int) -> float:
        return self.E0 * (2 * PI) ** (n / 3.0) * S26

    def generate_shells(self, N_shells: int = 26, N_pts: int = 30) -> List[Dict]:
        """Generate representative points on each resonance shell."""
        shells = []
        for n in range(1, N_shells + 1):
            R_n = self.shell_radius(n)
            E_n = self.shell_energy(n)
            for i in range(N_pts):
                theta = PI * (i + 0.5) / N_pts
                phi = 2 * PI * i * (1 + 5 ** 0.5) / 2  # golden angle
                x = R_n * math.sin(theta) * math.cos(phi)
                y = R_n * math.sin(theta) * math.sin(phi)
                z = R_n * math.cos(theta)
                shells.append({
                    "shell": n, "x": x, "y": y, "z": z,
                    "R_m": R_n, "E_eV": E_n / 1.602e-19,
                })
        return shells


# ── §4  Full 3D Simulation Generator ──────────────────────────────────────

class MUGEMagnetar3DSim:
    """Full 3D MUGE Magnetar simulation: core + vortex + phonon shells.

    Combines SCm core, magnetic vortex tubes, and 26-state phonon shells
    into a single 3D dataset suitable for rendering.
    Session 215.
    """

    def __init__(self, B_surface: float = 1e14, M_Msun: float = 1.4):
        self.core = SCmCoreModel()
        self.vortex = MagneticVortexModel(B_surface)
        self.shells = PhononResonanceShells()
        self.M = M_Msun * M_SUN

    def compute(self, dataset: dict) -> dict:
        B = float(dataset.get("B_surface", 1e14))
        N_tubes = int(dataset.get("N_tubes", 12))

        self.vortex.B = B
        core_pts = self.core.generate_shell(N_theta=10, N_phi=20)
        vortex_pts = self.vortex.generate_tubes(N_tubes, N_z=10)
        shell_pts = self.shells.generate_shells(N_shells=26, N_pts=10)

        return {
            "B_surface": B,
            "core_points": len(core_pts),
            "vortex_points": len(vortex_pts),
            "shell_points": len(shell_pts),
            "total_points": len(core_pts) + len(vortex_pts) + len(shell_pts),
            "vortex_density_m2": self.vortex.vortex_density(),
            "primary_equations": [
                "SCm Core: Δ(r) = Δ_0·(1-(r/R)²)",
                f"Vortex density: n_v = B/Φ_0 = {self.vortex.vortex_density():.3e} /m²",
                "Phonon shells: R_n = R_NS·(1+0.05n), E_n = E_0·(2π)^{n/3}·S₂₆",
                f"Total 3D points: {len(core_pts) + len(vortex_pts) + len(shell_pts)}",
            ],
        }


# ── §5  Self-Tests ─────────────────────────────────────────────────────────

def _run_tests() -> bool:
    ok = True

    # Test 1: Core gap profile
    core = SCmCoreModel()
    d0 = core.gap_profile(0)
    dR = core.gap_profile(core.R_core)
    if d0 <= 0 or dR != 0:
        print("[FAIL] Core gap profile invalid"); ok = False
    else:
        print(f"[ OK ] Core Δ(0) = {d0:.4e} J, Δ(R) = {dR}")

    # Test 2: Vortex density
    vm = MagneticVortexModel(1e14)
    nv = vm.vortex_density()
    if nv <= 0:
        print("[FAIL] Vortex density should be positive"); ok = False
    else:
        print(f"[ OK ] Vortex density = {nv:.3e} /m²")

    # Test 3: Phonon shells
    ps = PhononResonanceShells()
    R1 = ps.shell_radius(1)
    R26 = ps.shell_radius(26)
    if R26 <= R1:
        print("[FAIL] Shell radii should increase with n"); ok = False
    else:
        print(f"[ OK ] R_1 = {R1:.0f} m, R_26 = {R26:.0f} m")

    # Test 4: Full simulation
    sim = MUGEMagnetar3DSim()
    result = sim.compute({"B_surface": 1e14})
    if result["total_points"] <= 0:
        print("[FAIL] Should generate points"); ok = False
    else:
        print(f"[ OK ] Full sim: {result['total_points']} points")

    return ok


if __name__ == "__main__":
    import sys
    success = _run_tests()
    sys.exit(0 if success else 1)
