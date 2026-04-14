#!/usr/bin/env python3
"""
muge_3d_volumetric.py — 3D Volumetric MUGE Field Generator

Session 224 | Star Magic UQFF Framework
────────────────────────────────────────────────────────────────────────────────
Generates 3D volumetric MUGE gravitational field cubes on meshgrids.
Extends CondensedPhysics4.py MUGECluster3DSimCalc (radial-only output)
to full 3D density/velocity/gravity tensor field export.

Core equations:
    ρ_NFW(r)  = ρ_s / [(r/r_s)(1 + r/r_s)²]            [NFW dark matter halo]
    g_MUGE(r) = G·M(<r)/r² + Σ corrections               [8-term MUGE]
    v_circ(r) = sqrt(r · g_MUGE(r))                      [circular velocity]

Corrections: Expansion, SuperFreq, Envelope, Ug_sum, Cosmological,
             Quantum, Fluid (Navier-Stokes), Perturbation (dark matter)

Multi-system support: NGC 6302, Orion, Lagoon, Saturn, Crab, Andromeda, Sombrero

Architecture: Pure calculator. No hardcoded systems. Tier 2 compute.
────────────────────────────────────────────────────────────────────────────────
"""

import math
import time
from typing import Any, Dict, List, Optional, Tuple

# ── §0  CONSTANTS ──────────────────────────────────────────────────────────

PI        = math.pi
C         = 2.998e8
HBAR      = 1.055e-34
K_B       = 1.381e-23
G         = 6.674e-11
M_SUN     = 1.989e30
KPC       = 3.086e19
MPC       = 3.086e22

OMEGA_SCM = 2 * PI * 1.25e12
SSQ       = 0.57
BETA_I    = 0.603
KAPPA     = 0.0005 / 86400.0
H_0       = 2.195e-18          # s⁻¹  Hubble constant (67.4 km/s/Mpc)
LAMBDA_C  = 1.11e-52           # m⁻²  cosmological constant
S26_STATIC = sum(math.exp(-SSQ * k / 26.0) for k in range(1, 27))


# ── §1  NFW DARK MATTER PROFILE ───────────────────────────────────────────

class NFWProfile:
    """
    Navarro-Frenk-White dark matter density profile.

    ρ(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]

    Variables (from dataset):
        rho_s:  characteristic density (kg/m³)
        r_s:    scale radius (m)
        r:      radius (m) or array of radii
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        rho_s = dataset.get('rho_s', 1e-21)   # kg/m³
        r_s   = dataset.get('r_s', 20 * KPC)  # m
        r     = dataset.get('r', 10 * KPC)     # m

        if isinstance(r, (list, tuple)):
            densities = []
            for ri in r:
                x = max(ri / r_s, 1e-10)
                rho = rho_s / (x * (1 + x) ** 2)
                densities.append(rho)
            return {
                'primary_equations': ["ρ_NFW(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]"],
                'densities_kg_m3': densities,
                'rho_s': rho_s,
                'r_s_m': r_s,
            }
        else:
            x = max(r / r_s, 1e-10)
            rho = rho_s / (x * (1 + x) ** 2)
            return {
                'primary_equations': ["ρ_NFW(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]"],
                'density_kg_m3': rho,
                'rho_s': rho_s,
                'r_s_m': r_s,
                'r_m': r,
            }


# ── §2  EIGHT-TERM MUGE GRAVITY ───────────────────────────────────────────

class MUGEGravityField:
    """
    Eight-term MUGE gravitational acceleration.

    g_MUGE(r) = g_Newton + g_Expansion + g_Super + g_Envelope
              + g_Ug_sum + g_Cosmological + g_Quantum + g_Fluid

    Variables (from dataset):
        M:            enclosed mass (kg)
        r:            radius (m)
        Omega_rot:    rotation rate (rad/s, default 0)
        B_field:      magnetic field (T, default 0)
        rho_local:    local density (kg/m³, default 1e-21)
        v_bulk:       bulk velocity (m/s, default 0)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        M       = dataset.get('M', M_SUN)
        r       = max(dataset.get('r', KPC), 1.0)
        Omega   = dataset.get('Omega_rot', 0.0)
        B       = dataset.get('B_field', 0.0)
        rho     = dataset.get('rho_local', 1e-21)
        v_bulk  = dataset.get('v_bulk', 0.0)

        # 1. Newtonian base
        g_N = G * M / r ** 2

        # 2. Hubble expansion correction
        g_exp = -H_0 ** 2 * r

        # 3. SuperFreq magnetic suppression
        mu_0 = 4 * PI * 1e-7
        g_super = -B ** 2 / (2 * mu_0 * rho * r) if rho > 0 else 0.0

        # 4. Envelope (rotation)
        g_env = Omega ** 2 * r

        # 5. Ug_sum (26-layer buoyancy)
        g_ug = sum(
            G * M / r ** 2 * SSQ * i / 26 * BETA_I
            for i in range(1, 27)
        )

        # 6. Cosmological constant
        g_cosm = -LAMBDA_C * C ** 2 * r / 3

        # 7. Quantum correction
        g_quant = HBAR / (M * r ** 2) if M > 0 else 0.0

        # 8. Fluid (Navier-Stokes viscous)
        nu = 1e-6  # kinematic viscosity (m²/s, ISM)
        g_fluid = -nu * v_bulk / r ** 2 if r > 0 else 0.0

        g_total = g_N + g_exp + g_super + g_env + g_ug + g_cosm + g_quant + g_fluid

        return {
            'primary_equations': [
                "g_MUGE = g_N + g_exp + g_super + g_env + g_Ug + g_cosm + g_quant + g_fluid",
                f"g_Newton    = {g_N:.6e} m/s²",
                f"g_Expansion = {g_exp:.6e} m/s²",
                f"g_SuperFreq = {g_super:.6e} m/s²",
                f"g_Envelope  = {g_env:.6e} m/s²",
                f"g_Ug_sum    = {g_ug:.6e} m/s²",
                f"g_Cosm      = {g_cosm:.6e} m/s²",
                f"g_Quantum   = {g_quant:.6e} m/s²",
                f"g_Fluid     = {g_fluid:.6e} m/s²",
                f"g_TOTAL     = {g_total:.6e} m/s²",
            ],
            'g_total_m_s2': g_total,
            'g_newton': g_N,
            'g_expansion': g_exp,
            'g_super': g_super,
            'g_envelope': g_env,
            'g_ug_sum': g_ug,
            'g_cosmological': g_cosm,
            'g_quantum': g_quant,
            'g_fluid': g_fluid,
            'components': {
                'newton': g_N, 'expansion': g_exp, 'super': g_super,
                'envelope': g_env, 'ug_sum': g_ug, 'cosmological': g_cosm,
                'quantum': g_quant, 'fluid': g_fluid,
            },
        }


# ── §3  3D VOLUMETRIC FIELD GENERATOR ─────────────────────────────────────

class MUGEVolumetricField:
    """
    Generate 3D volumetric MUGE gravity/density/velocity fields on a meshgrid.

    Creates a cube of (nx, ny, nz) points centered on the system,
    computing at each point:
        - ρ_NFW(r)      dark matter density
        - g_MUGE(r)     gravitational acceleration magnitude
        - g_vec(x,y,z)  gravitational acceleration vector (radial)
        - v_circ(r)     circular velocity

    Variables (from dataset):
        M:          total enclosed mass (kg)
        rho_s:      NFW characteristic density (kg/m³)
        r_s:        NFW scale radius (m)
        box_size:   half-side of cube (m, default 50 kpc)
        nx, ny, nz: grid resolution (default 32 each)
        B_field:    magnetic field (T)
        Omega_rot:  rotation rate (rad/s)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        M       = dataset.get('M', 1e11 * M_SUN)
        rho_s   = dataset.get('rho_s', 1e-21)
        r_s     = dataset.get('r_s', 20 * KPC)
        box     = dataset.get('box_size', 50 * KPC)
        nx      = dataset.get('nx', 32)
        ny      = dataset.get('ny', 32)
        nz      = dataset.get('nz', 32)
        B       = dataset.get('B_field', 0.0)
        Omega   = dataset.get('Omega_rot', 0.0)

        nfw = NFWProfile()
        muge = MUGEGravityField()

        t_start = time.perf_counter()

        # Generate grid coordinates
        xs = [(-box + 2 * box * i / max(nx - 1, 1)) for i in range(nx)]
        ys = [(-box + 2 * box * j / max(ny - 1, 1)) for j in range(ny)]
        zs = [(-box + 2 * box * k / max(nz - 1, 1)) for k in range(nz)]

        # Output cubes (flat lists for JSON export)
        density_cube = []
        g_mag_cube = []
        gx_cube = []
        gy_cube = []
        gz_cube = []
        v_circ_cube = []

        r_min_used = float('inf')
        r_max_used = 0.0

        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    x = xs[ix]
                    y = ys[iy]
                    z = zs[iz]
                    r = math.sqrt(x ** 2 + y ** 2 + z ** 2)
                    r = max(r, 1.0)  # avoid singularity

                    r_min_used = min(r_min_used, r)
                    r_max_used = max(r_max_used, r)

                    # NFW density
                    nfw_res = nfw.compute({'rho_s': rho_s, 'r_s': r_s, 'r': r})
                    rho_val = nfw_res['density_kg_m3']

                    # MUGE gravity
                    muge_res = muge.compute({
                        'M': M, 'r': r,
                        'B_field': B, 'Omega_rot': Omega,
                        'rho_local': rho_val,
                    })
                    g_mag = muge_res['g_total_m_s2']

                    # Direction: radial inward (toward center)
                    inv_r = 1.0 / r
                    gx = -g_mag * x * inv_r
                    gy = -g_mag * y * inv_r
                    gz = -g_mag * z * inv_r

                    # Circular velocity
                    v_c = math.sqrt(abs(r * g_mag)) if g_mag > 0 else 0.0

                    density_cube.append(rho_val)
                    g_mag_cube.append(g_mag)
                    gx_cube.append(gx)
                    gy_cube.append(gy)
                    gz_cube.append(gz)
                    v_circ_cube.append(v_c)

        elapsed_ms = (time.perf_counter() - t_start) * 1000
        n_points = nx * ny * nz

        return {
            'primary_equations': [
                "ρ_NFW(r) = ρ_s / [(r/r_s)(1 + r/r_s)²]",
                "g_MUGE(r) = g_N + 7 corrections (8-term)",
                "v_circ(r) = sqrt(r · |g_MUGE(r)|)",
                f"Grid: {nx}×{ny}×{nz} = {n_points} points",
                f"Box: ±{box:.3e} m ({box/KPC:.1f} kpc)",
            ],
            'density_cube': density_cube,
            'g_magnitude_cube': g_mag_cube,
            'gx_cube': gx_cube,
            'gy_cube': gy_cube,
            'gz_cube': gz_cube,
            'v_circ_cube': v_circ_cube,
            'grid_shape': (nx, ny, nz),
            'box_size_m': box,
            'box_size_kpc': box / KPC,
            'n_points': n_points,
            'r_min_m': r_min_used,
            'r_max_m': r_max_used,
            'elapsed_ms': elapsed_ms,
        }


# ── §4  MULTI-SYSTEM VOLUMETRIC BATCH ─────────────────────────────────────

class MUGEMultiSystemVolume:
    """
    Batch volumetric MUGE fields across multiple astrophysical systems.

    Variables (from dataset):
        systems: list of dicts, each with:
            name:     system name
            M:        total mass (kg)
            rho_s:    NFW density (kg/m³)
            r_s:      NFW scale radius (m)
            box_size: half-cube side (m)
            nx, ny, nz: resolution
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        systems = dataset.get('systems', [])
        if not systems:
            return {'primary_equations': [], 'error': 'No systems provided'}

        engine = MUGEVolumetricField()
        results = []

        for sys in systems:
            res = engine.compute(sys)
            res['system_name'] = sys.get('name', 'unnamed')
            # Summary stats
            n = len(res['density_cube'])
            if n > 0:
                res['density_min'] = min(res['density_cube'])
                res['density_max'] = max(res['density_cube'])
                res['g_mag_min'] = min(res['g_magnitude_cube'])
                res['g_mag_max'] = max(res['g_magnitude_cube'])
                res['v_circ_max'] = max(res['v_circ_cube'])
            # Remove heavy cubes from batch summary (keep for individual access)
            summary = {k: v for k, v in res.items()
                       if k not in ('density_cube', 'g_magnitude_cube',
                                    'gx_cube', 'gy_cube', 'gz_cube', 'v_circ_cube')}
            results.append(summary)

        return {
            'primary_equations': [
                f"Multi-system 3D MUGE volumetric fields: {len(results)} systems",
                "Each system: NFW + 8-term MUGE + v_circ on 3D meshgrid",
            ],
            'system_summaries': results,
            'n_systems': len(results),
        }


# ── §5  ROTATION CURVE EXTRACTOR ──────────────────────────────────────────

class MUGERotationCurve:
    """
    Extract MUGE rotation curve v_circ(r) along the midplane.

    Variables (from dataset):
        M:        total mass (kg)
        rho_s:    NFW density (kg/m³)
        r_s:      NFW scale radius (m)
        r_min:    inner radius (m, default 1 kpc)
        r_max:    outer radius (m, default 100 kpc)
        n_radii:  number of radial bins (default 200)
    """

    def compute(self, dataset: dict) -> Dict[str, Any]:
        M     = dataset.get('M', 1e11 * M_SUN)
        rho_s = dataset.get('rho_s', 1e-21)
        r_s   = dataset.get('r_s', 20 * KPC)
        r_min = dataset.get('r_min', 1 * KPC)
        r_max = dataset.get('r_max', 100 * KPC)
        n_r   = dataset.get('n_radii', 200)

        muge = MUGEGravityField()

        radii = []
        v_circs = []
        g_totals = []
        g_newtons = []

        for i in range(n_r):
            r = r_min + (r_max - r_min) * i / max(n_r - 1, 1)
            radii.append(r)

            res = muge.compute({'M': M, 'r': r})
            g = res['g_total_m_s2']
            g_n = res['g_newton']

            v_c = math.sqrt(abs(r * g)) if g > 0 else 0.0
            v_n = math.sqrt(abs(r * g_n)) if g_n > 0 else 0.0

            v_circs.append(v_c)
            g_totals.append(g)
            g_newtons.append(g_n)

        # Flat rotation curve indicator: ratio of max to min v_circ in outer half
        outer_v = v_circs[n_r // 2:]
        v_max_outer = max(outer_v) if outer_v else 0
        v_min_outer = min(outer_v) if outer_v else 0
        flatness = v_min_outer / v_max_outer if v_max_outer > 0 else 0

        return {
            'primary_equations': [
                "v_circ(r) = sqrt(r · g_MUGE(r))",
                f"Radial range: {r_min/KPC:.1f} – {r_max/KPC:.1f} kpc",
                f"Flatness (v_min/v_max outer): {flatness:.4f}",
            ],
            'radii_m': radii,
            'radii_kpc': [r / KPC for r in radii],
            'v_circ_m_s': v_circs,
            'v_circ_km_s': [v / 1e3 for v in v_circs],
            'g_total_m_s2': g_totals,
            'g_newton_m_s2': g_newtons,
            'flatness_ratio': flatness,
            'n_radii': n_r,
        }


# ══════════════════════════════════════════════════════════════════════════════
# §6  SELF-TESTS
# ══════════════════════════════════════════════════════════════════════════════

def _run_tests():
    print("=" * 72)
    print("muge_3d_volumetric.py — Self-Tests")
    print("=" * 72)

    passed = 0
    failed = 0

    # Test 1: NFW profile
    nfw = NFWProfile()
    res = nfw.compute({'rho_s': 1e-21, 'r_s': 20 * KPC, 'r': 20 * KPC})
    expected_rho = 1e-21 / (1.0 * (1 + 1.0) ** 2)
    assert abs(res['density_kg_m3'] - expected_rho) < 1e-30
    print(f"  [PASS] Test 1: NFW ρ(r_s) = {res['density_kg_m3']:.4e} kg/m³")
    passed += 1

    # Test 2: MUGE gravity — Newtonian dominance at small r
    muge = MUGEGravityField()
    res2 = muge.compute({'M': M_SUN, 'r': 6.96e8})
    assert res2['g_total_m_s2'] > 0
    assert abs(res2['g_newton'] - 274.0) < 5.0  # Solar surface gravity ~274 m/s²
    print(f"  [PASS] Test 2: Solar surface g = {res2['g_newton']:.2f} m/s²")
    passed += 1

    # Test 3: MUGE has 8 components
    assert len(res2['components']) == 8
    print(f"  [PASS] Test 3: MUGE 8-term components present")
    passed += 1

    # Test 4: Volumetric field (small grid)
    vol = MUGEVolumetricField()
    vres = vol.compute({
        'M': 1e11 * M_SUN, 'rho_s': 1e-21, 'r_s': 20 * KPC,
        'box_size': 50 * KPC, 'nx': 8, 'ny': 8, 'nz': 8,
    })
    assert vres['n_points'] == 512
    assert len(vres['density_cube']) == 512
    assert len(vres['g_magnitude_cube']) == 512
    assert len(vres['v_circ_cube']) == 512
    print(f"  [PASS] Test 4: Volumetric 8³ = {vres['n_points']} points, {vres['elapsed_ms']:.1f} ms")
    passed += 1

    # Test 5: Density cube has expected range
    rho_max = max(vres['density_cube'])
    rho_min = min(vres['density_cube'])
    assert rho_max > rho_min
    assert rho_max > 0
    print(f"  [PASS] Test 5: Density range [{rho_min:.3e}, {rho_max:.3e}] kg/m³")
    passed += 1

    # Test 6: Gravity vector components present
    assert len(vres['gx_cube']) == 512
    assert len(vres['gy_cube']) == 512
    assert len(vres['gz_cube']) == 512
    print(f"  [PASS] Test 6: 3D gravity vector cubes present")
    passed += 1

    # Test 7: Multi-system batch
    multi = MUGEMultiSystemVolume()
    mres = multi.compute({'systems': [
        {'name': 'NGC 6302', 'M': 0.64 * M_SUN, 'rho_s': 1e-18, 'r_s': 0.5 * KPC,
         'box_size': 2 * KPC, 'nx': 4, 'ny': 4, 'nz': 4},
        {'name': 'Crab Nebula', 'M': 1.4 * M_SUN, 'rho_s': 1e-19, 'r_s': 1 * KPC,
         'box_size': 5 * KPC, 'nx': 4, 'ny': 4, 'nz': 4},
    ]})
    assert mres['n_systems'] == 2
    for s in mres['system_summaries']:
        assert 'density_max' in s
        assert 'v_circ_max' in s
    print(f"  [PASS] Test 7: Multi-system batch: {mres['n_systems']} systems")
    passed += 1

    # Test 8: Rotation curve
    rc = MUGERotationCurve()
    rcres = rc.compute({
        'M': 1e11 * M_SUN, 'r_min': 1 * KPC, 'r_max': 100 * KPC, 'n_radii': 50,
    })
    assert rcres['n_radii'] == 50
    assert len(rcres['v_circ_km_s']) == 50
    assert max(rcres['v_circ_km_s']) > 0
    print(f"  [PASS] Test 8: Rotation curve, v_max = {max(rcres['v_circ_km_s']):.1f} km/s")
    passed += 1

    # Test 9: Flatness ratio meaningful
    assert 0 < rcres['flatness_ratio'] <= 1.0
    print(f"  [PASS] Test 9: Flatness ratio = {rcres['flatness_ratio']:.4f}")
    passed += 1

    # Test 10: Primary equations present everywhere
    assert len(vres['primary_equations']) >= 1
    assert len(mres['primary_equations']) >= 1
    assert len(rcres['primary_equations']) >= 1
    print(f"  [PASS] Test 10: Primary equations in all outputs")
    passed += 1

    print("-" * 72)
    print(f"Results: {passed}/{passed + failed} passed, {failed} failed")
    return passed, failed


if __name__ == '__main__':
    _run_tests()
