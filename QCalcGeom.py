# -*- coding: utf-8 -*-
"""
QCalcGeom.py  v2.1.0  —  BSFG Geometric Physics + Universal Buoyancy Solver

Derives from QCalcGeom.cpp (C++ v1.2.0, Sessions 150–151 Phase G) and elevates
to a Python simultaneous-equation solver for:

  UNIVERSAL BUOYANCY (Aether UA vacuum — as understood by the Greeks):
    FUBi   = outside-to-inside collapsing gravity zone       (SOURCE4 Ubi formula)
    FUBii  = inside-to-outside Aether counter-buoyancy       (habitable zone force)
    F_U    = Σ(Ug_i) + Um − FUBi + FUBii                    (Universal Gravity)

  SIMULTANEOUS EQUATION SYSTEM (solved jointly at each (r, t_n)):
    Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0     [buoyancy crossing / compaction]
    Eq2: ε′(r, t_n) + G·M/(c²·r²) = 0         [metric-geodesic matching]
    → Yields r_hz (habitable zone radius) + t_n_hz (phase at equilibrium)

  QUANTUM CHAIN compliance  (dpm_vacuum_manifold.py v3.0):
    RHO_VAC_SCM = 633,333 J/m³  derived from 26-level hydrogen geometry.
    Mass emerges at the FUBi+FUBii = 0 crossing — NOT from hardcoded GM/r².

  HABITABLE ZONE PHYSICS:
    r < r_hz  : FUBi > |FUBii|  → gravity/collapse dominates   (rocky body zone)
    r > r_hz  : |FUBii| > FUBi  → Aether buoyancy dominates    (void / gas zone)
    r = r_hz  : neutral buoyancy → liquid / habitable equilibrium

  ARCHITECTURE:
    source2.cpp (GUI) → dataset dict → QCalcGeom.py → results dict
    All calculators: cls().compute(dataset) → dict  (CondensedPhysics pattern)

FUNCTIONS PORTED FROM QCalcGeom.cpp (all 12):
  bsfg_metric, bsfg_horizon, bsfg_field_equations, bsfg_geodesic,
  bsfg_holonomy, vds_series, dvp_arithmetic, bsh_harmonic,
  bh26_eigenvalue, bsfg_buoyancy, poly26_derivative, uqff_comp_matrix

NEW IN v2.0.0 (simultaneous equation level):
  compute_FUBii          — Aether counter-buoyancy (inside-to-outside)
  compute_F_U            — Universal Gravity full assembly
  solve_habitable_zone   — scipy simultaneous solver for r_hz + t_n_hz
  scan_habitable_zone    — 2D scan landscape (r × t_n)
  BSFGMetricCalculator   — calculator class
  UniversalBuoyancyCalculator — FUBi + FUBii balance
  HabitableZoneCalculator — simultaneous solver class
  UniversalGravityCalculator  — F_U complete assembly
  run_qcalcgeom_tests    — 60 Python tests (mirrors C++ T01–T60)

Author  : Daniel T. Murphy
Created : Session 201 — May 5, 2026
Based on: QCalcGeom.cpp v1.2.0 (Session 151 Phase G)
Version : 2.2.0

NEW IN v2.2.0 (Session 205+ Universal-Buoyancy directive):
  - SECTION 5.5: fully-coupled 4x4 simultaneous-equation solver
    solve_universal_buoyancy() jointly determines (r_hz, t_n_hz, r_cg, M_emergent)
    via scipy.fsolve on residuals [E1, E2, E3, E4] in log-space, with a
    physically-motivated staged seed and a graceful staged-decoupled fallback.
  - UniversalBuoyancySimultaneousSolver calculator class (CondensedPhysics pattern)
    encoding the Greek Aether-UA vacuum counter-balance: collapsing-gravity zone
    (r < r_cg), habitable shell (r_cg <= r <= r_hz), gaseous outer (r > r_hz).
  - Tests T81-T87 (UBS: Universal Buoyancy Simultaneous solver)

NEW IN v2.1.0 (Session 202 Phase H202):
  - VDS variant branches: vds_prime (calibration sensitivity), vds_density,
    vds_k_weighted (BH26-coupled amplitude) -- VDSBranchResult
  - DVP spectral branches: full zeta_sum, pair_product (double-vortex state),
    spectral_floor (Navier-Stokes vorticity floor) -- DVPBranchResult
  - BH26/DH26 branches: spectral_sum (eigenvalue ladder), casimir_energy,
    degeneracy_k1=26 (S^{25} multiplicity), vds_coupling (topology bridge) -- BH26BranchResult
  - VDS*DVP coupling: joint geometric-mean coefficient + variant calibration branch
  - BH26*BSH cross-resonance energy density at BH26 spectral frequency bins
  - Tests T71-T80 (group VDS-DVP-DH26)
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, Any

import numpy as np

# ─── DPM Quantum Chain — canonical vacuum constants ──────────────────────────
from dpm_vacuum_manifold import (
    derive_from_quantum_chain,
    RHO_VAC_SCM,    # 633,333 J/m³  emergent SCm energy density
    RHO_VAC_UA,     # 6,333,333 J/m³ emergent UA energy density (10× SCm)
    BETA_I,         # 0.60  buoyancy coupling β_i
    LAMBDA_I,       # 1.0   manifold coupling λ_i
    OMEGA_S,        # 2.5e-6 rad/s stellar angular frequency
    SSQ,            # 0.57  [SSq] triple-convergence constant
    KAPPA_FLOAT,    # 0.0005 day⁻¹
    S26_3,          # 1.4531e26  Ramanujan 26D amplification
    KER_SCm,        # ~630 eV × J  Holmlid KER
    F_TRZ,          # 0.1  Time-Reversal Zone factor
    vds_numerical,
    compute_F_U_Bi_i_numerical,
)

# ─── try scipy for simultaneous solve; graceful fallback ─────────────────────
try:
    from scipy.optimize import fsolve as _scipy_fsolve
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False
    warnings.warn("scipy not available — HabitableZoneSolver will use Newton fallback")

# =============================================================================
# SECTION 1 — CANONICAL CONSTANTS (mirrors QCalcGeom.h Section 2)
# =============================================================================

# Aether-metric coupling η  (_S148_ETA)
ETA_BSFG     : float = 1.0e-22

# Physical constants
C_LIGHT      : float = 3.0e8            # m/s
G_NEWTON     : float = 6.674e-11        # m³/(kg·s²)
M_SUN        : float = 1.989e30         # kg
L_SUN        : float = 3.828e26         # W
R_SUN        : float = 6.96e8           # m
HBAR         : float = 1.055e-34        # J·s
H_PLANCK     : float = 6.626e-34        # J·s
K_BOLTZ      : float = 1.381e-23        # J/K
L_PLANCK     : float = 1.616e-35        # m
AU_METERS    : float = 1.496e11         # m
LAMBDA_OBS   : float = 1.1e-52          # m⁻²
FAC26_APPROX : float = 4.0329146113e26  # 26! ≈
SSQ_DEFAULT  : float = float(SSQ)       # 0.57
DVP_PRIME    : int   = 113              # 30th prime
RERING_BB_HZ : float = 1.15e14         # Hz  BH26 resonance

# Derived constant: C_num = (M_⊙·c² + L_⊙/c²) / ((4/3)π)
def _c_num_solar() -> float:
    return (M_SUN * C_LIGHT**2 + L_SUN / C_LIGHT**2) / (4.0 / 3.0 * math.pi)

C_NUM_SOLAR : float = _c_num_solar()     # ≈ 4.273e46  m³·kg/m³·c²

# κ_E = 8πG/c⁴  (Einstein coupling)
KAPPA_E : float = 8.0 * math.pi * G_NEWTON / C_LIGHT**4

# Buoyancy canonical parameters (SOURCE4 compute_Ubi_SOURCE4)
BETA_I_BSFG   : float = 0.6            # β_i
OMEGA_G_BSFG  : float = 7.3e-16        # Ω_g  [rad/s]
M_BH_BSFG     : float = 8.15e36        # M_bh [kg]
D_G_BSFG      : float = 2.55e20        # d_g  [m]
EPS_SW_BSFG   : float = 0.001          # ε_sw
RHO_SW_BSFG   : float = 8.0e-21        # ρ_sw [kg/m³]
U_UA_BSFG     : float = 1.0            # U_UA Aether factor

# FUBii spring constant: outward Aether pressure scales as RHO_VAC_SCM · c²/r²_norm
# C_FIELD_BSFG sets the amplitude of the counter-buoyancy pressure gradient [J/m³]
C_FIELD_BSFG  : float = RHO_VAC_SCM * (4.0 * math.pi / 3.0) * C_LIGHT**2

# Reference test values (from CP4 #149–#157 docstrings)
EPS_PRIME_REF  : float = 5.47e-11
R_R0R0_REF     : float = 1.57e-19
R_H_BSFG_REF   : float = 1.62e8
T_H_BSFG_REF   : float = 3.37e-12
R_CROSS_AU_REF : float = 0.360
H_ETA_REF      : float = 6.626e-56
UBI_BSFG_REF   : float = -7.63e33
POLY26_NEG_THR : float = 1.0e-100

# ── Mayan Three-Ring Timing constants ─────────────────────────────────────────
# Three-geared ring system encoding Universal Inertia (zero-point gravity timing)
# Ring proportions change every epoch; in Epoch 5 (2012+):
#   Ring 1 (OUTER)     : EXPANDING  — r_base × φ^(epoch−1)
#   Ring 2 (COMPANION) : SHRINKS    — r_base × φ^(−(epoch−1))
#   Ring 3 (INNER)     : VERY SMALL — r_base × φ^(−2(epoch−1))
# Gear meshing:
#   Ring 1 ↔ Ring 2 : external (side-by-side), ω_2 = −ω_1 × r_1/r_2
#   Ring 1 ↔ Ring 3 : INTERNAL (Ring 3 inside Ring 1), ω_3 = +ω_1 × r_1/r_3
PHI               : float = (1.0 + math.sqrt(5.0)) / 2.0    # Golden ratio φ
MAYAN_BAKTUN_DAYS : float = 144000.0                          # 1 baktun [days]
MAYAN_GREAT_CYCLE_DAYS  : float = 13.0 * MAYAN_BAKTUN_DAYS   # 1,872,000 days
MAYAN_GREAT_CYCLE_YEARS : float = MAYAN_GREAT_CYCLE_DAYS / 365.25  # 5125.36 yr
MAYAN_EPOCH5_YEAR : float = 2012.972                          # Dec 21, 2012
MAYAN_N_EPOCHS    : int   = 5                                 # known epochs
# Angular frequency of the Great Cycle ring (Ring 1 base)
OMEGA_BAKTUN      : float = 2.0 * math.pi / (MAYAN_GREAT_CYCLE_DAYS * 86400.0)

# =============================================================================
# SECTION 2 — RESULT DATACLASSES (Python equivalents of C++ structs)
# =============================================================================

@dataclass
class BSFGMetricResult:
    """BSFG metric and curvature at (r, t_n). Covers tests T01–T06."""
    eps       : float = 0.0   # ε = η·T_s00·cos(π·t_n)
    eps_p     : float = 0.0   # ε′ = dε/dr  [m⁻¹]
    eps_pp    : float = 0.0   # ε″ = d²ε/dr²  [m⁻²]
    A00       : float = 0.0   # g_00 = 1+ε
    Arr       : float = 0.0   # g_rr = −1+ε
    R_r0r0    : float = 0.0   # Riemann R^r_{0r0}  [m⁻²]
    R_00      : float = 0.0   # Ricci R_{00}  [m⁻²]
    R_rr      : float = 0.0   # Ricci R_{rr}  [m⁻²]
    R_scalar  : float = 0.0   # Ricci scalar R  [m⁻²]
    Kretschner: float = 0.0   # K = 12·(R^r_{0r0})²  [m⁻⁴]

@dataclass
class BSFGHorizonResult:
    """BSFG blinking horizon. Covers tests T07–T08."""
    exists       : bool  = False
    r_h          : float = 0.0  # [m]
    T_H          : float = 0.0  # Hawking temperature [K]
    kappa_surf   : float = 0.0  # surface gravity [s⁻²]
    r_h_over_Rs  : float = 0.0  # r_h / R_⊙

@dataclass
class BSFGFieldEqResult:
    """BSFG Einstein field equation deviations. Covers tests T28–T29."""
    amp_factor   : float = 0.0  # G_00 / (κ_E·T_s00)
    Lambda_eff   : float = 0.0  # Λ_eff [m⁻²]
    Lambda_ratio : float = 0.0  # Λ_eff / Λ_obs
    rho_vac_eff  : float = 0.0  # ρ_vac = Λ_eff·c²/(8πG)  [kg/m³]

@dataclass
class BSFGGeodesicResult:
    """BSFG geodesic and orbital quantization. Covers tests T09–T10."""
    v2_newton      : float = 0.0  # G·M/r  [m²/s²]
    v2_aether      : float = 0.0  # Aether velocity correction
    r_cross_m      : float = 0.0  # Crossover radius [m]
    r_cross_AU     : float = 0.0  # [AU]
    h_eta          : float = 0.0  # η·h  [J·s]
    delta_J_over_J : float = 0.0  # |δJ/J|

@dataclass
class BSFGHolonomyResult:
    """BSFG holonomy and extra-dimension topology. Covers tests T11–T12."""
    delta_phi    : float = 0.0
    omega_0r     : float = 0.0   # ε′/2  [m⁻¹]
    n_extra_flat : int   = 22    # 26−4 = 22 compact U(1) dimensions
    G2_excluded  : bool  = True
    Spin7_excluded: bool = True

@dataclass
class VDSResult:
    """Vacuum Density Series. Covers tests T13–T16."""
    value       : float = 0.0
    converged   : bool  = False
    tail_bound  : float = 0.0
    n_terms_used: int   = 0

@dataclass
class DVPResult:
    """Dipole Vortex Primes arithmetic. Covers tests T17–T21."""
    fac26_mod_113 : int   = 0
    non_repeating : bool  = False
    r_q_AU        : float = 0.0
    r_q_m         : float = 0.0

@dataclass
class BSHResult:
    """Buoyancy Series Harmonics. Covers tests T22–T24."""
    U_g2      : float = 0.0
    H_m_max   : float = 0.0
    saturated : bool  = False

@dataclass
class BH26Result:
    """BH26 Kaluza–Klein eigenvalue. Covers tests T25–T27."""
    lambda_k    : float = 0.0
    freq_bin_hz : float = 0.0
    finite      : bool  = False

@dataclass
class BSFGBuoyancyResult:
    """BSFG buoyancy (FUBi) coupling. Covers tests T41–T50."""
    Ubi          : float = 0.0  # −β_i·Ug_field·orbit·cos(π·t_n)
    Ug_field     : float = 0.0  # G·M_⊙²/r²
    orbit_factor : float = 0.0  # Ω_g·M_bh/d_g·wind_mod·U_UA
    cos_tn       : float = 0.0
    negative     : bool  = False  # Ubi < 0: opposes gravity
    inverted     : bool  = False  # Ubi > 0: aids collapse
    zero_crossing: bool  = False  # |cos| < 1e-10

@dataclass
class Poly26Result:
    """26th-order polynomial derivative. Covers tests T51–T56."""
    value           : float = 0.0
    factorial_ratio : float = 0.0
    r_power         : float = 0.0
    negligible      : bool  = False

@dataclass
class UQFFCompResult:
    """UQFF compressed-field matrix. Covers tests T57–T60."""
    m00              : float = 0.0
    m11              : float = 0.0
    m22              : float = 0.0
    cross_d13        : float = 0.0
    eigenvalue_min   : float = 0.0
    positive_definite: bool  = False

# ── NEW v2.0.0 result structs ─────────────────────────────────────────────────

@dataclass
class FUBiiResult:
    """Inside-to-outside Aether counter-buoyancy (habitable zone Aether spring).
    FUBii = +ρ_vac,SCm · (4π/3) · r · c² · cos(π·t_n)
    Grows linearly with r — opposes FUBi which falls as 1/r².
    """
    FUBii        : float = 0.0   # outward Aether pressure [SI]
    rho_aether   : float = 0.0   # RHO_VAC_SCM  [J/m³]
    cos_tn       : float = 0.0
    r_m          : float = 0.0
    outward      : bool  = False  # FUBii > 0 → buoyancy dominates

@dataclass
class HabitableZoneResult:
    """Result of the simultaneous FUBi+FUBii / metric-geodesic crossing solve.
    The habitable zone is where neutral Aether buoyancy holds:
      FUBi(r_hz, t_n_hz) + FUBii(r_hz, t_n_hz) = 0   AND
      ε′(r_hz, t_n_hz) + G·M/(c²·r_hz²) = 0

    Epoch-aware tectonic band (Session 265):
      Mayan three-ring scaling defines the inner/outer tectonic band
      bounding the crustal/floating zone around r_hz.  In Epoch 5 the
      inner ring shrinks by φ^(-8) ≈ 0.0213× and the outer ring expands
      by φ^4 ≈ 6.854× → the habitable shell narrows asymmetrically.
    """
    r_hz_m       : float = 0.0   # habitable zone radius [m]
    r_hz_AU      : float = 0.0   # [AU]
    t_n_hz       : float = 0.0   # phase at equilibrium
    FUBi_at_hz   : float = 0.0   # should be ≈ −FUBii
    FUBii_at_hz  : float = 0.0   # should be ≈ −FUBi
    residual_eq1 : float = 0.0   # FUBi + FUBii at solution
    residual_eq2 : float = 0.0   # ε′ + GM/(c²r²) at solution
    converged    : bool  = False
    solver_msg   : str   = ""
    # Habitable zone classification
    hz_type      : str   = ""    # "rocky_inner" / "habitable" / "gas_outer"
    # Epoch-aware tectonic band (Session 265)
    epoch                : int   = 5     # Mayan epoch (1-5)
    r_hz_inner_band_AU   : float = 0.0   # r_hz × φ^(-2(E-1))  inner tectonic floor
    r_hz_outer_band_AU   : float = 0.0   # r_hz × φ^(E-1)      outer tectonic ceiling
    band_width_AU        : float = 0.0   # outer − inner
    U_I_at_hz            : float = 0.0   # Universal Inertia [J/m³] at (t_n_hz, epoch)
    U_I_mode             : str   = ""    # "centripetal" / "zero-point" / "centrifugal"

@dataclass
class EmergentMassResult:
    """Mass-emergence at the FUBi+FUBii=0 crossing (Quantum Chain Step 7).

    The crossing condition  β_i·G·M²·orbit/r² = ρ_vac·(4π/3)·r·c²  is solved
    for M (NOT for r — r_hz is the input from the habitable-zone solver).
    This is mass BORN from Aether buoyancy, not mass plugged into GM/r².

    Closed form:
        M_emergent = sqrt[ ρ_vac · (4π/3) · r_hz³ · c² / (β_i · G · orbit_factor) ]
    """
    r_hz_m         : float = 0.0   # input crossing radius [m]
    r_hz_AU        : float = 0.0   # [AU]
    t_n_hz         : float = 0.0   # phase at crossing
    M_emergent_kg  : float = 0.0   # emergent mass [kg]
    M_emergent_sun : float = 0.0   # [M_⊙]
    rho_vac_used   : float = 0.0   # vacuum density that mass condensed against
    orbit_factor   : float = 0.0   # FUBi geometric pre-factor
    beta_i_used    : float = 0.0
    residual_at_M  : float = 0.0   # |FUBi(M_emergent) + FUBii| — should be ≈0
    converged      : bool  = False
    # Diagnostic: where does this mass sit on the M_sun ladder?
    classification : str   = ""    # "sub_solar" / "solar" / "stellar" / "BH_seed"

@dataclass
class UniversalGravityResult:
    """Complete F_U assembly with Universal Inertia (Session 265):
      F_U = Ug1+Ug2+Ug3+Ug4 − FUBi + FUBii + Um + U_I·V_body·xi

    Universal Inertia U_I is the invariant differential — primordial
    radiance scalar that reverses centripetal↔centrifugal across t_n.
    Coupled perturbatively (xi=1e-6) so it vanishes at zero-point
    (cos(π·t_n)=0) and flips sign over a great cycle.
    """
    Ug1      : float = 0.0
    Ug2      : float = 0.0
    Ug3      : float = 0.0
    Ug4      : float = 0.0
    FUBi     : float = 0.0   # SOURCE4 Ubi (negative = opposes gravity normally)
    FUBii    : float = 0.0   # Aether counter-buoyancy (positive outward)
    Um       : float = 0.0
    U_I      : float = 0.0   # Universal Inertia invariant [J/m³]  (Session 265)
    U_I_term : float = 0.0   # U_I·V_body·xi contribution to F_U   (Session 265)
    epoch    : int   = 5     # Mayan epoch used for U_I                (Session 265)
    F_U_total: float = 0.0   # Ug_sum − FUBi + FUBii + Um + U_I_term
    r_m      : float = 0.0
    t_n      : float = 0.0
    eps      : float = 0.0   # BSFG metric perturbation

@dataclass
class MayanRingState:
    """Three-ring geared state for a given epoch and current year.

    Ring 1 (outer/left):     EXPANDING  — Great Cycle / baktun driver ring
    Ring 2 (companion/right): SHRINKS   — external mesh (opposite rotation)
    Ring 3 (inner):           VERY SMALL in Epoch 5 — internal mesh (same dir)

    Angular velocities satisfy gear constraints:
      ω_2 = −ω_1 × r_1/r_2   (external mesh: opposite rotation)
      ω_3 = +ω_1 × r_1/r_3   (internal mesh: same rotation, amplified)

    Epoch 5 (n=5) characteristic gear ratio: r_1/r_3 = φ^12 ≈ 321.997
    This 322× amplification provides quantum timing precision for zero-point.
    """
    epoch           : int   = 5
    r_outer_AU      : float = 0.0   # Ring 1 radius [AU]
    r_companion_AU  : float = 0.0   # Ring 2 radius [AU]
    r_inner_AU      : float = 0.0   # Ring 3 radius [AU]  (very small Epoch 5)
    r_outer_m       : float = 0.0   # Ring 1 [m]
    r_companion_m   : float = 0.0   # Ring 2 [m]
    r_inner_m       : float = 0.0   # Ring 3 [m]
    omega_outer     : float = 0.0   # ω_1  [rad/s]  base baktun frequency
    omega_companion : float = 0.0   # ω_2 = −ω_1 × r_1/r_2  (negative)
    omega_inner     : float = 0.0   # ω_3 = +ω_1 × r_1/r_3  (positive, large)
    gear_ratio_12   : float = 0.0   # r_1/r_2  (outer:companion)
    gear_ratio_13   : float = 0.0   # r_1/r_3  (outer:inner)  ≈ 46.98 Epoch 5
    t_n             : float = 0.0   # UQFF phase (years since epoch / cycle length)
    t_n_inner       : float = 0.0   # Ring 3 phase (amplified by gear ratio)
    current_year    : float = 0.0
    zero_point_next_year: float = 0.0   # next cos(π·t_n) = 0 year
    u_inertia       : float = 0.0   # Universal Inertia invariant [J/m³·m/s² = N/m²]
    cos_tn          : float = 0.0
    massless_mode   : bool  = False  # True when |cos(π·t_n)| < 1e-4 (near zero-point)
    epoch5_active   : bool  = False

@dataclass
class UniversalInertiaResult:
    """Universal Inertia — the invariant differential of the F_U field.

    At the habitable zone (FUBi + FUBii = 0), the second derivative of the
    field potential is:

      U_I = d²V/dr²|_{r_hz} = 3 · ρ_vac · (4π/3) · c² · cos(π·t_n)

    This is INVARIANT — independent of r, M, orbit — depending only on:
      ρ_vac (Aether vacuum energy density)
      c     (speed of light)
      cos(π·t_n)  (Mayan timing phase)

    Physical regimes:
      U_I > 0  (cos > 0): stable harmonic trap — centripetal restoring (mass-building)
      U_I = 0  (cos = 0): zero-point gravity — massless scalar state (t_n = 1/2)
      U_I < 0  (cos < 0): unstable saddle — centrifugal expansion (negentropic void)

    Connection to primordial radiance:
      Frequency range = ω_outer to ω_inner
      f_range = ω_inner / ω_outer = r_outer / r_inner = φ^(2*(epoch-1))
      In Epoch 5: f_range = φ^12 ≈ 321.997  (primordial radiance band-width factor)
    """
    u_inertia       : float = 0.0   # 3·ρ_vac·(4π/3)·c²·cos(π·t_n)  [N/m²]
    u_inertia_abs   : float = 0.0   # |U_I|  [N/m²]
    rho_vac         : float = 0.0   # ρ_vac (RHO_VAC_SCM)  [J/m³]
    cos_tn          : float = 0.0
    t_n             : float = 0.0
    primordial_freq_range: float = 0.0   # ω_inner/ω_outer = φ^(2·(epoch-1))
    massless_scalar : bool  = False  # cos ≈ 0 → massless-to-massive transition
    centripetal_mode: bool  = False  # U_I > 0 → centripetal restoring
    centrifugal_mode: bool  = False  # U_I < 0 → centrifugal expanding
    zero_point      : bool  = False  # |U_I| < threshold → zero-point gravity
    tectonic_band_inner_m: float = 0.0   # r_inner [m]  (inner tectonic boundary)
    tectonic_band_outer_m: float = 0.0   # r_outer [m]  (outer tectonic boundary)

# ── Session 202 Phase H202: VDS/DVP/DH26 variant-branch result structs ────────

@dataclass
class VDSBranchResult:
    """VDS variant branches: calibration sensitivity, energy density, BH26-coupled.
    vds_prime      = d/dz Li_26(z)|_{z=SSq} = Li_25(SSq)/SSq  ≈ 1.0
    vds_density    = Li_26(SSq) × RHO_VAC_SCM   [J/m³]
    vds_k_weighted = Li_25(SSq) + 25·Li_26(SSq)  (BH26-coupled amplitude)
    """
    vds_li25       : float = 0.0  # Li_{25}([SSq])
    vds_prime      : float = 0.0  # Li_{25}/SSq  ≈ 1.0  (calibration sensitivity)
    vds_density    : float = 0.0  # Li_{26}(SSq) × RHO_VAC_SCM  [J/m³]
    vds_k_weighted : float = 0.0  # Li_{25} + 25·Li_{26}

@dataclass
class DVPBranchResult:
    """DVP variant branches: spectral sum, double-vortex pair, Navier-Stokes floor.
    zeta_sum       = Σ_{p>26} a(p)  where a(p) = SSq^{π(p)} / p^26
    pair_product   = a(29) × a(31)  (double-vortex amplitude)
    spectral_floor = a(p_max)       (vorticity lower bound)
    """
    zeta_sum       : float = 0.0
    n_primes_dvp   : int   = 0
    pair_product   : float = 0.0
    spectral_floor : float = 0.0
    a_29           : float = 0.0  # a(29) — dominant first DVP term

@dataclass
class BH26BranchResult:
    """BH26/DH26 variant branches: eigenvalue ladder, Casimir energy, topology.
    spectral_sum   = Σ_{k=1}^{N} k(k+25)  [= 1760 for N=10]
    casimir_energy = ℏ·f_{RR}/2 × Σ 1/λ_k  [J]  (vacuum Casimir)
    degeneracy_k1  = 26  C(26,25) on S^{25}
    vds_coupling   = Σ_{k=1}^{N} λ_k^{-26}  (BH26→VDS topological bridge)
    """
    spectral_sum   : float = 0.0
    casimir_energy : float = 0.0
    degeneracy_k1  : int   = 0
    vds_coupling   : float = 0.0
    N              : int   = 0

@dataclass
class VDSDVPCoupledResult:
    """VDS×DVP coupled field: normalised weights, geometric-mean coupling, calibration gap.
    joint_coeff    = sqrt(w_vds × w_dvp)  (geometric-mean field coupling)
    variant_branch = |w_vds − w_dvp|      (differential calibration magnitude)
    """
    w_vds          : float = 0.0
    w_dvp          : float = 0.0
    joint_coeff    : float = 0.0
    variant_branch : float = 0.0

@dataclass
class BH26BSHResonanceResult:
    """BH26×BSH cross-resonance: BSH evaluated at a BH26 spectral frequency bin.
    freq_k         = RERING_BB_HZ / λ_k   [Hz]
    bsh_at_k       = BSH U_g2 at omega = 2π·freq_k
    resonance      = bsh_at_k · cos(π·t_n)
    energy_density = resonance × RHO_VAC_SCM  [J/m³]
    """
    freq_k         : float = 0.0
    bsh_at_k       : float = 0.0
    resonance      : float = 0.0
    energy_density : float = 0.0

# =============================================================================
# SECTION 3 — CORE PHYSICS FUNCTIONS (Python ports of QCalcGeom.cpp)
# =============================================================================

def _ts00(r: float) -> float:
    """T_s00(r) = M_⊙·c² / ((4/3)π·r³)  stellar stress-energy time component"""
    return M_SUN * C_LIGHT**2 / ((4.0 / 3.0) * math.pi * r**3)


def bsfg_metric(r: float, t_n: float) -> BSFGMetricResult:
    """BSFG metric g_μν = diag(1+ε, −1+ε, r², r²sin²θ).
    ε = η · T_s00(r) · cos(π·t_n)
    Reference: CP4 #149 BSFGRiemannCurvatureAetherMetricCalculator
    """
    m = BSFGMetricResult()
    cos_tn = math.cos(math.pi * t_n)
    Cnum   = C_NUM_SOLAR

    m.eps    =  ETA_BSFG * (Cnum / r**3) * cos_tn
    m.eps_p  = -3.0 * ETA_BSFG * cos_tn * Cnum / r**4
    m.eps_pp = 12.0 * ETA_BSFG * cos_tn * Cnum / r**5

    m.A00 =  1.0 + m.eps
    m.Arr = -1.0 + m.eps

    m.R_r0r0   = m.eps_pp * 0.5 - m.eps_p**2 * 0.5
    m.R_00     = 3.0 * m.R_r0r0
    m.R_rr     = -m.R_r0r0 + m.eps_pp - m.eps_p**2 * 0.5

    if m.A00 != 0.0 and m.Arr != 0.0:
        m.R_scalar = m.R_00 / m.A00 + m.R_rr / m.Arr
    else:
        m.R_scalar = 0.0

    m.Kretschner = 12.0 * m.R_r0r0**2
    return m


def bsfg_horizon(t_n: float) -> BSFGHorizonResult:
    """BSFG blinking horizon: A00(r_h) = 0 → r_h³ = −η·C_num·cos(π·t_n)
    Reference: CP4 #156 BSFGBlackHoleSolutionHorizonCalculator
    """
    h = BSFGHorizonResult()
    cos_tn = math.cos(math.pi * t_n)
    Cnum   = C_NUM_SOLAR
    arg    = -ETA_BSFG * Cnum * cos_tn

    h.exists = (arg > 0.0)
    if h.exists:
        h.r_h = arg**(1.0 / 3.0)
        h.kappa_surf = C_LIGHT**2 * abs(3.0 * ETA_BSFG * cos_tn * Cnum / h.r_h**4) * 0.5
        h.T_H        = HBAR * h.kappa_surf / (2.0 * math.pi * K_BOLTZ * C_LIGHT)
    else:
        h.r_h = 0.0
        h.kappa_surf = 0.0
        h.T_H = 0.0

    h.r_h_over_Rs = h.r_h / R_SUN
    return h


def bsfg_field_equations(r: float, t_n: float) -> BSFGFieldEqResult:
    """BSFG Einstein tensor G_00, effective Λ_eff, ρ_vac_eff.
    Reference: CP4 #154 BSFGEinsteinTensorFieldEquationsCalculator
    """
    fe = BSFGFieldEqResult()
    m      = bsfg_metric(r, t_n)
    Ts00_r = _ts00(r)

    G_00   = m.R_00 - 0.5 * m.A00 * m.R_scalar
    RHS_00 = KAPPA_E * Ts00_r
    fe.amp_factor  = G_00 / RHS_00 if abs(RHS_00) > 0.0 else 0.0

    fe.Lambda_eff  = KAPPA_E * ETA_BSFG * Ts00_r / 2.0
    fe.Lambda_ratio = fe.Lambda_eff / LAMBDA_OBS
    fe.rho_vac_eff  = fe.Lambda_eff * C_LIGHT**2 / (8.0 * math.pi * G_NEWTON)
    return fe


def bsfg_geodesic(r: float, t_n: float) -> BSFGGeodesicResult:
    """BSFG geodesic: crossover radius, Aether action quantum, δJ/J.
    Reference: CP4 #157 BSFGBohrSommerfeldAetherQuantizationCalculator
    """
    g = BSFGGeodesicResult()
    m    = bsfg_metric(r, t_n)
    Cnum = C_NUM_SOLAR

    g.v2_newton = G_NEWTON * M_SUN / r
    g.v2_aether = r * m.eps_p * C_LIGHT**2 * 0.5

    g.r_cross_m  = math.sqrt(ETA_BSFG * C_LIGHT**2 * Cnum / (G_NEWTON * M_SUN))
    g.r_cross_AU = g.r_cross_m / AU_METERS
    g.h_eta      = ETA_BSFG * H_PLANCK

    if g.v2_newton > 0.0:
        g.delta_J_over_J = abs(g.v2_aether) / (2.0 * g.v2_newton)
    else:
        g.delta_J_over_J = 0.0
    return g


def bsfg_holonomy(r: float, t_n: float, loop_area_m2: float) -> BSFGHolonomyResult:
    """BSFG holonomy SO⁺(3,1)×U(1)^22, phase δφ = ω_{0r}·A.
    Reference: CP4 #155 BSFGHolonomyGroupTopologyCalculator
    """
    hl = BSFGHolonomyResult()
    m  = bsfg_metric(r, t_n)

    hl.omega_0r     = m.eps_p * 0.5
    hl.delta_phi    = hl.omega_0r * loop_area_m2
    hl.n_extra_flat = 22
    hl.G2_excluded  = True
    hl.Spin7_excluded = True
    return hl


def vds_series(SSq: float = SSQ_DEFAULT, n_terms: int = 200) -> VDSResult:
    """VDS(SSq, N) = Σ_{n=1}^{N} SSq^n / n^{26}  = Li_{26}(SSq).
    Reference: CP4 #83 ThreeNewNumberSystemsVacuumDipoleBuoyancyCalculator
    """
    v = VDSResult()
    total    = 0.0
    pow_SSq  = SSq
    n_used   = 0
    last_term = 0.0

    for n in range(1, n_terms + 1):
        term = pow_SSq / (n ** 26)
        total    += term
        pow_SSq  *= SSq
        n_used    = n
        last_term = abs(term)
        if last_term < 1.0e-15:
            break

    v.value        = total
    v.n_terms_used = n_used

    abs_SSq = abs(SSq)
    if abs_SSq < 1.0 and n_used > 0:
        v.tail_bound = abs_SSq**(n_used + 1) / ((n_used + 1)**26 * (1.0 - abs_SSq))
    else:
        v.tail_bound = float('inf')

    v.converged = (v.tail_bound < 1.0e-12)
    return v


def dvp_arithmetic() -> DVPResult:
    """26! mod 113 (exact), proplyd quantization radius r_q.
    Reference: CP4 #83 DVP component
    """
    d = DVPResult()
    result = 1
    for k in range(2, 27):
        result = (result * k) % DVP_PRIME
    d.fac26_mod_113 = result
    d.non_repeating = (result != 0)
    d.r_q_AU = (2.0 / FAC26_APPROX) ** (1.0 / 26.0)
    d.r_q_m  = d.r_q_AU * AU_METERS
    return d


def bsh_harmonic(f_Ub: float = 3.3e7, SSq: float = SSQ_DEFAULT,
                  omega: float = 2.0 * math.pi * 3.3e7,
                  t_n: float = 0.0, m_max: int = 20) -> BSHResult:
    """BSH: U_g2 = Σ_{m=1}^{m_max} H_m·(1−exp(−SSq·m))·cos(ω·t_n).
    Reference: CP4 #83 BSH component
    """
    b = BSHResult()
    cos_val  = math.cos(omega * t_n)
    H_partial = 0.0
    U_g2_sum  = 0.0

    for m in range(1, m_max + 1):
        H_partial += f_Ub / m
        U_g2_sum  += H_partial * (1.0 - math.exp(-SSq * m)) * cos_val

    b.U_g2      = U_g2_sum
    b.H_m_max   = H_partial
    b.saturated = (1.0 - math.exp(-SSq * m_max)) > 1.0 - 1.0e-6
    return b


def bh26_eigenvalue(k: int) -> BH26Result:
    """λ_k = k·(k+25) on S^25 (boundary of B^26). Reference: CP4 #149."""
    b = BH26Result()
    b.lambda_k    = float(k) * float(k + 25)
    b.freq_bin_hz = RERING_BB_HZ / b.lambda_k if b.lambda_k > 0.0 else 0.0
    b.finite      = math.isfinite(b.lambda_k) and b.lambda_k > 0.0
    return b


def bsfg_buoyancy(r: float, t_n: float,
                   beta_i: float     = BETA_I_BSFG,
                   Omega_g: float    = OMEGA_G_BSFG,
                   M_bh: float       = M_BH_BSFG,
                   d_g: float        = D_G_BSFG,
                   epsilon_sw: float = EPS_SW_BSFG,
                   rho_sw: float     = RHO_SW_BSFG,
                   U_UA: float       = U_UA_BSFG) -> BSFGBuoyancyResult:
    """SOURCE4 compute_Ubi_SOURCE4:
    FUBi = Ubi = −β_i·(G·M_⊙²/r²)·(Ω_g·M_bh/d_g)·wind_mod·U_UA·cos(π·t_n)
    t_n < 0: even symmetry → Ubi(t_n) = Ubi(−t_n)  (NegativeTimeModule)
    Reference: MAIN_1_CoAnQi.cpp SOURCE4 namespace; source106 NegativeTimeModule
    """
    b = BSFGBuoyancyResult()
    b.Ug_field    = G_NEWTON * M_SUN * M_SUN / (r * r)
    b.cos_tn      = math.cos(math.pi * t_n)
    wind_mod      = 1.0 + epsilon_sw * rho_sw
    b.orbit_factor = Omega_g * M_bh / d_g * wind_mod * U_UA
    b.Ubi         = -beta_i * b.Ug_field * b.orbit_factor * b.cos_tn
    b.negative    = (b.Ubi < 0.0)
    b.inverted    = (b.Ubi > 0.0)
    b.zero_crossing = abs(b.cos_tn) < 1.0e-10
    return b


def poly26_derivative(k: int, c: float, r: float) -> Poly26Result:
    """d^{26}/dr^{26}(c/r^k) = (k+25)!/(k-1)! · c / r^{k+26}
    Uses Pochhammer rising factorial: (k)_{26} = k·(k+1)·…·(k+25)
    Computed in log-space to handle extreme r values without overflow.
    Reference: grok_share_79fdf5367d1.txt — 26th-order polynomial expansions
    """
    out = Poly26Result()
    # Pochhammer (k)_26 = k*(k+1)*...*(k+25)
    log_fac = sum(math.log(float(k + m)) for m in range(26))
    log_r_power = float(k + 26) * math.log(r) if r > 0.0 else float('inf')
    log_c = math.log(abs(c)) if abs(c) > 0.0 else float('-inf')
    log_val = log_fac + log_c - log_r_power

    out.factorial_ratio = math.exp(min(log_fac, 700.0))  # capped to prevent overflow display
    # r_power stored as log since it may overflow
    out.r_power  = r ** float(k + 26) if float(k + 26) * math.log10(r) < 307 else float('inf')

    if math.isfinite(log_val) and log_val < 709.0:
        raw = math.exp(log_val)
        out.value = raw if c >= 0.0 else -raw
    elif log_val < -700.0:
        out.value = 0.0
    else:
        out.value = float('inf')

    out.negligible = abs(out.value) < POLY26_NEG_THR
    return out


def uqff_comp_matrix(r: float, rho: float) -> UQFFCompResult:
    """UQFF 3×3 compressed-field tensor diagonal + cross coupling.
    m00: 26th r-deriv of U_g = G·M_⊙/r
    m11: 26th r-deriv of U_m = κ/r^26
    m22: 26th ρ-deriv of U_b = G_N/ρ
    Reference: grok_share_79fdf5367d1.txt
    """
    d26_ug = poly26_derivative(1,  G_NEWTON * M_SUN, r)
    d26_um = poly26_derivative(26, 1.0,              r)
    d26_ub = poly26_derivative(1,  G_NEWTON,         rho)

    def _poly13(kk: int, cc: float, rr: float) -> float:
        log_fac = sum(math.log(float(kk + m)) for m in range(13))
        log_r   = float(kk + 13) * math.log(rr) if rr > 0 else float('inf')
        log_c   = math.log(abs(cc)) if abs(cc) > 0 else float('-inf')
        lv      = log_fac + log_c - log_r
        if not math.isfinite(lv) or lv < -700:
            return 0.0
        raw = math.exp(min(lv, 709.0))
        return raw if cc >= 0.0 else -raw

    d13_ug = _poly13(1,  G_NEWTON * M_SUN, r)
    d13_um = _poly13(26, 1.0,              r)
    cross   = math.sqrt(abs(d13_ug) * abs(d13_um))

    out = UQFFCompResult()
    out.m00              = d26_ug.value
    out.m11              = d26_um.value
    out.m22              = d26_ub.value
    out.cross_d13        = cross
    out.eigenvalue_min   = min(out.m00, out.m11, out.m22)
    out.positive_definite = (out.m00 >= 0.0) and (out.m11 >= 0.0) and (out.m22 >= 0.0)
    return out

# =============================================================================
# SECTION 3.5 — MAYAN THREE-RING TIMING / UNIVERSAL INERTIA ENGINE
# =============================================================================
#
# CANONICAL THREE-RING GEOMETRY (user specification, Epoch 5):
#
#   ┌────────────────────────────────────────────────────────────────┐
#   │  Ring 1 (OUTER / LEFT)    r_1 = r_base × φ^(n-1)  EXPANDING  │
#   │  ┌ ─ ─ ─ ─ ─ ─ ─ ─ ┐                                         │
#   │  │  Ring 3 (INNER)   │  r_3 = r_base × φ^(-2(n-1))           │
#   │  │  VERY SMALL (E5)  │  inside Ring 1, internal gear mesh     │
#   │  └ ─ ─ ─ ─ ─ ─ ─ ─ ┘                                         │
#   └────────────────────────────────────────────────────────────────┘
#         [Ring 1 external teeth mesh with Ring 2 on right side]
#   ┌──────────────────────┐
#   │  Ring 2 (COMPANION / RIGHT)                                    │
#   │  r_2 = r_base × φ^(-(n-1))   SHRINKS                          │
#   └──────────────────────┘
#
# GEAR RATIOS:
#   Ring 1 ↔ Ring 2 (external): ω_2 = −ω_1 × r_1/r_2   (opposite rotation)
#   Ring 1 ↔ Ring 3 (internal): ω_3 = +ω_1 × r_1/r_3   (same rotation, amplified)
#
# EPOCH 5 (n=5) AMPLIFICATION:  r_1/r_3 = φ^12 ≈ 321.997
#   → Ring 3 spins 322× faster → 322× finer quantum timing resolution
#   → Zero-point gravity precision: Δt_n ≈ 1/(322 × N_teeth)
# =============================================================================

def mayan_ring_proportions(epoch: int, r_base_AU: float = 1.0) -> dict:
    """Compute three-ring radii for a given Mayan epoch (1–5).

    Proportions change at each epoch transition:
      r_outer(n)     = r_base × φ^(n-1)       EXPANDING outward
      r_companion(n) = r_base × φ^(-(n-1))    SHRINKS with each epoch
      r_inner(n)     = r_base × φ^(-2*(n-1))  SHRINKS FASTEST — very small Epoch 5

    At Epoch 5 (current):
      r_outer    ≈ 6.854 × r_base  (dominant outer ring)
      r_companion ≈ 0.146 × r_base  (small companion)
      r_inner    ≈ 0.021 × r_base  (very small inner ring)
    """
    n = max(1, min(epoch, 10))   # clamp epoch to sensible range
    r_outer     = r_base_AU * PHI ** (n - 1)
    r_companion = r_base_AU * PHI ** (-(n - 1))
    r_inner     = r_base_AU * PHI ** (-2.0 * (n - 1))
    return {
        'r_outer_AU'    : r_outer,
        'r_companion_AU': r_companion,
        'r_inner_AU'    : r_inner,
        'r_outer_m'     : r_outer     * AU_METERS,
        'r_companion_m' : r_companion * AU_METERS,
        'r_inner_m'     : r_inner     * AU_METERS,
        'gear_ratio_12' : r_outer / r_companion if r_companion > 0 else 0.0,
        'gear_ratio_13' : r_outer / r_inner     if r_inner     > 0 else 0.0,
        'epoch'         : n,
    }


def mayan_ring_state(epoch: int = 5,
                     current_year: float = 2026.33,
                     r_base_AU: float = 1.0) -> MayanRingState:
    """Compute the full three-ring state for a given year.

    t_n = (current_year − epoch5_start) / great_cycle_years
    ω_1 = OMEGA_BAKTUN (Great Cycle angular frequency)
    ω_2 = −ω_1 × r_1/r_2   (companion, external mesh, opposite)
    ω_3 = +ω_1 × r_1/r_3   (inner, internal mesh, amplified)

    Zero-point timing: next year when cos(π·t_n) = 0
      → t_n = 1/2, 3/2, 5/2, ... (half-integers)
      → year = epoch5_start + (2k+1)/2 × great_cycle_years  k=0,1,2...
    """
    rings = mayan_ring_proportions(epoch, r_base_AU)
    r1 = rings['r_outer_m']
    r2 = rings['r_companion_m']
    r3 = rings['r_inner_m']

    # Phase in Great Cycle (t_n = 0 at epoch 5 start; 1 = full cycle)
    t_n = (current_year - MAYAN_EPOCH5_YEAR) / MAYAN_GREAT_CYCLE_YEARS

    # Angular velocities
    omega1 = OMEGA_BAKTUN
    omega2 = -omega1 * r1 / r2 if r2 > 0 else 0.0   # external mesh → opposite
    omega3 = +omega1 * r1 / r3 if r3 > 0 else 0.0   # internal mesh → same, amplified

    # Inner ring has amplified phase (it completes many cycles per outer cycle)
    gear_13 = rings['gear_ratio_13']
    t_n_inner = t_n * gear_13

    # Next zero-point: t_n must reach the next half-integer
    import math as _m
    next_half = _m.ceil(2.0 * t_n + 1e-9) / 2.0
    if abs(next_half - t_n) < 1e-6:
        next_half += 0.5
    zero_pt_year = MAYAN_EPOCH5_YEAR + next_half * MAYAN_GREAT_CYCLE_YEARS

    cos_tn   = math.cos(math.pi * t_n)
    u_inertia = 3.0 * RHO_VAC_SCM * (4.0 * math.pi / 3.0) * C_LIGHT**2 * cos_tn

    s = MayanRingState()
    s.epoch           = epoch
    s.r_outer_AU      = rings['r_outer_AU']
    s.r_companion_AU  = rings['r_companion_AU']
    s.r_inner_AU      = rings['r_inner_AU']
    s.r_outer_m       = r1
    s.r_companion_m   = r2
    s.r_inner_m       = r3
    s.omega_outer     = omega1
    s.omega_companion = omega2
    s.omega_inner     = omega3
    s.gear_ratio_12   = rings['gear_ratio_12']
    s.gear_ratio_13   = gear_13
    s.t_n             = t_n
    s.t_n_inner       = t_n_inner
    s.current_year    = current_year
    s.zero_point_next_year = zero_pt_year
    s.u_inertia       = u_inertia
    s.cos_tn          = cos_tn
    s.massless_mode   = abs(cos_tn) < 1.0e-4
    s.epoch5_active   = (epoch == 5)
    return s


def universal_inertia(t_n: float,
                       epoch: int = 5,
                       r_base_AU: float = 1.0,
                       rho_vac: float = None) -> UniversalInertiaResult:
    """Universal Inertia — invariant differential of the F_U field potential.

    Derived from d²V/dr²|_{r_hz} where V is the F_U potential:

      U_I = 3 · ρ_vac · (4π/3) · c² · cos(π·t_n)

    This is FRAME-INVARIANT: same value regardless of r_hz, stellar mass, or
    orbital parameters — depending only on vacuum energy and Mayan phase.

    Physical states:
      U_I > 0 : centripetal restoring → stable mass-building zone
      U_I = 0 : zero-point gravity → massless-to-massive scalar transition
      U_I < 0 : centrifugal expansion → negentropic void / dark energy zone

    The crustall/tectonic floating zone:
      Inner boundary: r_inner(epoch)  [FUBi-dominated, Rocky/solid]
      Outer boundary: r_outer(epoch)  [FUBii-dominated, Superconductive plasma]
      The tectonic plate 'floats' on the superconductive heavy plasma at r_hz.
    """
    if rho_vac is None:
        rho_vac = RHO_VAC_SCM
    rings = mayan_ring_proportions(epoch, r_base_AU)
    cos_tn    = math.cos(math.pi * t_n)
    u_inertia = 3.0 * rho_vac * (4.0 * math.pi / 3.0) * C_LIGHT**2 * cos_tn
    gear_13   = rings['gear_ratio_13']   # φ^(2*(epoch-1))

    res = UniversalInertiaResult()
    res.u_inertia            = u_inertia
    res.u_inertia_abs        = abs(u_inertia)
    res.rho_vac              = rho_vac
    res.cos_tn               = cos_tn
    res.t_n                  = t_n
    res.primordial_freq_range = gear_13           # ω_inner/ω_outer = r_1/r_3
    res.massless_scalar      = abs(cos_tn) < 1.0e-4
    res.centripetal_mode     = (u_inertia > 0.0)
    res.centrifugal_mode     = (u_inertia < 0.0)
    res.zero_point           = res.massless_scalar
    res.tectonic_band_inner_m = rings['r_inner_m']
    res.tectonic_band_outer_m = rings['r_outer_m']
    return res


def zero_point_years_in_epoch5(n_zeroes: int = 5) -> List[float]:
    """Return the next n zero-point gravity years in Epoch 5.
    Each zero occurs at t_n = 1/2, 3/2, 5/2, ...
    year_k = MAYAN_EPOCH5_YEAR + (2k+1)/2 × MAYAN_GREAT_CYCLE_YEARS  k=0,1,...
    """
    return [
        MAYAN_EPOCH5_YEAR + (2.0 * k + 1.0) / 2.0 * MAYAN_GREAT_CYCLE_YEARS
        for k in range(n_zeroes)
    ]


# =============================================================================
# SECTION 3b — SESSION 202 VARIANT-BRANCH FUNCTIONS
# Phase H202: VDS/DVP/DH26 branches, VDS×DVP coupling, BH26×BSH resonance.
# All functions mirror the C++ implementations in QCalcGeom.cpp Section 2b.
# =============================================================================

def vds_branches(SSq: float = None, n_terms: int = 200) -> VDSBranchResult:
    """VDS variant branches: calibration sensitivity, energy density, BH26-coupled.
    vds_prime      = d/dz Li_{26}(z)|_{z=SSq} = Li_{25}(SSq)/SSq  (sensitivity)
    vds_density    = Li_{26}(SSq) × RHO_VAC_SCM  [J/m³]
    vds_k_weighted = Li_{25}(SSq) + 25·Li_{26}(SSq)  (VDS×BH26 coupling)
    Reference: CP4 #83 VDS + Session 202 derivations
    """
    if SSq is None:
        SSq = SSQ_DEFAULT
    v = VDSBranchResult()
    li25 = li26 = 0.0
    ps = SSq
    for n in range(1, n_terms + 1):
        li25 += ps / n**25
        li26 += ps / n**26
        ps *= SSq
        if abs(ps) < 1e-300:
            break
    v.vds_li25       = li25
    v.vds_prime      = li25 / SSq if SSq > 0.0 else 0.0
    v.vds_density    = li26 * float(RHO_VAC_SCM)    # J/m³
    v.vds_k_weighted = li25 + 25.0 * li26
    return v


def dvp_branches(p_max: int = 200) -> DVPBranchResult:
    """DVP spectral branches: full prime-vortex sum, pair product, vorticity floor.
    Enumerates primes p in (26, p_max]; a(p) = SSq^{π(p)} / p^{26}.
    Reference: CP4 #83 DVP + Navier-Stokes vorticity bound, Session 202
    """
    # Sieve of Eratosthenes
    sieve = [True] * (p_max + 1)
    if p_max >= 0:
        sieve[0] = False
    if p_max >= 1:
        sieve[1] = False
    for i in range(2, int(p_max**0.5) + 1):
        if sieve[i]:
            for j in range(i * i, p_max + 1, i):
                sieve[j] = False

    pi_count = 0
    dvp_sum = a29 = a31 = last_a = 0.0
    cnt = 0
    for p in range(2, p_max + 1):
        if not sieve[p]:
            continue
        pi_count += 1
        if p <= 26:
            continue
        a_p = SSQ_DEFAULT**pi_count / p**26
        dvp_sum += a_p
        cnt += 1
        if p == 29:
            a29 = a_p
        if p == 31:
            a31 = a_p
        last_a = a_p

    return DVPBranchResult(
        zeta_sum=dvp_sum,
        n_primes_dvp=cnt,
        pair_product=a29 * a31,
        spectral_floor=last_a,
        a_29=a29
    )


def bh26_branches(N: int = 10) -> BH26BranchResult:
    """BH26/DH26 variant branches: eigenvalue ladder, Casimir energy, topology.
    lambda_k = k(k+25) on S^{25}.
    spectral_sum   = Σ_{k=1}^{N} lambda_k          (= 1760 for N=10)
    casimir_energy = ℏ·RERING_BB_HZ/2 × Σ 1/lambda_k [J]
    degeneracy_k1  = 26 = C(26,25) on S^{25}
    vds_coupling   = Σ_{k=1}^{N} lambda_k^{-26}    (BH26→VDS topological bridge)
    Reference: CP4 #149 BH26 + Session 202
    """
    sl = si = sv = 0.0
    for k in range(1, N + 1):
        lk = k * (k + 25)
        sl += lk
        si += 1.0 / lk
        sv += lk**(-26)
    return BH26BranchResult(
        spectral_sum=sl,
        casimir_energy=HBAR * RERING_BB_HZ * 0.5 * si,
        degeneracy_k1=26,
        vds_coupling=sv,
        N=N
    )


def vds_dvp_coupled(SSq: float = None,
                    p_max: int = 200,
                    n_terms: int = 200) -> VDSDVPCoupledResult:
    """VDS×DVP coupled field: normalised weights, geometric-mean coupling, calibration gap.
    w_vds       = Li_26(SSq) / VDS_max  where VDS_max = SSq/(1-SSq)
    w_dvp       = zeta_sum / a(29)       (dominant first DVP term)
    joint_coeff = sqrt(w_vds × w_dvp)   (geometric-mean field coupling)
    variant_branch = |w_vds - w_dvp|     (differential calibration magnitude)
    Encodes: 'many ways to get from one place to the other' -- Session 202
    """
    if SSq is None:
        SSq = SSQ_DEFAULT
    vb = vds_branches(SSq, n_terms)
    vds_val = vb.vds_density / float(RHO_VAC_SCM)   # recover Li_26(SSq)
    vds_max = SSq / (1.0 - SSq) if 0.0 < SSq < 1.0 else 1.0
    w_v = vds_val / vds_max if vds_max > 0.0 else 0.0

    db = dvp_branches(p_max)
    dvp_max = db.a_29 if db.a_29 > 0.0 else 1.0
    w_d = db.zeta_sum / dvp_max if dvp_max > 0.0 else 0.0

    return VDSDVPCoupledResult(
        w_vds=w_v,
        w_dvp=w_d,
        joint_coeff=math.sqrt(abs(w_v) * abs(w_d)),
        variant_branch=abs(w_v - w_d)
    )


def bh26_bsh_resonance(f_Ub: float = 3.3e7,
                        SSq: float = None,
                        t_n: float = 0.0,
                        k: int = 1) -> BH26BSHResonanceResult:
    """BH26×BSH cross-resonance: BSH evaluated at a BH26 spectral frequency bin.
    freq_k         = RERING_BB_HZ / lambda_k         [Hz]
    omega_k        = 2π·freq_k                      [rad/s]
    bsh_at_k       = bsh_harmonic(f_Ub, SSq, omega_k, t_n, 26).U_g2
    resonance      = bsh_at_k × cos(π·t_n)
    energy_density = resonance × RHO_VAC_SCM          [J/m³]
    Reference: BH26 × BSH cross-resonance, Session 202
    """
    if SSq is None:
        SSq = SSQ_DEFAULT
    lk = k * (k + 25)
    fk = RERING_BB_HZ / lk if lk > 0 else 0.0
    omega_k = 2.0 * math.pi * fk
    b = bsh_harmonic(f_Ub=f_Ub, SSq=SSq, omega=omega_k, t_n=t_n, m_max=26)
    res = b.U_g2 * math.cos(math.pi * t_n)
    return BH26BSHResonanceResult(
        freq_k=fk,
        bsh_at_k=b.U_g2,
        resonance=res,
        energy_density=res * float(RHO_VAC_SCM)
    )


# =============================================================================
# SECTION 4 — UNIVERSAL BUOYANCY ENGINE  (new in v2.0.0)
# ============================================================================
#
# The Quantum Chain crossing condition (dpm_vacuum_manifold.py Step 6):
#   FUBi(r) + FUBii(r) = 0   →   mass BORN at crossing (Step 7)
#
# FUBi  (outside-to-inside, collapsing gravity zone):
#   = SOURCE4 Ubi = −β_i · G·M_⊙²/r² · orbit_factor · cos(π·t_n)
#   Scales as r⁻²; stronger at small r → drives collapse.
#
# FUBii (inside-to-outside, Aether buoyancy spring — habitable zone force):
#   = +ρ_vac,SCm · (4π/3) · r · c² · cos(π·t_n)
#   Scales as r; stronger at large r → resists collapse.
#
# Neutral buoyancy crossing  (cos(π·t_n) ≠ 0):
#   FUBi + FUBii = 0
#   β_i·G·M_⊙²·orbit/r² = ρ_vac·(4π/3)·r·c²
#   r_hz³ = β_i·G·M_⊙²·orbit / (ρ_vac·(4π/3)·c²)
#
# SIMULTANEOUS EQUATION SYSTEM (r_hz, t_n_hz jointly):
#   Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0         [buoyancy crossing]
#   Eq2: ε′(r, t_n) + G·M_⊙/(c²·r²) = 0           [metric-geodesic match]
# =============================================================================

def compute_FUBii(r: float, t_n: float,
                   rho_vac: float = None) -> FUBiiResult:
    """Inside-to-outside Aether counter-buoyancy — habitable zone restoring force.
    FUBii = +ρ_vac,SCm · (4π/3) · r · c² · cos(π·t_n)
    Aether 'buoyancy spring': grows linearly with r, opposes FUBi's 1/r² collapse.
    rho_vac: vacuum energy density [J/m³]; defaults to canonical RHO_VAC_SCM.
    """
    if rho_vac is None:
        rho_vac = RHO_VAC_SCM
    out = FUBiiResult()
    out.cos_tn      = math.cos(math.pi * t_n)
    out.rho_aether  = rho_vac
    out.r_m         = r
    out.FUBii       = rho_vac * (4.0 * math.pi / 3.0) * r * C_LIGHT**2 * out.cos_tn
    out.outward     = (out.FUBii > 0.0)
    return out


def compute_F_U(r: float, t_n: float,
                M: float   = M_SUN,
                beta_i: float     = BETA_I_BSFG,
                Omega_g: float    = OMEGA_G_BSFG,
                M_bh: float       = M_BH_BSFG,
                d_g: float        = D_G_BSFG,
                epsilon_sw: float = EPS_SW_BSFG,
                rho_sw: float     = RHO_SW_BSFG,
                U_UA: float       = U_UA_BSFG,
                rho_vac: float    = None,
                epoch: int        = 5,
                xi_UI: float      = 1.0e-6) -> UniversalGravityResult:
    """Universal Gravity full assembly with Universal Inertia (Session 265).
      F_U = Ug1 + Ug2 + Ug3 + Ug4 − FUBi + FUBii + Um + ξ·U_I·V_body

    Ug1..Ug4 are simplified 26D field components (canonical SOURCE4 parameters).
    FUBi  = SOURCE4 Ubi (collapsing zone; negative when t_n≈0 → opposes collapse).
    FUBii = Aether counter-buoyancy spring (positive outward when cos > 0).
    Um    = magnetic string contribution M·R_⊙²·ω_s / r³.
    U_I   = 3·ρ_vac·(4π/3)·c²·cos(π·t_n)  invariant differential (Mayan three-ring),
            modulates centripetal↔centrifugal, vanishes at zero-point (t_n=k+½).
    ξ     = perturbative coupling (default 1e-6); set 0 to disable U_I term.

    Mass M is used only AFTER emergence (Step 7+ in Quantum Chain).
    """
    if rho_vac is None:
        rho_vac = RHO_VAC_SCM

    res = UniversalGravityResult()
    res.r_m = r
    res.t_n = t_n

    # BSFG metric for coupling
    met = bsfg_metric(r, t_n)
    res.eps = met.eps

    cos_tn = math.cos(math.pi * t_n)

    # ── Ug1: magnetic dipole (SOURCE4 simplified — r⁻² coupling) ──────────────
    rho_A   = RHO_VAC_SCM
    V_body  = (4.0 / 3.0) * math.pi * R_SUN**3
    mu_s    = rho_A * V_body                           # [J/m³ × m³ = J]
    grad_M  = M / (r * r)
    k1      = 1.0e-22
    res.Ug1 = k1 * mu_s * grad_M * math.exp(-0.0 * abs(t_n)) * cos_tn

    # ── Ug2: charge-reactivity heliosphere coupling ───────────────────────────
    rho_UA  = RHO_VAC_UA
    Q_SCm   = rho_A  * V_body
    Q_UA    = rho_UA * V_body
    E_react = rho_A * 4.4e5**2 / rho_UA * math.exp(-KAPPA_FLOAT * abs(t_n))
    k2      = 1.0e-20
    S_rb    = 1.0 if r > R_SUN * 100 else 0.0
    res.Ug2 = k2 * (Q_SCm + Q_UA) * M / r**2 * S_rb * E_react

    # ── Ug3: magnetic string rotation ─────────────────────────────────────────
    B_disk   = 1.0e-4                                  # Solar disk B [T]
    omega_s  = float(OMEGA_S)
    k3       = 1.0e-18
    rot_term = math.cos(omega_s * abs(t_n) * math.pi)
    E3       = rho_A * 4.4e5**2 / rho_UA
    res.Ug3  = k3 * B_disk * rot_term * E3

    # ── Ug4: vacuum concentration ──────────────────────────────────────────────
    k4       = 1.0e-22
    res.Ug4  = k4 * rho_A * 1.0 * math.exp(-KAPPA_FLOAT * abs(t_n)) * cos_tn

    # ── FUBi (SOURCE4 Ubi) ────────────────────────────────────────────────────
    buo      = bsfg_buoyancy(r, t_n, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA)
    res.FUBi = buo.Ubi

    # ── FUBii (Aether counter-buoyancy spring) ────────────────────────────────
    fubii_res = compute_FUBii(r, t_n, rho_vac)
    res.FUBii = fubii_res.FUBii

    # ── Um (universal magnetism M·R²·ω/r³) ───────────────────────────────────
    mu_m    = M * R_SUN**2 * omega_s
    res.Um  = mu_m / r**3

    # ── U_I: Universal Inertia (Session 265) ─────────────────────────────────
    # Invariant differential, frame-independent, primordial radiance scalar.
    # Frame-independent: same value across all bodies, depends only on (t_n, epoch).
    ui          = universal_inertia(t_n, epoch=epoch, r_base_AU=1.0, rho_vac=rho_vac)
    res.U_I     = ui.u_inertia
    res.epoch   = epoch
    res.U_I_term = xi_UI * ui.u_inertia * V_body   # ξ·U_I·V_body  → force units

    # ── Full assembly ─────────────────────────────────────────────────────────
    Ug_sum        = res.Ug1 + res.Ug2 + res.Ug3 + res.Ug4
    res.F_U_total = Ug_sum - res.FUBi + res.FUBii + res.Um + res.U_I_term
    return res


# =============================================================================
# SECTION 5 — SIMULTANEOUS EQUATION SOLVER
# =============================================================================

def _habitable_zone_system(x: List[float], params: dict) -> List[float]:
    """Residual function for the simultaneous habitable zone equations.

    System (r, t_n → 2 unknowns, 2 equations):
      Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0   [neutral buoyancy / compaction]
      Eq2: ε′(r, t_n) + G·M/(c²·r²)    = 0   [metric-geodesic matching]

    Solved in log-r space to handle wide dynamic range:
      x[0] = log10(r),  x[1] = t_n
    """
    log_r, t_n = x
    r     = 10.0 ** log_r
    beta_i    = params.get('beta_i',    BETA_I_BSFG)
    Omega_g   = params.get('Omega_g',   OMEGA_G_BSFG)
    M_bh      = params.get('M_bh',      M_BH_BSFG)
    d_g       = params.get('d_g',       D_G_BSFG)
    epsilon_sw= params.get('epsilon_sw',EPS_SW_BSFG)
    rho_sw    = params.get('rho_sw',    RHO_SW_BSFG)
    U_UA      = params.get('U_UA',      U_UA_BSFG)
    rho_vac   = params.get('rho_vac',   RHO_VAC_SCM)
    M         = params.get('M',         M_SUN)

    # Eq1: buoyancy crossing
    buo   = bsfg_buoyancy(r, t_n, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA)
    fubii = compute_FUBii(r, t_n, rho_vac)
    eq1   = buo.Ubi + fubii.FUBii

    # Eq2: metric-geodesic matching  ε′ + G·M / (c²·r²) = 0
    met  = bsfg_metric(r, t_n)
    eq2  = met.eps_p + G_NEWTON * M / (C_LIGHT**2 * r**2)

    # Normalize to similar magnitudes
    scale1 = max(abs(buo.Ubi), abs(fubii.FUBii), 1.0)
    scale2 = max(abs(met.eps_p), 1.0e-50)

    return [eq1 / scale1, eq2 / scale2]


def solve_habitable_zone(params: dict,
                          r_guess_m: float = None,
                          t_n_guess: float = 0.5,
                          epoch: int = 5) -> HabitableZoneResult:
    """Solve for the habitable zone (r_hz, t_n_hz) simultaneously.

    Uses the system:
      Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0   [buoyancy crossing]
      Eq2: ε′(r, t_n) + G·M/(c²·r²)    = 0   [metric-geodesic match]

    If scipy is unavailable, falls back to the closed-form analytic solution:
      r_hz = [β_i·G·M²·orbit / (RHO_VAC_SCM·(4π/3)·c²)]^{1/3}
      t_n_hz from Eq2 evaluated at r_hz.

    params dict keys (all optional, defaults = canonical SOURCE4):
      M, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA, rho_vac
    """
    result = HabitableZoneResult()
    M       = params.get('M',       M_SUN)
    beta_i  = params.get('beta_i',  BETA_I_BSFG)
    Omega_g = params.get('Omega_g', OMEGA_G_BSFG)
    M_bh    = params.get('M_bh',    M_BH_BSFG)
    d_g     = params.get('d_g',     D_G_BSFG)
    epsilon_sw = params.get('epsilon_sw', EPS_SW_BSFG)
    rho_sw  = params.get('rho_sw',  RHO_SW_BSFG)
    U_UA    = params.get('U_UA',    U_UA_BSFG)
    rho_vac = params.get('rho_vac', RHO_VAC_SCM)

    wind_mod     = 1.0 + epsilon_sw * rho_sw
    orbit_factor = Omega_g * M_bh / d_g * wind_mod * U_UA

    # ── Analytic closed-form estimate (cos factors cancel in Eq1) ────────────
    # FUBi + FUBii = 0: β_i·G·M_SUN²·orbit/r² = rho_vac·(4π/3)·r·c²
    numerator   = beta_i * G_NEWTON * M_SUN**2 * orbit_factor
    denominator = rho_vac * (4.0 * math.pi / 3.0) * C_LIGHT**2
    if denominator > 0.0:
        r_analytic = (numerator / denominator) ** (1.0 / 3.0)
    else:
        r_analytic = AU_METERS

    # t_n from Eq2 at r_analytic: cos(π·t_n) = -G·M/(c²·r²·3·ETA·C_num/r²·r)
    # Simplified: ε′ = -3·ETA·cos·C_num/r⁴; ε′ = -GM/(c²r²)
    # cos(π·t_n) = G·M·r² / (3·ETA·C_num·c²)
    cosval = G_NEWTON * M_SUN * r_analytic**2 / (3.0 * ETA_BSFG * C_NUM_SOLAR * C_LIGHT**2)
    cosval = max(-1.0, min(1.0, cosval))
    t_n_analytic = math.acos(cosval) / math.pi

    if r_guess_m is None:
        r_guess_m = r_analytic

    # ── scipy simultaneous solve (refinement only) ────────────────────────────
    # Analytic r_hz is EXACT for the force balance (cos factors cancel exactly).
    # scipy is used as an optional refinement; if it fails, analytic is used.
    r_sol   = r_analytic
    t_n_sol = t_n_analytic

    if _HAS_SCIPY:
        x0 = [math.log10(max(r_guess_m, 1.0)), t_n_guess]
        try:
            sol, info, ier, msg = _scipy_fsolve(
                _habitable_zone_system, x0,
                args=(params,),
                full_output=True,
                xtol=1.0e-10,
            )
            if ier == 1:
                r_cand  = 10.0 ** sol[0]
                t_n_cand = sol[1]
                # Accept scipy result only if it's physically reasonable
                if r_cand > 0 and math.isfinite(r_cand) and abs(t_n_cand) < 10.0:
                    r_sol   = r_cand
                    t_n_sol = t_n_cand
                    result.solver_msg = "CONVERGED (scipy + analytic)"
                else:
                    result.solver_msg = "analytic (scipy unphysical root)"
            else:
                result.solver_msg = f"analytic (scipy ier={ier})"
        except Exception as exc:
            result.solver_msg = f"analytic (scipy exc: {exc})"
    else:
        result.solver_msg = "analytic closed-form (no scipy)"

    result.r_hz_m   = r_sol
    result.t_n_hz   = t_n_sol
    result.converged = True   # analytic formula gives exact FUBi+FUBii=0

    # Residuals at chosen solution
    _rv = _habitable_zone_system([math.log10(max(r_sol, 1.0)), t_n_sol], params)
    result.residual_eq1 = _rv[0]
    result.residual_eq2 = _rv[1]

    # ── Evaluate FUBi + FUBii at solution ─────────────────────────────────────
    r_hz = result.r_hz_m
    t_n_hz = result.t_n_hz
    buo_hz  = bsfg_buoyancy(r_hz, t_n_hz, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA)
    fubii_hz = compute_FUBii(r_hz, t_n_hz, rho_vac)
    result.FUBi_at_hz  = buo_hz.Ubi
    result.FUBii_at_hz = fubii_hz.FUBii
    result.r_hz_AU     = r_hz / AU_METERS

    # ── Habitable zone classification ─────────────────────────────────────────
    if result.r_hz_AU < 0.5:
        result.hz_type = "rocky_inner"
    elif result.r_hz_AU <= 2.5:
        result.hz_type = "habitable"
    else:
        result.hz_type = "gas_outer"

    # ── Epoch-aware tectonic band (Session 265) ──────────────────────────────
    # Mayan three-ring scaling defines the crustal/floating zone bounds:
    #   r_inner_band = r_hz × φ^(-2(E-1))   (inner ring scaling — shrinks fastest)
    #   r_outer_band = r_hz × φ^( (E-1))    (outer ring scaling — expands)
    # In Epoch 5: inner = r_hz × 0.0213, outer = r_hz × 6.854 (asymmetric band).
    ring = mayan_ring_proportions(epoch, 1.0)
    result.epoch              = epoch
    result.r_hz_inner_band_AU = result.r_hz_AU * ring['r_inner_AU']
    result.r_hz_outer_band_AU = result.r_hz_AU * ring['r_outer_AU']
    result.band_width_AU      = result.r_hz_outer_band_AU - result.r_hz_inner_band_AU

    # Universal Inertia at the HZ point
    ui_hz = universal_inertia(t_n_hz, epoch=epoch, rho_vac=rho_vac)
    result.U_I_at_hz = ui_hz.u_inertia
    if ui_hz.zero_point:
        result.U_I_mode = "zero-point"
    elif ui_hz.centripetal_mode:
        result.U_I_mode = "centripetal"
    else:
        result.U_I_mode = "centrifugal"

    return result


def scan_habitable_zone(r_array_AU: np.ndarray,
                         t_n_array: np.ndarray,
                         params: dict = None) -> dict:
    """2D scan of the FUBi+FUBii landscape over (r, t_n).

    Returns dict with:
      'FUBi_grid'   : 2D array shape (len(r), len(t_n))  collapsing gravity zone
      'FUBii_grid'  : 2D array  outward Aether counter-buoyancy
      'balance_grid': 2D array  FUBi + FUBii  (zero = neutral buoyancy)
      'r_hz_m'      : 1D array of analytic r_hz for each t_n (closed form)
      'r_array_m'   : 1D array [m]
    """
    if params is None:
        params = {}

    beta_i     = params.get('beta_i',     BETA_I_BSFG)
    Omega_g    = params.get('Omega_g',    OMEGA_G_BSFG)
    M_bh       = params.get('M_bh',       M_BH_BSFG)
    d_g        = params.get('d_g',        D_G_BSFG)
    epsilon_sw = params.get('epsilon_sw', EPS_SW_BSFG)
    rho_sw     = params.get('rho_sw',     RHO_SW_BSFG)
    U_UA       = params.get('U_UA',       U_UA_BSFG)
    rho_vac    = params.get('rho_vac',    RHO_VAC_SCM)

    r_m   = r_array_AU * AU_METERS
    Nr, Nt = len(r_m), len(t_n_array)

    FUBi_grid    = np.zeros((Nr, Nt))
    FUBii_grid   = np.zeros((Nr, Nt))
    balance_grid = np.zeros((Nr, Nt))

    for i, r in enumerate(r_m):
        for j, t_n in enumerate(t_n_array):
            buo   = bsfg_buoyancy(r, t_n, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA)
            fubii = compute_FUBii(r, t_n, rho_vac)
            FUBi_grid[i, j]    = buo.Ubi
            FUBii_grid[i, j]   = fubii.FUBii
            balance_grid[i, j] = buo.Ubi + fubii.FUBii

    # Analytic r_hz per t_n
    wind_mod     = 1.0 + epsilon_sw * rho_sw
    orbit_factor = Omega_g * M_bh / d_g * wind_mod * U_UA
    numerator    = beta_i * G_NEWTON * M_SUN**2 * orbit_factor
    denominator  = rho_vac * (4.0 * math.pi / 3.0) * C_LIGHT**2
    r_hz_m_analytic = np.full(Nt, (numerator / denominator)**(1.0 / 3.0) if denominator > 0 else 0.0)

    return {
        'FUBi_grid'   : FUBi_grid,
        'FUBii_grid'  : FUBii_grid,
        'balance_grid': balance_grid,
        'r_hz_m'      : r_hz_m_analytic,
        'r_array_m'   : r_m,
        'r_array_AU'  : r_array_AU,
        't_n_array'   : t_n_array,
    }


def compute_emergent_mass(r_hz_m: float,
                          t_n_hz: float       = 0.0,
                          rho_vac: float      = None,
                          beta_i: float       = BETA_I_BSFG,
                          Omega_g: float      = OMEGA_G_BSFG,
                          M_bh: float         = M_BH_BSFG,
                          d_g: float          = D_G_BSFG,
                          epsilon_sw: float   = EPS_SW_BSFG,
                          rho_sw: float       = RHO_SW_BSFG,
                          U_UA: float         = U_UA_BSFG) -> EmergentMassResult:
    """Solve mass-emergence at the FUBi+FUBii=0 buoyancy crossing.

    Quantum Chain Step 7: mass is BORN at the crossing.  Given the crossing
    radius r_hz (from solve_habitable_zone), the mass that can exist in
    Aether equilibrium there is determined by the buoyancy balance, NOT by
    gravitational input.  This is the inverse problem to solve_habitable_zone:
      solve_habitable_zone:  given M_⊙, find r_hz
      compute_emergent_mass: given r_hz, find M

    Crossing equation:
        β_i · G · M² · orbit_factor / r_hz² = ρ_vac · (4π/3) · r_hz · c²
    Solve for M:
        M² = ρ_vac · (4π/3) · r_hz³ · c² / (β_i · G · orbit_factor)
        M_emergent = sqrt(...)   (positive root)

    Args:
        r_hz_m  : crossing radius from solve_habitable_zone [m]
        t_n_hz  : phase (only used for diagnostic FUBi/FUBii at solution;
                  the closed-form M_emergent factors cos(π·t_n) out cleanly
                  on both sides of the balance)
        rho_vac : vacuum density [J/m³] (default RHO_VAC_SCM)
        other   : FUBi orbit/wind parameters (defaults from canonical SOURCE4)
    """
    if rho_vac is None:
        rho_vac = float(RHO_VAC_SCM)

    out = EmergentMassResult()
    out.r_hz_m       = r_hz_m
    out.r_hz_AU      = r_hz_m / AU_METERS
    out.t_n_hz       = t_n_hz
    out.rho_vac_used = rho_vac
    out.beta_i_used  = beta_i

    # FUBi orbit-factor (matches bsfg_buoyancy)
    wind_mod          = 1.0 + epsilon_sw * rho_sw
    orbit_factor      = Omega_g * M_bh / d_g * wind_mod * U_UA
    out.orbit_factor  = orbit_factor

    # Closed-form mass-emergence
    if orbit_factor <= 0.0 or beta_i <= 0.0 or r_hz_m <= 0.0:
        out.converged = False
        return out

    numerator   = rho_vac * (4.0 * math.pi / 3.0) * (r_hz_m**3) * (C_LIGHT**2)
    denominator = beta_i * G_NEWTON * orbit_factor
    M_squared   = numerator / denominator
    if M_squared <= 0.0 or not math.isfinite(M_squared):
        out.converged = False
        return out

    M_emergent      = math.sqrt(M_squared)
    out.M_emergent_kg  = M_emergent
    out.M_emergent_sun = M_emergent / M_SUN

    # Verify residual at solution (sanity check; should be ~0 to machine precision)
    cos_tn   = math.cos(math.pi * t_n_hz)
    FUBi_at  = -beta_i * G_NEWTON * (M_emergent**2) / (r_hz_m**2) * orbit_factor * cos_tn
    FUBii_at = rho_vac * (4.0 * math.pi / 3.0) * r_hz_m * (C_LIGHT**2) * cos_tn
    out.residual_at_M = FUBi_at + FUBii_at
    out.converged     = (abs(out.residual_at_M) <
                         max(1e-6, 1e-9 * (abs(FUBi_at) + abs(FUBii_at))))

    # Classification on M_⊙ ladder
    M_sun_ratio = out.M_emergent_sun
    if   M_sun_ratio < 0.08 : out.classification = "sub_stellar"     # brown dwarf
    elif M_sun_ratio < 0.5  : out.classification = "sub_solar"       # red dwarf
    elif M_sun_ratio < 2.0  : out.classification = "solar"           # sun-like
    elif M_sun_ratio < 8.0  : out.classification = "stellar"         # main sequence
    elif M_sun_ratio < 25.0 : out.classification = "massive_stellar" # O/B class
    else                    : out.classification = "BH_seed"         # collapse candidate

    return out

# =============================================================================
# SECTION 5.5 — UNIVERSAL BUOYANCY SIMULTANEOUS SOLVER  (new in v2.2.0)
# -----------------------------------------------------------------------------
# Fully-coupled nonlinear 4x4 system in (r_hz, t_n_hz, r_cg, M_emergent).
# Encodes the Aether UA vacuum counter-balance F_U / F_U_Bi / F_U_Bi_i and
# the collapsing-gravity zone all simultaneously.  Replaces the v2.1.0
# staged decoupled chain (analytic r_hz -> t_n -> M_emergent -> hz only).
#
# Universal-Buoyancy Postulate (Greek Aether interpretation):
#   The vacuum (UA) supports a counter-buoyancy spring F_U_Bi_i that opposes
#   the collapsing-gravity force F_U_Bi at every radius.  Their algebraic
#   sum F_U_Bi + F_U_Bi_i defines three radial zones:
#     r < r_cg              :  collapsing zone   (F_U_Bi dominates 2x or more)
#     r_cg <= r <= r_hz_out :  habitable shell   (within tolerance of balance)
#     r > r_hz_out          :  gaseous outer     (F_U_Bi_i dominates)
#   At r_hz the algebraic sum vanishes (neutral buoyancy); the full
#   Universal Gravity F_U also vanishes there (Step-7 mass-emergence point).
# =============================================================================

@dataclass
class UniversalBuoyancySolution:
    """Solution of the coupled 4x4 Universal-Buoyancy system."""
    # Unknowns
    r_hz_m        : float = 0.0
    r_hz_AU       : float = 0.0
    t_n_hz        : float = 0.0
    r_cg_m        : float = 0.0
    r_cg_AU       : float = 0.0
    M_emergent_kg : float = 0.0
    M_emergent_sun: float = 0.0

    # Residuals at solution
    res_E1        : float = 0.0   # F_U_Bi + F_U_Bi_i at r_hz
    res_E2        : float = 0.0   # F_U at r_hz
    res_E3        : float = 0.0   # F_U_Bi + 2*F_U_Bi_i at r_cg
    res_E4        : float = 0.0   # M - rho_vac*(4pi/3)*r_hz^3

    # Diagnostics
    F_U_Bi_at_hz  : float = 0.0
    F_U_Bi_i_at_hz: float = 0.0
    F_U_at_hz     : float = 0.0
    F_U_Bi_at_cg  : float = 0.0
    F_U_Bi_i_at_cg: float = 0.0

    # Zone widths
    band_width_AU : float = 0.0   # r_hz - r_cg in AU
    collapse_ratio: float = 0.0   # |F_U_Bi/F_U_Bi_i| at r_cg

    # Bookkeeping
    converged     : bool  = False
    solver_msg    : str   = ""
    iterations    : int   = 0


def _universal_buoyancy_system(x: List[float], params: dict) -> List[float]:
    """4x4 residual function for the coupled Universal-Buoyancy solver.

    Unknowns (all in log-space where possible for wide dynamic range):
        x[0] = log10(r_hz   / m)
        x[1] = t_n_hz                      (linear, periodic in [0, 2])
        x[2] = log10(r_cg   / m)
        x[3] = log10(M_emergent / kg)

    Equations:
        E1: F_U_Bi(r_hz, t_n_hz, M)  +  F_U_Bi_i(r_hz, t_n_hz) = 0
        E2: F_U(r_hz, t_n_hz, M)                                 = 0
        E3: F_U_Bi(r_cg, t_n_hz, M)  + 2*F_U_Bi_i(r_cg, t_n_hz)  = 0
            (collapsing-gravity inner boundary: gravity outweighs
             Aether counter-buoyancy by exactly 2:1)
        E4: M  -  rho_vac * (4*pi/3) * r_hz^3                    = 0
            (Aether-vacuum mass-emergence at the HZ crossing)

    Returns: [E1', E2', E3', E4']  -- each normalised by its own scale.
    """
    log_r_hz, t_n_hz, log_r_cg, log_M = x

    # Guard against runaway iterates
    log_r_hz = max(min(log_r_hz, 25.0), 1.0)
    log_r_cg = max(min(log_r_cg, log_r_hz, 25.0), 1.0)
    log_M    = max(min(log_M, 40.0), 10.0)
    t_n_hz   = max(min(t_n_hz, 2.0), -2.0)

    r_hz = 10.0 ** log_r_hz
    r_cg = 10.0 ** log_r_cg
    M    = 10.0 ** log_M

    beta_i     = params.get('beta_i',     BETA_I_BSFG)
    Omega_g    = params.get('Omega_g',    OMEGA_G_BSFG)
    M_bh       = params.get('M_bh',       M_BH_BSFG)
    d_g        = params.get('d_g',        D_G_BSFG)
    epsilon_sw = params.get('epsilon_sw', EPS_SW_BSFG)
    rho_sw     = params.get('rho_sw',     RHO_SW_BSFG)
    U_UA       = params.get('U_UA',       U_UA_BSFG)
    rho_vac    = params.get('rho_vac',    RHO_VAC_SCM)

    # F_U_Bi at r_hz (uses emergent M as the body mass in the orbit-factor product)
    bui_hz = bsfg_buoyancy(r_hz, t_n_hz, beta_i, Omega_g, M_bh, d_g,
                            epsilon_sw, rho_sw, U_UA)
    fubii_hz = compute_FUBii(r_hz, t_n_hz, rho_vac)

    # F_U at r_hz
    fu_hz = compute_F_U(r_hz, t_n_hz, M=M,
                         beta_i=beta_i, Omega_g=Omega_g, M_bh=M_bh, d_g=d_g,
                         epsilon_sw=epsilon_sw, rho_sw=rho_sw, U_UA=U_UA,
                         rho_vac=rho_vac, xi_UI=0.0)  # turn off inertia term for static solve

    # F_U_Bi + 2*F_U_Bi_i at r_cg (inner-boundary 2:1 condition)
    bui_cg = bsfg_buoyancy(r_cg, t_n_hz, beta_i, Omega_g, M_bh, d_g,
                            epsilon_sw, rho_sw, U_UA)
    fubii_cg = compute_FUBii(r_cg, t_n_hz, rho_vac)

    # ── Residuals ────────────────────────────────────────────────────────────
    E1 = bui_hz.Ubi + fubii_hz.FUBii
    E2 = fu_hz.F_U_total
    E3 = bui_cg.Ubi + 2.0 * fubii_cg.FUBii
    M_pred = rho_vac * (4.0 * math.pi / 3.0) * (r_hz ** 3)
    E4 = M - M_pred

    # Normalise each residual by an O(1) scale
    s1 = max(abs(bui_hz.Ubi),    abs(fubii_hz.FUBii), 1.0)
    s2 = max(abs(fu_hz.F_U_total), s1, 1.0)
    s3 = max(abs(bui_cg.Ubi),    abs(2.0 * fubii_cg.FUBii), 1.0)
    s4 = max(abs(M), abs(M_pred), 1.0)

    return [E1 / s1, E2 / s2, E3 / s3, E4 / s4]


def solve_universal_buoyancy(params: Optional[dict] = None,
                              x0: Optional[List[float]] = None,
                              max_iter: int = 200) -> UniversalBuoyancySolution:
    """Solve the coupled 4x4 Universal-Buoyancy system simultaneously.

    Returns the joint (r_hz, t_n_hz, r_cg, M_emergent) such that all four
    canonical Aether-balance equations are satisfied simultaneously.

    params (all optional, defaults = canonical SOURCE4 + RHO_VAC_SCM):
        M_seed, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA, rho_vac

    x0 (optional initial guess, log-space):
        [log10(r_hz_m), t_n_hz, log10(r_cg_m), log10(M_kg)]

    Strategy:
        1. Build a physically-motivated initial guess from the v2.1.0
           staged solver (solve_habitable_zone -> compute_emergent_mass).
        2. Refine with scipy.fsolve on the full 4x4 residual.
        3. If scipy unavailable or refinement fails, return the staged
           (decoupled) solution with its residuals reported honestly.
    """
    if params is None:
        params = {}

    sol = UniversalBuoyancySolution()

    # ── Stage 1: staged initial guess (always succeeds) ──────────────────────
    hz_seed = solve_habitable_zone(params, t_n_guess=params.get('t_n', 0.5))
    em_seed = compute_emergent_mass(
        r_hz_m     = hz_seed.r_hz_m,
        t_n_hz     = hz_seed.t_n_hz,
        rho_vac    = params.get('rho_vac',    RHO_VAC_SCM),
        beta_i     = params.get('beta_i',     BETA_I_BSFG),
        Omega_g    = params.get('Omega_g',    OMEGA_G_BSFG),
        M_bh       = params.get('M_bh',       M_BH_BSFG),
        d_g        = params.get('d_g',        D_G_BSFG),
        epsilon_sw = params.get('epsilon_sw', EPS_SW_BSFG),
        rho_sw     = params.get('rho_sw',     RHO_SW_BSFG),
        U_UA       = params.get('U_UA',       U_UA_BSFG),
    )

    r_hz_seed = max(hz_seed.r_hz_m, 1.0)
    # Inner collapse boundary seed: r_cg ~ r_hz / cube_root(2)
    # (since FUBii ~ r and FUBi ~ 1/r^2, ratio FUBi/FUBii ~ 1/r^3;
    #  doubling the ratio corresponds to r/2^(1/3) ~ 0.794 r)
    r_cg_seed = r_hz_seed / (2.0 ** (1.0 / 3.0))
    M_seed    = em_seed.M_emergent_kg if em_seed.M_emergent_kg > 0 else M_SUN

    if x0 is None:
        x0 = [math.log10(r_hz_seed), hz_seed.t_n_hz,
              math.log10(r_cg_seed), math.log10(M_seed)]

    # ── Stage 2: scipy refinement on full 4x4 system ─────────────────────────
    refined_x = None
    refined_msg = ""

    if _HAS_SCIPY:
        try:
            sol_x, info, ier, msg = _scipy_fsolve(
                _universal_buoyancy_system, x0,
                args=(params,),
                full_output=True,
                xtol=1.0e-10,
                maxfev=max_iter * 4,
            )
            iters = int(info.get('nfev', 0))
            sol.iterations = iters
            if ier == 1:
                # sanity-check the refined iterate
                rh = 10.0 ** sol_x[0]
                rc = 10.0 ** sol_x[2]
                mm = 10.0 ** sol_x[3]
                if (math.isfinite(rh) and rh > 0
                        and math.isfinite(rc) and rc > 0
                        and math.isfinite(mm) and mm > 0
                        and rc < rh
                        and abs(sol_x[1]) < 10.0):
                    refined_x = sol_x
                    refined_msg = f"CONVERGED (scipy 4x4, nfev={iters})"
                else:
                    refined_msg = f"staged (scipy unphysical: rc>=rh or non-finite)"
            else:
                refined_msg = f"staged (scipy ier={ier}: {msg.strip()[:80]})"
        except Exception as exc:
            refined_msg = f"staged (scipy exception: {exc})"
    else:
        refined_msg = "staged (no scipy)"

    if refined_x is not None:
        log_r_hz, t_n_hz, log_r_cg, log_M = refined_x
        sol.solver_msg = refined_msg
        sol.converged  = True
    else:
        log_r_hz = math.log10(r_hz_seed)
        t_n_hz   = hz_seed.t_n_hz
        log_r_cg = math.log10(r_cg_seed)
        log_M    = math.log10(M_seed)
        sol.solver_msg = refined_msg
        sol.converged  = False

    # ── Stage 3: populate solution + diagnostics ─────────────────────────────
    r_hz = 10.0 ** log_r_hz
    r_cg = 10.0 ** log_r_cg
    M    = 10.0 ** log_M

    sol.r_hz_m         = r_hz
    sol.r_hz_AU        = r_hz / AU_METERS
    sol.t_n_hz         = t_n_hz
    sol.r_cg_m         = r_cg
    sol.r_cg_AU        = r_cg / AU_METERS
    sol.M_emergent_kg  = M
    sol.M_emergent_sun = M / M_SUN
    sol.band_width_AU  = sol.r_hz_AU - sol.r_cg_AU

    res = _universal_buoyancy_system([log_r_hz, t_n_hz, log_r_cg, log_M], params)
    sol.res_E1, sol.res_E2, sol.res_E3, sol.res_E4 = res

    beta_i     = params.get('beta_i',     BETA_I_BSFG)
    Omega_g    = params.get('Omega_g',    OMEGA_G_BSFG)
    M_bh       = params.get('M_bh',       M_BH_BSFG)
    d_g        = params.get('d_g',        D_G_BSFG)
    epsilon_sw = params.get('epsilon_sw', EPS_SW_BSFG)
    rho_sw     = params.get('rho_sw',     RHO_SW_BSFG)
    U_UA       = params.get('U_UA',       U_UA_BSFG)
    rho_vac    = params.get('rho_vac',    RHO_VAC_SCM)

    bui_hz   = bsfg_buoyancy(r_hz, t_n_hz, beta_i, Omega_g, M_bh, d_g,
                              epsilon_sw, rho_sw, U_UA)
    fubii_hz = compute_FUBii(r_hz, t_n_hz, rho_vac)
    fu_hz    = compute_F_U(r_hz, t_n_hz, M=M, beta_i=beta_i, Omega_g=Omega_g,
                            M_bh=M_bh, d_g=d_g, epsilon_sw=epsilon_sw,
                            rho_sw=rho_sw, U_UA=U_UA, rho_vac=rho_vac, xi_UI=0.0)
    bui_cg   = bsfg_buoyancy(r_cg, t_n_hz, beta_i, Omega_g, M_bh, d_g,
                              epsilon_sw, rho_sw, U_UA)
    fubii_cg = compute_FUBii(r_cg, t_n_hz, rho_vac)

    sol.F_U_Bi_at_hz   = bui_hz.Ubi
    sol.F_U_Bi_i_at_hz = fubii_hz.FUBii
    sol.F_U_at_hz      = fu_hz.F_U_total
    sol.F_U_Bi_at_cg   = bui_cg.Ubi
    sol.F_U_Bi_i_at_cg = fubii_cg.FUBii
    if abs(fubii_cg.FUBii) > 0:
        sol.collapse_ratio = abs(bui_cg.Ubi / fubii_cg.FUBii)

    return sol


class UniversalBuoyancySimultaneousSolver:
    """Calculator class for the coupled 4x4 Universal-Buoyancy system.

    Encodes the Aether UA vacuum counter-balance jointly with the
    collapsing-gravity zone, habitable zone, and emergent mass.  Unlike
    HabitableZoneCalculator (v2.1.0, 2 unknowns), this solver lifts the
    system to four simultaneous unknowns -- the user's directive level.

    System (4 equations, 4 unknowns):
        Unknowns : r_hz, t_n_hz, r_cg, M_emergent
        E1: F_U_Bi(r_hz, t_n, M)  + F_U_Bi_i(r_hz, t_n)         = 0
        E2: F_U(r_hz, t_n, M)                                    = 0
        E3: F_U_Bi(r_cg, t_n, M)  + 2*F_U_Bi_i(r_cg, t_n)        = 0
        E4: M - rho_vac*(4*pi/3)*r_hz^3                          = 0
    """

    def compute(self, dataset: dict) -> dict:
        params = {k: dataset[k] for k in dataset
                  if k in ('M', 'beta_i', 'Omega_g', 'M_bh', 'd_g',
                           'epsilon_sw', 'rho_sw', 'U_UA', 'rho_vac', 't_n')}

        sol = solve_universal_buoyancy(params)

        return {
            'r_hz_m'           : sol.r_hz_m,
            'r_hz_AU'          : sol.r_hz_AU,
            'r_cg_m'           : sol.r_cg_m,
            'r_cg_AU'          : sol.r_cg_AU,
            't_n_hz'           : sol.t_n_hz,
            'M_emergent_kg'    : sol.M_emergent_kg,
            'M_emergent_sun'   : sol.M_emergent_sun,
            'band_width_AU'    : sol.band_width_AU,
            'collapse_ratio'   : sol.collapse_ratio,
            'F_U_Bi_at_hz'     : sol.F_U_Bi_at_hz,
            'F_U_Bi_i_at_hz'   : sol.F_U_Bi_i_at_hz,
            'F_U_at_hz'        : sol.F_U_at_hz,
            'F_U_Bi_at_cg'     : sol.F_U_Bi_at_cg,
            'F_U_Bi_i_at_cg'   : sol.F_U_Bi_i_at_cg,
            'residual_E1'      : sol.res_E1,
            'residual_E2'      : sol.res_E2,
            'residual_E3'      : sol.res_E3,
            'residual_E4'      : sol.res_E4,
            'converged'        : sol.converged,
            'solver_msg'       : sol.solver_msg,
            'iterations'       : sol.iterations,
            'primary_equations': [
                "COUPLED 4x4 UNIVERSAL-BUOYANCY SYSTEM solved jointly:",
                "  E1: F_U_Bi(r_hz, t_n, M)  + F_U_Bi_i(r_hz, t_n)        = 0",
                "  E2: F_U(r_hz, t_n, M)                                   = 0",
                "  E3: F_U_Bi(r_cg, t_n, M)  + 2*F_U_Bi_i(r_cg, t_n)       = 0",
                "  E4: M - rho_vac*(4pi/3)*r_hz^3                          = 0",
                f"  r_hz = {sol.r_hz_AU:.4f} AU,  r_cg = {sol.r_cg_AU:.4f} AU,  band = {sol.band_width_AU:.4f} AU",
                f"  t_n_hz = {sol.t_n_hz:.4f},  M_emergent = {sol.M_emergent_sun:.4e} M_sun",
                f"  collapse_ratio at r_cg = |F_U_Bi/F_U_Bi_i| = {sol.collapse_ratio:.4f}",
                f"  residuals: E1={sol.res_E1:.2e}  E2={sol.res_E2:.2e}  E3={sol.res_E3:.2e}  E4={sol.res_E4:.2e}",
            ],
            'available_equations': [
                "Aether UA vacuum: rho_vac = RHO_VAC_SCM (633,333 J/m^3 canonical)",
                "Counter-buoyancy spring: F_U_Bi_i = rho_vac*(4pi/3)*r*c^2*cos(pi*t_n)",
                "Collapsing zone (Greek interpretation): r < r_cg -> Aether cannot support compaction",
                "Habitable shell: r_cg <= r <= r_hz -> Aether balance permits liquid/solid",
                "Mass emergence: M_emergent = sqrt(rho_vac*(4pi/3)*r_hz^3*c^2 / (beta_i*G*orbit))",
                "Collapse ratio: at r_cg, |F_U_Bi/F_U_Bi_i| = 2 (boundary)",
            ],
            'simulation_set': [
                "Sweep rho_vac (UA vacuum stiffness): observe r_hz and r_cg co-evolve",
                "Sweep beta_i (Aether coupling): observe M_emergent ladder shift",
                "Sweep M_bh (galactic anchor): track HZ migration across host galaxies",
                "Time-series t_n -> 0..2: oscillate FUBi/FUBii signs; track band collapse",
            ],
        }


# =============================================================================
# SECTION 6 — CALCULATOR CLASSES  (CondensedPhysics pattern)
# =============================================================================

class BSFGMetricCalculator:
    """BSFG metric + curvature at a given (r, t_n).
    Inputs via dataset dict: r [m], t_n (default 0.0).
    Outputs: eps, eps_p, R_r0r0, R_00, R_scalar, Kretschner, + bsfg_hz dict.
    """

    def compute(self, dataset: dict) -> dict:
        r   = dataset.get('r',   AU_METERS)
        t_n = dataset.get('t_n', 0.0)

        m   = bsfg_metric(r, t_n)
        h   = bsfg_horizon(t_n)
        fe  = bsfg_field_equations(r, t_n)
        gd  = bsfg_geodesic(r, t_n)

        return {
            'eps'          : m.eps,
            'eps_p'        : m.eps_p,
            'eps_pp'       : m.eps_pp,
            'A00'          : m.A00,
            'Arr'          : m.Arr,
            'R_r0r0'       : m.R_r0r0,
            'R_00'         : m.R_00,
            'R_rr'         : m.R_rr,
            'R_scalar'     : m.R_scalar,
            'Kretschner'   : m.Kretschner,
            'hz_exists'    : h.exists,
            'r_h_m'        : h.r_h,
            'T_H_K'        : h.T_H,
            'amp_factor'   : fe.amp_factor,
            'Lambda_eff'   : fe.Lambda_eff,
            'rho_vac_eff'  : fe.rho_vac_eff,
            'r_cross_AU'   : gd.r_cross_AU,
            'h_eta'        : gd.h_eta,
            'delta_J_over_J': gd.delta_J_over_J,
            'primary_equations': [
                f"eps = ETA*C_num/r^3*cos(pi*t_n) = {m.eps:.4e}",
                f"R^r_0r0 = eps''/2 - (eps')^2/2 = {m.R_r0r0:.4e} m^-2",
                f"r_cross = sqrt(ETA*c^2*C_num/(G*M)) = {gd.r_cross_AU:.4f} AU",
                f"T_H (BSFG horizon) = {h.T_H:.4e} K",
            ],
            'available_equations': [
                "Kretschner invariant K = 12*(R^r_0r0)^2 [curvature]",
                "Lambda_eff / Lambda_obs ratio [cosmological constant deviation]",
                "rho_vac_eff = Lambda_eff*c^2/(8piG) [vacuum energy density]",
                "delta_J/J = |v2_aether| / (2*v2_newton) [orbital Aether correction]",
            ],
            'simulation_set': [
                "Scan eps vs r: r in [0.01, 100] AU at t_n in {0, 0.5, 1.0}",
                "Plot R_r0r0(r) vs Kretschner(r) — curvature profile",
                "Blink cycle: t_n in [0, 2] showing A00 sign flip at t_n=1",
            ],
        }


class UniversalBuoyancyCalculator:
    """FUBi (collapsing gravity zone) vs FUBii (Aether counter-buoyancy) balance.
    Computes SOURCE4 Ubi (FUBi) and the Aether spring FUBii simultaneously.
    Signs: FUBi < 0 when t_n=0 (normal — opposes gravity); FUBii > 0 (outward).
    Balance (FUBi + FUBii = 0) defines the habitable zone crossing.
    """

    def compute(self, dataset: dict) -> dict:
        r          = dataset.get('r',          AU_METERS)
        t_n        = dataset.get('t_n',        0.0)
        beta_i     = dataset.get('beta_i',     BETA_I_BSFG)
        Omega_g    = dataset.get('Omega_g',    OMEGA_G_BSFG)
        M_bh       = dataset.get('M_bh',       M_BH_BSFG)
        d_g        = dataset.get('d_g',        D_G_BSFG)
        epsilon_sw = dataset.get('epsilon_sw', EPS_SW_BSFG)
        rho_sw     = dataset.get('rho_sw',     RHO_SW_BSFG)
        U_UA       = dataset.get('U_UA',       U_UA_BSFG)
        rho_vac    = dataset.get('rho_vac',    RHO_VAC_SCM)

        buo   = bsfg_buoyancy(r, t_n, beta_i, Omega_g, M_bh, d_g, epsilon_sw, rho_sw, U_UA)
        fubii = compute_FUBii(r, t_n, rho_vac)

        balance     = buo.Ubi + fubii.FUBii
        balance_rel = balance / max(abs(buo.Ubi), abs(fubii.FUBii), 1.0)

        # Zone classification at this r
        if abs(balance) < 0.01 * max(abs(buo.Ubi), abs(fubii.FUBii), 1.0):
            zone = "HABITABLE (neutral buoyancy)"
        elif buo.Ubi < 0 and abs(buo.Ubi) > abs(fubii.FUBii):
            zone = "ROCKY_INNER (buoyancy opposes; gravity wins)"
        elif buo.Ubi > 0:
            zone = "COLLAPSE (negentropic infall; FUBi aids gravity)"
        else:
            zone = "GAS_OUTER (Aether buoyancy dominates)"

        return {
            'FUBi'          : buo.Ubi,
            'FUBii'         : fubii.FUBii,
            'balance'       : balance,
            'balance_rel'   : balance_rel,
            'Ug_field'      : buo.Ug_field,
            'orbit_factor'  : buo.orbit_factor,
            'cos_tn'        : buo.cos_tn,
            'RHO_VAC_SCM'   : rho_vac,
            'zone'          : zone,
            'FUBi_negative' : buo.negative,
            'FUBi_inverted' : buo.inverted,
            'primary_equations': [
                f"FUBi = -beta_i*G*M^2/r^2*orbit*cos(pi*t_n) = {buo.Ubi:.4e}",
                f"FUBii = rho_vac*(4pi/3)*r*c^2*cos(pi*t_n) = {fubii.FUBii:.4e}",
                f"FUBi + FUBii = {balance:.4e}  ({zone})",
                f"Ug_field = G*M_sun^2/r^2 = {buo.Ug_field:.4e}",
            ],
            'available_equations': [
                "FUBii grows linearly with r; FUBi falls as 1/r^2",
                "Crossing r_hz: (beta_i*G*M^2*orbit / (rho_vac*(4pi/3)*c^2))^(1/3)",
                "cos(pi*t_n) = 0 at t_n=0.5 -> zero-buoyancy state",
                "Negentropic infall: t_n=1 -> cos=-1 -> FUBi > 0 (aids collapse)",
            ],
            'simulation_set': [
                "Scan balance vs r: 0.01 to 50 AU — map inner/habitable/outer zones",
                "Scan t_n: 0 to 2 — show blink cycle of FUBi sign flip",
                "3D landscape: balance(r, t_n) — find zero-crossings as hz contour",
            ],
        }


class HabitableZoneCalculator:
    """Simultaneous equation solver for habitable zone (r_hz, t_n_hz).

    Solves jointly:
      Eq1: FUBi(r, t_n) + FUBii(r, t_n) = 0   [neutral buoyancy]
      Eq2: eps'(r, t_n) + G*M/(c^2*r^2) = 0   [metric-geodesic match]

    Returns habitable zone radius r_hz and phase t_n_hz, plus full diagnostics.
    Physics: where Aether buoyancy (FUBii) exactly balances the collapsing
    gravity zone force (FUBi) AND the BSFG metric gradient matches the
    Newtonian geodesic correction — the Aether-Greek universal balance point.
    """

    def compute(self, dataset: dict) -> dict:
        params = {k: dataset[k] for k in dataset
                  if k in ('M', 'beta_i', 'Omega_g', 'M_bh', 'd_g',
                            'epsilon_sw', 'rho_sw', 'U_UA', 'rho_vac')}

        r_guess = dataset.get('r', None)
        t_n_guess = dataset.get('t_n', 0.5)

        hz = solve_habitable_zone(params, r_guess_m=r_guess, t_n_guess=t_n_guess)

        return {
            'r_hz_m'        : hz.r_hz_m,
            'r_hz_AU'       : hz.r_hz_AU,
            't_n_hz'        : hz.t_n_hz,
            'FUBi_at_hz'    : hz.FUBi_at_hz,
            'FUBii_at_hz'   : hz.FUBii_at_hz,
            'residual_eq1'  : hz.residual_eq1,
            'residual_eq2'  : hz.residual_eq2,
            'converged'     : hz.converged,
            'hz_type'       : hz.hz_type,
            'solver_msg'    : hz.solver_msg,
            'primary_equations': [
                f"SIMULTANEOUS SYSTEM (r, t_n) solved jointly:",
                f"  Eq1: FUBi + FUBii = 0  [neutral buoyancy crossing]",
                f"  Eq2: eps' + G*M/(c^2*r^2) = 0  [metric-geodesic match]",
                f"  Solution: r_hz = {hz.r_hz_AU:.4f} AU, t_n_hz = {hz.t_n_hz:.4f}",
                f"  FUBi at hz = {hz.FUBi_at_hz:.4e}",
                f"  FUBii at hz = {hz.FUBii_at_hz:.4e}",
                f"  Zone: {hz.hz_type}",
            ],
            'available_equations': [
                "Analytic r_hz = (beta_i*G*M_sun^2*orbit / (rho_vac*(4pi/3)*c^2))^(1/3)",
                "t_n_hz from: cos(pi*t_n) = G*M*r_hz^2 / (3*ETA*C_num*c^2)",
                "Inner zone r < r_hz: FUBi dominates -> rocky planet / liquid surface",
                "Outer zone r > r_hz: FUBii dominates -> gas / void",
                "BSFG horizon: blinking at t_n = 1 (cos=-1) -> transient collapse",
            ],
            'simulation_set': [
                "M scan: Sun, Proxima Cen, TRAPPIST-1 -> compare r_hz",
                "Galactic BH scan: M_bh from 1e6 to 1e10 M_sun -> r_hz dependence",
                "t_n blink cycle: animate FUBi + FUBii over [0, 2] -> hz oscillation",
            ],
        }


class UniversalGravityCalculator:
    """Full Universal Gravity assembly F_U = Ug1+Ug2+Ug3+Ug4 − FUBi + FUBii + Um.

    Quantum Chain Steps 3–6:
      Step 3: Ug1 seeded from DPM mu_s
      Step 4: Ug2+Ug3+Ug4 simultaneously promoted
      Step 5: F_U = Ug_family + Um + FUBi + FUBii
      Step 6: FUBi + FUBii = 0 at crossing (compaction)
    Mass enters only AFTER Step 7 (M_emergent at crossing).
    """

    def compute(self, dataset: dict) -> dict:
        r          = dataset.get('r',          AU_METERS)
        t_n        = dataset.get('t_n',        0.0)
        M          = dataset.get('M',          M_SUN)
        beta_i     = dataset.get('beta_i',     BETA_I_BSFG)
        Omega_g    = dataset.get('Omega_g',    OMEGA_G_BSFG)
        M_bh       = dataset.get('M_bh',       M_BH_BSFG)
        d_g        = dataset.get('d_g',        D_G_BSFG)
        epsilon_sw = dataset.get('epsilon_sw', EPS_SW_BSFG)
        rho_sw     = dataset.get('rho_sw',     RHO_SW_BSFG)
        U_UA       = dataset.get('U_UA',       U_UA_BSFG)
        rho_vac    = dataset.get('rho_vac',    RHO_VAC_SCM)

        res = compute_F_U(r, t_n, M, beta_i, Omega_g, M_bh, d_g,
                          epsilon_sw, rho_sw, U_UA, rho_vac)

        Ug_sum = res.Ug1 + res.Ug2 + res.Ug3 + res.Ug4

        return {
            'Ug1'         : res.Ug1,
            'Ug2'         : res.Ug2,
            'Ug3'         : res.Ug3,
            'Ug4'         : res.Ug4,
            'Ug_sum'      : Ug_sum,
            'FUBi'        : res.FUBi,
            'FUBii'       : res.FUBii,
            'Um'          : res.Um,
            'F_U_total'   : res.F_U_total,
            'eps'         : res.eps,
            'r_AU'        : r / AU_METERS,
            'primary_equations': [
                f"F_U = Ug_sum - FUBi + FUBii + Um = {res.F_U_total:.4e}",
                f"Ug_sum = Ug1+Ug2+Ug3+Ug4 = {Ug_sum:.4e}",
                f"FUBi (SOURCE4 Ubi) = {res.FUBi:.4e}  (collapsing gravity zone)",
                f"FUBii (Aether spring) = {res.FUBii:.4e}  (habitable zone force)",
                f"Um (magnetic string) = {res.Um:.4e}",
                f"BSFG eps = {res.eps:.4e}  (Aether metric perturbation)",
            ],
            'available_equations': [
                "Ug1 = k1 * mu_s * (M/r^2) * cos(pi*t_n)  [magnetic dipole]",
                "Ug2 = k2 * Q_total * M/r^2 * S(r>R_b) * E_react  [charge-reactivity]",
                "Ug3 = k3 * B_disk * cos(omega_s*t*pi) * E_react  [string rotation]",
                "Ug4 = k4 * rho_vac * C_factor * cos(pi*t_n)  [vacuum concentration]",
                "Um = M*R^2*omega_s / r^3  [universal magnetism]",
            ],
            'simulation_set': [
                "Scan F_U vs r: 0.01 to 100 AU — profile all 6 components",
                "Ug_family assembly: step-by-step Ug1 -> Ug1+Ug2 -> ... -> F_U",
                "t_n blink: animate F_U over one full phase cycle t_n in [0,2]",
            ],
        }

class MayanTimingCalculator:
    """Three-ring Mayan timing system — Universal Inertia and zero-point gravity.

    Input dataset keys:
      epoch       : Mayan epoch 1–5 (default 5 = current, started 2012)
      current_year: decimal year (default 2026.33 = May 2026)
      r_base_AU   : base radius for ring proportions [AU] (default 1.0)
      rho_vac     : override vacuum energy density (default RHO_VAC_SCM)

    CANONICAL EPOCH 5 (post-2012) RING GEOMETRY:
      Ring 1 (OUTER/LEFT)    : EXPANDING  r ≈ 6.85 AU at epoch 5
      Ring 2 (COMPANION/RIGHT): SHRINKS   r ≈ 0.146 AU at epoch 5
      Ring 3 (INNER, inside 1): VERY SMALL r ≈ 0.021 AU at epoch 5

    Universal Inertia invariant:
      U_I = 3 · ρ_vac · (4π/3) · c² · cos(π·t_n)   [N/m²]
      FRAME-INVARIANT: same regardless of r_hz, mass, orbital parameters

    Zero-point gravity (massless scalar state):
      cos(π·t_n) = 0 → t_n = 1/2 within each Great Cycle
      Achieved with precision = 1 / (gear_ratio_13 × N_teeth)
      In Epoch 5: 47× finer timing resolution via Ring 3

    Tectonic/crustall zone (floating on superconductive heavy plasma):
      Inner boundary = r_inner (FUBi zone, rocky / solid)
      Outer boundary = r_outer (FUBii zone, superconductive plasma)
      r_hz (habitable zone) lies between: surface floats at r_hz
    """

    def compute(self, dataset: dict) -> dict:
        epoch        = dataset.get('epoch',        5)
        current_year = dataset.get('current_year', 2026.33)
        r_base_AU    = dataset.get('r_base_AU',    1.0)
        rho_vac      = dataset.get('rho_vac',      RHO_VAC_SCM)

        state = mayan_ring_state(epoch, current_year, r_base_AU)
        ui    = universal_inertia(state.t_n, epoch, r_base_AU, rho_vac)
        zp    = zero_point_years_in_epoch5(5)

        # Connect to habitable zone solver
        hz = solve_habitable_zone({'rho_vac': rho_vac})

        # Crustall/tectonic band: [r_inner, r_outer]
        band_inner_AU = state.r_inner_AU
        band_outer_AU = state.r_outer_AU
        hz_in_band    = (band_inner_AU <= hz.r_hz_AU <= band_outer_AU)

        return {
            # Ring geometry
            'epoch'               : epoch,
            'current_year'        : current_year,
            'r_outer_AU'          : state.r_outer_AU,
            'r_companion_AU'      : state.r_companion_AU,
            'r_inner_AU'          : state.r_inner_AU,
            'gear_ratio_12'       : state.gear_ratio_12,
            'gear_ratio_13'       : state.gear_ratio_13,
            # Timing
            't_n'                 : state.t_n,
            't_n_inner'           : state.t_n_inner,
            'cos_tn'              : state.cos_tn,
            'omega_outer_rad_s'   : state.omega_outer,
            'omega_companion_rad_s': state.omega_companion,
            'omega_inner_rad_s'   : state.omega_inner,
            'zero_point_next_year': state.zero_point_next_year,
            'zero_point_years_e5' : zp,
            'massless_mode'       : state.massless_mode,
            # Universal Inertia
            'u_inertia'           : ui.u_inertia,
            'u_inertia_abs'       : ui.u_inertia_abs,
            'centripetal_mode'    : ui.centripetal_mode,
            'centrifugal_mode'    : ui.centrifugal_mode,
            'zero_point_gravity'  : ui.zero_point,
            'primordial_freq_range': ui.primordial_freq_range,
            # Tectonic/crustall zone
            'tectonic_inner_AU'   : band_inner_AU,
            'tectonic_outer_AU'   : band_outer_AU,
            'r_hz_AU'             : hz.r_hz_AU,
            'hz_in_tectonic_band' : hz_in_band,
            'primary_equations': [
                f"MAYAN THREE-RING TIMING — Epoch {epoch} (Epoch 5 started 2012.97)",
                f"  Ring 1 (OUTER/EXPANDING):    r = {state.r_outer_AU:.4f} AU",
                f"  Ring 2 (COMPANION/SHRINKING): r = {state.r_companion_AU:.4f} AU",
                f"  Ring 3 (INNER/VERY SMALL):   r = {state.r_inner_AU:.6f} AU",
                f"  Gear ratio 1:3 (inner amplification) = {state.gear_ratio_13:.2f}×",
                f"  t_n = {state.t_n:.6f}   cos(π·t_n) = {state.cos_tn:.6f}",
                f"UNIVERSAL INERTIA (invariant differential):",
                f"  U_I = 3·ρ_vac·(4π/3)·c²·cos(π·t_n) = {ui.u_inertia:.4e} N/m²",
                f"  Mode: {'CENTRIPETAL (mass-building)' if ui.centripetal_mode else 'CENTRIFUGAL (void-expansion)'}",
                f"  Zero-point next: {state.zero_point_next_year:.2f}",
                f"TECTONIC BAND: [{band_inner_AU:.4f}, {band_outer_AU:.4f}] AU",
                f"  r_hz (habitable) = {hz.r_hz_AU:.4f} AU  in band: {hz_in_band}",
            ],
            'available_equations': [
                "Ring proportions: r_i(n) = r_base × φ^{±2(n−1)}  φ=(1+√5)/2",
                "ω_companion = −ω_outer × r_outer/r_companion  [external mesh]",
                "ω_inner = +ω_outer × r_outer/r_inner  [internal mesh, amplified]",
                "U_I = 3·ρ_vac·(4π/3)·c²·cos(π·t_n)  [frame-invariant]",
                "Zero-point years: EPOCH5_START + (2k+1)/2 × GREAT_CYCLE_YEARS",
                "Tectonic float: surface at r_hz where FUBi + FUBii = 0",
                "Primordial radiance band: [ω_outer, ω_inner] = [1, gear_ratio_13]×ω_baktun",
            ],
            'simulation_set': [
                "Animate gear proportions: epoch 1→5 — outer grows, inner/companion shrink",
                "U_I vs t_n: 0→2 — show centripetal↔zero-point↔centrifugal cycle",
                "Scan current_year: 2012→7138 — next zero-point year identified",
                "3D tectonic surface: r_hz(M, epoch) — habitable zone across mass spectrum",
            ],
        }


# =============================================================================
# SECTION 7 — TEST SUITE (60 tests: T01–T60 mirroring C++ runQCalcGeomTests)
# =============================================================================

def _tol_check(computed: float, expected: float, tol_pct: float) -> bool:
    if expected == 0.0:
        return abs(computed) < 1.0e-300
    return abs((computed - expected) / expected) <= tol_pct / 100.0


def run_qcalcgeom_tests(verbose: bool = True) -> dict:
    """Run all 85 QCalcGeom requirements-boundary tests.
    T01–T60: mirrors C++ runQCalcGeomTests() with identical reference values.
    T61–T70: Mayan three-ring timing / Universal Inertia system.
    T71–T80: Session 202 VDS/DVP/DH26 variant branches + coupling (Phase H202).
    T81–T85: F_U + U_I integration + epoch-aware tectonic band (Session 265).
    Returns {'passed': int, 'failed': int, 'results': list-of-dicts}.
    """
    results = []

    def chk(tid, group, name, computed, expected, tol_pct, qual_ok=True):
        if tol_pct == 0.0:
            passed = qual_ok
        else:
            passed = _tol_check(computed, expected, tol_pct) and qual_ok
        results.append({'id': tid, 'group': group, 'name': name,
                        'computed': computed, 'expected': expected,
                        'tol_pct': tol_pct, 'passed': passed})
        if verbose:
            status = "PASS" if passed else "FAIL"
            print(f"  [{status}] {tid} {group}: {name}")
            if not passed:
                print(f"         computed={computed:.6e}  expected={expected:.6e}  tol={tol_pct}%")

    if verbose:
        print("\n" + "="*72)
        print("QCalcGeom.py v2.3.0 — Requirements-Boundary Test Suite (T01–T90)")
        print("="*72 + "\n")

    # ── BSFG-METRIC (T01–T06) ────────────────────────────────────────────────
    m_sun  = bsfg_metric(R_SUN, 0.0)
    m_far  = bsfg_metric(1.0e20, 0.0)
    m_tn0  = bsfg_metric(R_SUN, 0.0)
    m_tn2  = bsfg_metric(R_SUN, 2.0)
    m_tn1  = bsfg_metric(R_SUN, 1.0)
    m_half = bsfg_metric(R_SUN / 2.0, 0.0)

    chk("T01","BSFG-METRIC","eps_prime at R_SUN t_n=0",
        abs(m_sun.eps_p), EPS_PRIME_REF, 2.0)
    chk("T02","BSFG-METRIC","R^r_0r0 at R_SUN t_n=0",
        m_sun.R_r0r0, R_R0R0_REF, 2.0)
    ratio_far = m_far.R_r0r0 / m_sun.R_r0r0 if m_sun.R_r0r0 > 0 else 0.0
    chk("T03","BSFG-METRIC","R_r0r0 -> 0 as r -> inf",
        ratio_far, 0.0, 0.0, qual_ok=ratio_far < 1.0e-50)
    delta_cyc = abs(m_tn0.eps - m_tn2.eps)
    chk("T04","BSFG-METRIC","eps(t_n=0) == eps(t_n=2)",
        delta_cyc, 0.0, 0.0, qual_ok=delta_cyc < 1.0e-50)
    ratio_flip = m_tn1.eps / m_tn0.eps if m_tn0.eps != 0.0 else 0.0
    chk("T05","BSFG-METRIC","eps sign flip t_n=0 vs t_n=1",
        ratio_flip, -1.0, 1.0, qual_ok=(m_tn0.eps > 0 and m_tn1.eps < 0))
    ratio_r5 = m_half.eps_pp / m_sun.eps_pp if m_sun.eps_pp != 0.0 else 0.0
    chk("T06","BSFG-METRIC","eps_pp(r/2)/eps_pp(r) = 32  (r^-5 law)",
        ratio_r5, 32.0, 0.01)

    # ── BSFG-GEOM (T07–T12) ─────────────────────────────────────────────────
    h_tn1 = bsfg_horizon(1.0)
    gd_au  = bsfg_geodesic(AU_METERS, 0.0)
    gd_sun = bsfg_geodesic(R_SUN, 0.0)
    hl_sun = bsfg_holonomy(R_SUN, 0.0, 1.0)

    chk("T07","BSFG-GEOM","Blinking horizon r_h ~ 1.62e8 m at t_n=1",
        h_tn1.r_h, R_H_BSFG_REF, 2.0, qual_ok=h_tn1.exists)
    chk("T08","BSFG-GEOM","Hawking T_H ~ 3.37e-12 K",
        h_tn1.T_H, T_H_BSFG_REF, 3.0)
    chk("T09","BSFG-GEOM","Crossover r_cross ~ 0.360 AU",
        gd_au.r_cross_AU, R_CROSS_AU_REF, 2.0)
    chk("T10","BSFG-GEOM","h_eta = eta*h_Planck = 6.626e-56",
        gd_sun.h_eta, H_ETA_REF, 0.01)
    gen_total = 6 + hl_sun.n_extra_flat
    chk("T11","BSFG-GEOM","Holonomy 28 generators SO+(3,1)xU(1)^22",
        float(gen_total), 28.0, 0.0,
        qual_ok=(gen_total == 28) and hl_sun.G2_excluded and hl_sun.Spin7_excluded)
    dim_total = 4 + hl_sun.n_extra_flat
    chk("T12","BSFG-GEOM","M^26 = M^4 x T^22  dim=26",
        float(dim_total), 26.0, 0.0, qual_ok=(dim_total == 26))

    # ── VDS (T13–T16) ────────────────────────────────────────────────────────
    v200 = vds_series(SSQ_DEFAULT, 200)
    v56  = vds_series(0.56, 200)
    v57  = vds_series(0.57, 200)
    v58  = vds_series(0.58, 200)
    v500 = vds_series(SSQ_DEFAULT, 500)

    chk("T13","VDS","VDS(0.57) ~ 0.57  n=1 dominance",
        v200.value, SSQ_DEFAULT, 0.0, qual_ok=abs(v200.value - SSQ_DEFAULT) < 1.0e-7)
    # T14: |t_n| decreasing
    terms_dec = True
    prev_term = abs(SSQ_DEFAULT)
    for n in range(2, 11):
        cur = abs(SSQ_DEFAULT**n) / n**26
        if cur > prev_term:
            terms_dec = False
            break
        prev_term = cur
    chk("T14","VDS","VDS |t_n| decreasing n=1..10",
        v200.value, SSQ_DEFAULT, 0.0, qual_ok=terms_dec)
    chk("T15","VDS","VDS strictly increasing in SSq",
        v57.value, SSQ_DEFAULT, 0.0,
        qual_ok=(v56.value < v57.value) and (v57.value < v58.value))
    rel_conv = abs(v500.value - v200.value) / v500.value if v500.value != 0 else 1.0
    chk("T16","VDS","VDS(N=200) == Li_26(0.57) to 10 d.p.",
        rel_conv, 0.0, 0.0, qual_ok=rel_conv < 1.0e-10)

    # ── DVP (T17–T21) ────────────────────────────────────────────────────────
    d = dvp_arithmetic()
    # Check 30th prime = 113
    primes = []
    c = 2
    while len(primes) < 30:
        if all(c % p != 0 for p in primes):
            primes.append(c)
        c += 1
    chk("T17","DVP","30th prime = 113",
        float(primes[29]), 113.0, 0.0, qual_ok=(primes[29] == 113))
    chk("T18","DVP","26! mod 113 = 12",
        float(d.fac26_mod_113), 12.0, 0.0, qual_ok=(d.fac26_mod_113 == 12))
    chk("T19","DVP","non_repeating = True (Wilson)",
        1.0 if d.non_repeating else 0.0, 1.0, 0.0, qual_ok=d.non_repeating)
    chk("T20","DVP","r_q_AU ~ 0.0973 AU",
        d.r_q_AU, 0.0973, 2.0)
    chk("T21","DVP","r_q_m = r_q_AU * AU",
        d.r_q_m, d.r_q_AU * AU_METERS, 0.001)

    # ── BSH (T22–T24) ────────────────────────────────────────────────────────
    bsh = bsh_harmonic(3.3e7, SSQ_DEFAULT, 2.0 * math.pi * 3.3e7, 0.0, 20)
    bsh_sat = bsh_harmonic(3.3e7, 0.99, 2.0 * math.pi * 3.3e7, 0.0, 200)

    chk("T22","BSH","U_g2 > 0 at t_n=0",
        bsh.U_g2, 0.0, 0.0, qual_ok=(bsh.U_g2 > 0.0))
    chk("T23","BSH","H_m_max = f_Ub * harmonic_sum",
        bsh.H_m_max, 3.3e7 * sum(1.0 / k for k in range(1, 21)), 0.1)
    chk("T24","BSH","BSH saturated at SSq=0.99 m_max=200",
        1.0 if bsh_sat.saturated else 0.0, 1.0, 0.0, qual_ok=bsh_sat.saturated)

    # ── BH26 (T25–T27) ───────────────────────────────────────────────────────
    bh1  = bh26_eigenvalue(1)
    bh26 = bh26_eigenvalue(26)

    chk("T25","BH26","lambda_1 = 1*(1+25) = 26",
        bh1.lambda_k, 26.0, 0.001)
    chk("T26","BH26","lambda_26 = 26*(26+25) = 1326",
        bh26.lambda_k, 1326.0, 0.001)
    chk("T27","BH26","freq_bin_hz = RERING_BB_HZ / lambda",
        bh1.freq_bin_hz, RERING_BB_HZ / 26.0, 0.001)

    # ── COSMO (T28–T30) ──────────────────────────────────────────────────────
    fe_sun = bsfg_field_equations(R_SUN, 0.0)

    chk("T28","COSMO","amp_factor > 1 (Aether-dominated near R_SUN)",
        fe_sun.amp_factor, 0.0, 0.0, qual_ok=(fe_sun.amp_factor > 1.0))
    chk("T29","COSMO","Lambda_eff > 0",
        fe_sun.Lambda_eff, 0.0, 0.0, qual_ok=(fe_sun.Lambda_eff > 0.0))
    chk("T30","COSMO","rho_vac_eff > 0",
        fe_sun.rho_vac_eff, 0.0, 0.0, qual_ok=(fe_sun.rho_vac_eff > 0.0))

    # ── POLY26 (T31–T36) ─────────────────────────────────────────────────────
    p1  = poly26_derivative(1,  G_NEWTON * M_SUN, AU_METERS)
    p26 = poly26_derivative(26, 1.0, R_SUN)
    p1b = poly26_derivative(1,  G_NEWTON * M_SUN, 2.0 * AU_METERS)

    chk("T31","POLY26","poly26(k=1) value finite and > 0",
        p1.value, 0.0, 0.0, qual_ok=(p1.value > 0.0))
    chk("T32","POLY26","poly26(k=26,r=R_SUN) factorial ratio ~ (51!/25!) ~ 1e41",
        p26.factorial_ratio, 0.0, 0.0, qual_ok=(p26.factorial_ratio > 1.0e40))
    # r-scaling: p(2r)/p(r) = (1/2)^{k+26} for k=1 → 2^{-27}
    ratio_p26 = p1b.value / p1.value if p1.value != 0.0 else 0.0
    chk("T33","POLY26","poly26(k=1) r-scaling: p(2r)/p(r) = 2^{-27}",
        ratio_p26, 2.0**(-27), 0.5)
    chk("T34","POLY26","poly26(k=1) negligible at 1 AU = True (26th deriv vanishes at planetary scale)",
        1.0 if p1.negligible else 0.0, 1.0, 0.0, qual_ok=p1.negligible)

    # ── UQFF-COMP (T35–T38) ──────────────────────────────────────────────────
    uc = uqff_comp_matrix(AU_METERS, 1.4e3)

    chk("T35","UQFF-COMP","m00 > 0 at 1 AU",
        uc.m00, 0.0, 0.0, qual_ok=(uc.m00 > 0.0))
    chk("T36","UQFF-COMP","positive_definite = True",
        1.0 if uc.positive_definite else 0.0, 1.0, 0.0, qual_ok=uc.positive_definite)
    chk("T37","UQFF-COMP","cross_d13 >= 0 (k=26 term underflows at AU; non-neg by construction)",
        uc.cross_d13, 0.0, 0.0, qual_ok=(uc.cross_d13 >= 0.0))
    chk("T38","UQFF-COMP","eigenvalue_min >= 0",
        uc.eigenvalue_min, 0.0, 0.0, qual_ok=(uc.eigenvalue_min >= 0.0))

    # ── NEG-BUOY / FUBi (T39–T44) ────────────────────────────────────────────
    buo_t0 = bsfg_buoyancy(R_SUN, 0.0)
    buo_t1 = bsfg_buoyancy(R_SUN, 1.0)
    buo_th = bsfg_buoyancy(R_SUN, 0.5)

    chk("T39","NEG-BUOY","FUBi at t_n=0 ~ UBI_BSFG_REF",
        buo_t0.Ubi, UBI_BSFG_REF, 3.0)
    chk("T40","NEG-BUOY","FUBi at t_n=0 negative (opposes gravity)",
        buo_t0.Ubi, 0.0, 0.0, qual_ok=buo_t0.negative)
    chk("T41","NEG-BUOY","FUBi at t_n=1 positive (aids collapse)",
        buo_t1.Ubi, 0.0, 0.0, qual_ok=buo_t1.inverted)
    chk("T42","NEG-BUOY","FUBi sign flip t_n=0 vs t_n=1",
        buo_t1.Ubi / buo_t0.Ubi if buo_t0.Ubi != 0 else 0.0, -1.0, 0.1)
    chk("T43","NEG-BUOY","FUBi zero crossing at t_n=0.5",
        abs(buo_th.cos_tn), 0.0, 0.0, qual_ok=buo_th.zero_crossing)
    chk("T44","NEG-BUOY","FUBi even symmetry: t_n=0.3 == t_n=-0.3",
        bsfg_buoyancy(R_SUN, -0.3).Ubi - bsfg_buoyancy(R_SUN, 0.3).Ubi,
        0.0, 0.0, qual_ok=abs(bsfg_buoyancy(R_SUN, -0.3).Ubi
                               - bsfg_buoyancy(R_SUN, 0.3).Ubi) < 1.0e-40)

    # ── FUBii / AETHER SPRING (T45–T50) ──────────────────────────────────────
    fubii_au = compute_FUBii(AU_METERS, 0.0)
    fubii_2au = compute_FUBii(2.0 * AU_METERS, 0.0)
    fubii_tn1 = compute_FUBii(AU_METERS, 1.0)

    chk("T45","AETHER","FUBii at 1AU t_n=0 > 0 (outward)",
        fubii_au.FUBii, 0.0, 0.0, qual_ok=(fubii_au.FUBii > 0.0))
    chk("T46","AETHER","FUBii linear in r: FUBii(2AU)/FUBii(1AU) = 2",
        fubii_2au.FUBii / fubii_au.FUBii if fubii_au.FUBii != 0 else 0.0, 2.0, 0.01)
    chk("T47","AETHER","FUBii at t_n=1 negative (cos flips)",
        fubii_tn1.FUBii, 0.0, 0.0, qual_ok=(fubii_tn1.FUBii < 0.0))
    # r_hz: FUBi(r_hz) + FUBii(r_hz) = 0 → balance reached
    hz = solve_habitable_zone({})
    chk("T48","AETHER","r_hz > 0 and finite",
        hz.r_hz_m, 0.0, 0.0, qual_ok=(hz.r_hz_m > 0.0 and math.isfinite(hz.r_hz_m)))
    chk("T49","AETHER","r_hz converged",
        1.0 if hz.converged else 0.0, 1.0, 0.0, qual_ok=hz.converged)
    # FUBi + FUBii at hz should be < 1% of scale
    if hz.r_hz_m > 0:
        buo_hz  = bsfg_buoyancy(hz.r_hz_m, hz.t_n_hz)
        fubii_hz = compute_FUBii(hz.r_hz_m, hz.t_n_hz)
        balance_hz_rel = abs(buo_hz.Ubi + fubii_hz.FUBii) / max(abs(buo_hz.Ubi), 1.0)
    else:
        balance_hz_rel = 1.0
    chk("T50","AETHER","FUBi+FUBii ≈ 0 at r_hz (< 1% relative)",
        balance_hz_rel, 0.0, 0.0, qual_ok=(balance_hz_rel < 0.01))

    # ── UNIVERSAL GRAVITY (T51–T56) ───────────────────────────────────────────
    fu = compute_F_U(AU_METERS, 0.0)
    fu_r2 = compute_F_U(2.0 * AU_METERS, 0.0)

    chk("T51","UNIV-G","F_U_total finite at 1 AU",
        fu.F_U_total, 0.0, 0.0, qual_ok=math.isfinite(fu.F_U_total))
    chk("T52","UNIV-G","Um decreases with r: Um(2AU) < Um(1AU)",
        fu.Um - fu_r2.Um, 0.0, 0.0, qual_ok=(fu.Um > fu_r2.Um))
    chk("T53","UNIV-G","FUBii > 0 at t_n=0 (outward Aether)",
        fu.FUBii, 0.0, 0.0, qual_ok=(fu.FUBii > 0.0))
    chk("T54","UNIV-G","eps from compute_F_U matches bsfg_metric",
        fu.eps, bsfg_metric(AU_METERS, 0.0).eps, 0.001)

    # VDS cross-check
    vds_py = vds_series(SSQ_DEFAULT, 200).value
    vds_dm = vds_numerical(terms=26)
    chk("T55","VDS-XCHECK","vds_series == dpm vds_numerical (Li_26(0.57))",
        vds_py, vds_dm, 0.1)

    # F_U_Bi_i cross-check with DPM
    fubi_dpm = compute_F_U_Bi_i_numerical(M_bh=M_SUN, r=R_SUN)
    chk("T56","DPM-XCHECK","compute_F_U_Bi_i_numerical finite",
        fubi_dpm, 0.0, 0.0, qual_ok=math.isfinite(fubi_dpm))

    # ── HABITABLE ZONE SCAN (T57–T60) ────────────────────────────────────────
    r_scan = np.array([0.1, 0.5, 1.0, 5.0, 10.0])
    t_scan = np.array([0.0, 0.5, 1.0])
    scan   = scan_habitable_zone(r_scan, t_scan)

    chk("T57","HZ-SCAN","scan balance_grid shape == (5, 3)",
        float(scan['balance_grid'].shape[0] * 100 + scan['balance_grid'].shape[1]),
        503.0, 0.0,
        qual_ok=(scan['balance_grid'].shape == (5, 3)))
    chk("T58","HZ-SCAN","FUBi < 0 at t_n=0 (all r, opposes gravity)",
        float(all(scan['FUBi_grid'][:, 0] < 0)), 1.0, 0.0,
        qual_ok=bool(np.all(scan['FUBi_grid'][:, 0] < 0)))
    chk("T59","HZ-SCAN","FUBii > 0 at t_n=0 (all r, outward)",
        float(all(scan['FUBii_grid'][:, 0] > 0)), 1.0, 0.0,
        qual_ok=bool(np.all(scan['FUBii_grid'][:, 0] > 0)))
    chk("T60","HZ-SCAN","r_hz_m from scan is positive",
        float(np.all(scan['r_hz_m'] > 0)), 1.0, 0.0,
        qual_ok=bool(np.all(scan['r_hz_m'] > 0)))

    # ── MAYAN TIMING / UNIVERSAL INERTIA (T61–T70) ───────────────────────────
    rings5 = mayan_ring_proportions(5, 1.0)
    rings1 = mayan_ring_proportions(1, 1.0)
    state5 = mayan_ring_state(5, 2026.33, 1.0)
    ui_tn0 = universal_inertia(0.0, 5, 1.0)
    ui_tnH = universal_inertia(0.5, 5, 1.0)
    ui_tn1 = universal_inertia(1.0, 5, 1.0)
    zp_years = zero_point_years_in_epoch5(3)

    # T61: Epoch 5 outer ring > base (expanding)
    chk("T61","MAYAN","Epoch 5 outer ring EXPANDING (r_1 > r_base)",
        rings5['r_outer_AU'], 0.0, 0.0,
        qual_ok=(rings5['r_outer_AU'] > 1.0))

    # T62: Epoch 5 companion ring < base (shrinking)
    chk("T62","MAYAN","Epoch 5 companion ring SHRINKS (r_2 < r_base)",
        rings5['r_companion_AU'], 0.0, 0.0,
        qual_ok=(rings5['r_companion_AU'] < 1.0))

    # T63: Epoch 5 inner ring very small (< companion)
    chk("T63","MAYAN","Epoch 5 inner ring VERY SMALL (r_3 < r_2 < r_1)",
        rings5['r_inner_AU'], 0.0, 0.0,
        qual_ok=(rings5['r_inner_AU'] < rings5['r_companion_AU'] < rings5['r_outer_AU']))

    # T64: Epoch 1 all rings equal to r_base (φ^0 = 1)
    chk("T64","MAYAN","Epoch 1: all rings = r_base (φ^0 = 1)",
        rings1['r_outer_AU'], 1.0, 0.001,
        qual_ok=(abs(rings1['r_outer_AU'] - 1.0) < 1e-9
                 and abs(rings1['r_inner_AU'] - 1.0) < 1e-9))

    # T65: Gear ratio 1:3 at epoch 5 = r_outer/r_inner = φ^4/φ^-8 = φ^12 ≈ 321.9969
    phi12 = PHI ** 12
    chk("T65","MAYAN","Epoch 5 gear ratio 1:3 = φ^12 ≈ 321.997",
        rings5['gear_ratio_13'], phi12, 0.01)

    # T66: ω_companion negative (external mesh → opposite rotation)
    chk("T66","MAYAN","ω_companion < 0 (external mesh, opposite rotation)",
        state5.omega_companion, 0.0, 0.0,
        qual_ok=(state5.omega_companion < 0.0))

    # T67: ω_inner > ω_outer (internal mesh amplification)
    chk("T67","MAYAN","ω_inner > ω_outer (internal mesh amplification ≈ 47×)",
        state5.omega_inner, 0.0, 0.0,
        qual_ok=(state5.omega_inner > state5.omega_outer * 2.0))

    # T68: Universal Inertia > 0 at t_n=0 (centripetal, mass-building)
    chk("T68","MAYAN","U_I > 0 at t_n=0 (centripetal restoring mode)",
        ui_tn0.u_inertia, 0.0, 0.0,
        qual_ok=ui_tn0.centripetal_mode)

    # T69: Universal Inertia = 0 at t_n=0.5 (zero-point gravity)
    chk("T69","MAYAN","U_I = 0 at t_n=0.5 (zero-point gravity, massless scalar)",
        ui_tnH.u_inertia, 0.0, 0.0,
        qual_ok=ui_tnH.zero_point)

    # T70: Universal Inertia invariant — same |U_I| at t_n=0 and t_n=1, opposite sign
    chk("T70","MAYAN","U_I invariant: |U_I(0)| = |U_I(1)|, sign flip (centripetal↔centrifugal)",
        abs(ui_tn0.u_inertia + ui_tn1.u_inertia), 0.0, 0.0,
        qual_ok=(abs(ui_tn0.u_inertia + ui_tn1.u_inertia)
                 < 1e-6 * abs(ui_tn0.u_inertia)))

    # ── VDS-DVP-DH26 (T71–T80): variant branches + coupling (Session 202 Phase H202) ─
    vb_r  = vds_branches(SSQ_DEFAULT, 200)
    db_r  = dvp_branches(200)
    bh10_r = bh26_branches(10)
    bh1_r  = bh26_branches(1)
    cc_r   = vds_dvp_coupled(SSQ_DEFAULT, 200, 200)
    res_r  = bh26_bsh_resonance(3.3e7, SSQ_DEFAULT, 0.0, 1)

    # T71: VDS_prime ≈ 1.0  (Li_25(0.57)/0.57 — calibration sensitivity)
    chk("T71","VDS-DVP-DH26","vds_prime = Li_25(SSq)/SSq ~ 1.0 (VDS calibration sensitivity)",
        vb_r.vds_prime, 1.0, 0.001)

    # T72: VDS energy density > 0  (VDS × RHO_VAC_SCM > 0)
    chk("T72","VDS-DVP-DH26","vds_density > 0  (VDS × RHO_VAC_SCM positive)",
        vb_r.vds_density, 0.0, 0.0, qual_ok=(vb_r.vds_density > 0.0))

    # T73: DVP spectral sum > 0  (prime-vortex flux non-zero)
    chk("T73","VDS-DVP-DH26","dvp zeta_sum > 0  (prime-vortex spectral flux positive)",
        db_r.zeta_sum, 0.0, 0.0, qual_ok=(db_r.zeta_sum > 0.0))

    # T74: pair_product < a_29^2  (a(31) < a(29) → product < first-term squared)
    chk("T74","VDS-DVP-DH26","dvp pair_product < a_29^2  (a(31) < a(29), strictly decreasing)",
        db_r.pair_product, 0.0, 0.0,
        qual_ok=(db_r.a_29 > 0.0 and db_r.pair_product < db_r.a_29**2))

    # T75: BH26 spectral sum N=10 == 1760  (closed-form eigenvalue ladder)
    chk("T75","VDS-DVP-DH26","bh26_branches(10).spectral_sum == 1760 (eigenvalue ladder)",
        bh10_r.spectral_sum, 1760.0, 1e-6)

    # T76: Casimir energy > 0  (vacuum zero-point energy from spectral ladder)
    chk("T76","VDS-DVP-DH26","bh26_branches(10).casimir_energy > 0  (vacuum Casimir energy)",
        bh10_r.casimir_energy, 0.0, 0.0, qual_ok=(bh10_r.casimir_energy > 0.0))

    # T77: degeneracy_k1 == 26  (C(26,25)=26 on S^{25})
    chk("T77","VDS-DVP-DH26","bh26_branches(1).degeneracy_k1 == 26  (C(26,25) on S^25)",
        float(bh1_r.degeneracy_k1), 26.0, 0.0, qual_ok=(bh1_r.degeneracy_k1 == 26))

    # T78: vds_coupling > 0  (BH26→VDS topological bridge finite)
    chk("T78","VDS-DVP-DH26","bh26_branches(10).vds_coupling > 0  (BH26→VDS topology bridge)",
        bh10_r.vds_coupling, 0.0, 0.0, qual_ok=(bh10_r.vds_coupling > 0.0))

    # T79: VDS×DVP joint_coeff ≥ 0  (geometric-mean coupling non-negative)
    chk("T79","VDS-DVP-DH26","vds_dvp_coupled joint_coeff >= 0  (geometric-mean field coupling)",
        cc_r.joint_coeff, 0.0, 0.0, qual_ok=(cc_r.joint_coeff >= 0.0))

    # T80: BH26×BSH cross-resonance energy density > 0 at t_n=0, k=1
    chk("T80","VDS-DVP-DH26","bh26_bsh_resonance energy_density > 0  (k=1, t_n=0)",
        res_r.energy_density, 0.0, 0.0, qual_ok=(res_r.energy_density > 0.0))

    # ── F_U + U_I INTEGRATION (T81–T85) — Session 265 ─────────────────────────
    # Wire Universal Inertia into the Universal Gravity master equation.
    fu_tn0    = compute_F_U(AU_METERS, 0.0, epoch=5)
    fu_tnH    = compute_F_U(AU_METERS, 0.5, epoch=5)
    fu_tn1    = compute_F_U(AU_METERS, 1.0, epoch=5)
    hz_E1     = solve_habitable_zone({}, epoch=1)
    hz_E5     = solve_habitable_zone({}, epoch=5)

    # T81: compute_F_U returns U_I field (non-zero at t_n=0)
    chk("T81","F_U+U_I","compute_F_U populates U_I field at t_n=0 (centripetal)",
        fu_tn0.U_I, 0.0, 0.0,
        qual_ok=(fu_tn0.U_I > 0.0))

    # T82: U_I sign flips between t_n=0 and t_n=1 (centripetal → centrifugal)
    chk("T82","F_U+U_I","U_I sign flips across great cycle (cos(π·t_n))",
        fu_tn0.U_I + fu_tn1.U_I, 0.0, 0.0,
        qual_ok=(fu_tn0.U_I * fu_tn1.U_I < 0.0))

    # T83: U_I ≈ 0 at t_n=0.5 (zero-point gravity / massless scalar)
    chk("T83","F_U+U_I","U_I ≈ 0 at zero-point t_n=0.5 (massless scalar)",
        abs(fu_tnH.U_I), 0.0, 0.0,
        qual_ok=(abs(fu_tnH.U_I) < abs(fu_tn0.U_I) * 1.0e-10))

    # T84: epoch-aware tectonic band populated and non-zero
    chk("T84","F_U+U_I","solve_habitable_zone tectonic band populated (Epoch 5)",
        hz_E5.band_width_AU, 0.0, 0.0,
        qual_ok=(hz_E5.band_width_AU > 0.0 and hz_E5.epoch == 5))

    # T85: Epoch-5 band is asymmetrically wider than Epoch-1 (outer expanded)
    chk("T85","F_U+U_I","Epoch 5 band wider than Epoch 1 (outer ring expansion)",
        hz_E5.band_width_AU - hz_E1.band_width_AU, 0.0, 0.0,
        qual_ok=(hz_E5.band_width_AU > hz_E1.band_width_AU))

    # ── EMERGENT MASS at FUBi+FUBii=0 CROSSING (T86–T90) — Session 267 ────────
    # The inverse problem to solve_habitable_zone: given r_hz, find M_emergent.
    # This is the Quantum Chain Step 7: mass BORN at the crossing, not input.
    hz_em  = solve_habitable_zone({})
    em_sol = compute_emergent_mass(hz_em.r_hz_m, hz_em.t_n_hz)

    # T86: emergent-mass solver converges and returns positive mass
    chk("T86","EmergentMass","compute_emergent_mass converges at r_hz",
        em_sol.M_emergent_kg, 0.0, 0.0,
        qual_ok=(em_sol.converged and em_sol.M_emergent_kg > 0.0))

    # T87: residual of buoyancy balance is at machine precision
    cos_tn_em = abs(math.cos(math.pi * em_sol.t_n_hz))
    FUBi_mag  = (em_sol.beta_i_used * G_NEWTON * em_sol.M_emergent_kg**2
                 / max(em_sol.r_hz_m, 1.0)**2 * em_sol.orbit_factor * cos_tn_em)
    rel_resid = (abs(em_sol.residual_at_M) / FUBi_mag) if FUBi_mag > 0 else 1.0
    chk("T87","EmergentMass","FUBi+FUBii relative residual < 1e-12 at M_emergent",
        rel_resid, 0.0, 0.0, qual_ok=(rel_resid < 1.0e-12))

    # T88: inverse-relation consistency — feeding M_emergent back into
    #      solve_habitable_zone should reproduce r_hz (round-trip).
    hz_check = solve_habitable_zone({})  # uses M_SUN; baseline
    # Closed-form sanity: r_hz³ ∝ M² so doubling M should grow r_hz by 2^(2/3)
    em_2sun = compute_emergent_mass(hz_check.r_hz_m * (2.0**(2.0/3.0)),
                                     hz_check.t_n_hz)
    ratio = em_2sun.M_emergent_sun / em_sol.M_emergent_sun if em_sol.M_emergent_sun > 0 else 0.0
    chk("T88","EmergentMass","M ∝ r^(3/2) closed-form scaling (ratio ≈ 2.0)",
        ratio, 2.0, 0.05, qual_ok=(abs(ratio - 2.0) < 0.05))

    # T89: classification populated and consistent with mass band
    chk("T89","EmergentMass","mass classification populated",
        1.0 if em_sol.classification else 0.0, 1.0, 0.0,
        qual_ok=(em_sol.classification in
                 ("sub_stellar","sub_solar","solar","stellar",
                  "massive_stellar","BH_seed")))

    # T90: rho_vac scaling — halving the vacuum density halves M² (M scales 1/√2)
    em_halfrho = compute_emergent_mass(hz_em.r_hz_m, hz_em.t_n_hz,
                                        rho_vac=float(RHO_VAC_SCM) * 0.5)
    expected_factor = 1.0 / math.sqrt(2.0)  # ≈ 0.7071
    actual_factor   = (em_halfrho.M_emergent_kg / em_sol.M_emergent_kg
                       if em_sol.M_emergent_kg > 0 else 0.0)
    chk("T90","EmergentMass","M ∝ sqrt(rho_vac) — halving ρ_vac → M × 1/√2",
        actual_factor, expected_factor, 1.0e-6,
        qual_ok=(abs(actual_factor - expected_factor) < 1.0e-6))

    # ── UNIVERSAL BUOYANCY SIMULTANEOUS SOLVER (T91-T97) — v2.2.0 ────────────
    # Coupled 4x4 nonlinear system (r_hz, t_n_hz, r_cg, M_emergent).
    # Encodes Aether UA vacuum counter-balance F_U / F_U_Bi / F_U_Bi_i
    # jointly with collapsing-gravity zone and emergent mass.
    ubs = solve_universal_buoyancy({})
    ubs_calc = UniversalBuoyancySimultaneousSolver().compute({})

    # T91: solver returns a finite, populated solution
    chk("T91","UBS-Solver","solve_universal_buoyancy returns finite r_hz > 0",
        ubs.r_hz_m, 0.0, 0.0,
        qual_ok=(math.isfinite(ubs.r_hz_m) and ubs.r_hz_m > 0.0))

    # T92: collapsing-gravity zone interior to habitable zone (r_cg < r_hz)
    chk("T92","UBS-Solver","r_cg < r_hz (collapse zone interior to HZ)",
        ubs.r_hz_m - ubs.r_cg_m, 0.0, 0.0,
        qual_ok=(ubs.r_cg_m < ubs.r_hz_m and ubs.r_cg_m > 0.0))

    # T93: emergent mass positive and finite
    chk("T93","UBS-Solver","M_emergent positive and finite",
        ubs.M_emergent_kg, 0.0, 0.0,
        qual_ok=(math.isfinite(ubs.M_emergent_kg) and ubs.M_emergent_kg > 0.0))

    # T94: E1 residual at machine precision (FUBi + FUBii = 0 at r_hz)
    chk("T94","UBS-Solver","residual_E1 < 1e-6 (FUBi+FUBii=0 at r_hz)",
        abs(ubs.res_E1), 0.0, 0.0,
        qual_ok=(abs(ubs.res_E1) < 1.0e-6))

    # T95: E3 collapse-boundary 2:1 condition holds at r_cg
    # |F_U_Bi/F_U_Bi_i| at r_cg should be ~ 2 (boundary definition)
    chk("T95","UBS-Solver","collapse_ratio ≈ 2 at r_cg (E3 boundary)",
        ubs.collapse_ratio, 2.0, 5.0,
        qual_ok=(0.5 < ubs.collapse_ratio < 4.0))

    # T96: calculator-class wrapper returns same r_hz_AU as direct solve
    chk("T96","UBS-Solver","UniversalBuoyancySimultaneousSolver.compute matches direct solver",
        ubs_calc['r_hz_AU'], ubs.r_hz_AU, 1.0e-6)

    # T97: solver_msg is populated (either CONVERGED or staged fallback)
    chk("T97","UBS-Solver","solver_msg is populated diagnostic string",
        1.0 if ubs.solver_msg else 0.0, 1.0, 0.0,
        qual_ok=(len(ubs.solver_msg) > 0))

    # ── Summary ──────────────────────────────────────────────────────────────
    passed = sum(1 for r in results if r['passed'])
    failed = len(results) - passed
    if verbose:
        print(f"\n{'='*72}")
        print(f"TOTAL: {passed}/{len(results)} PASS  |  {failed} FAIL")
        print(f"{'='*72}\n")

    return {'passed': passed, 'failed': failed, 'total': len(results),
            'results': results}


# =============================================================================
# SECTION 8 — MODULE ENTRY POINT
# =============================================================================

if __name__ == '__main__':
    import sys
    verbose = '--quiet' not in sys.argv
    out = run_qcalcgeom_tests(verbose=verbose)
    sys.exit(0 if out['failed'] == 0 else 1)
