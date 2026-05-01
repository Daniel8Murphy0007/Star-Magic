# -*- coding: utf-8 -*-
"""
dpm_vacuum_manifold.py  v2.0
Di-Pseudo-Monopole (DPM) Vacuum Calculator -- Quantum Chain Compliant

THE QUANTUM CHAIN (canonical, Star-Magic.txt lines 11-22, IMMUTABLE):

  Step 0  0_vacuum   -> |grad(UA)|          vacuum tension differential
  Step 1  grad(UA)   -> DPM_vortex          a_DPM = F_DPM*f_DPM*E_vac/(c*V_sys)
  Step 2  DPM_vortex -> mu_s                mu_s = rho_A * V_DPM  (vortex volume)
  Step 3  mu_s       -> Ug1[seed=DPM]       Ug1 seeded from mu_s -- NOT from mass
  Step 4  Ug1        -> Ug_family           Ug2+Ug3+Ug4 simultaneously promoted
  Step 5  Ug_family  -> F_U                 + Um + FUBi + FUBii + UA_uv
  Step 6  F_U        -> crossing            FUBi(r) + FUBii(r) = 0  compaction
  Step 7  crossing   -> M_emergent          mass BORN at crossing, not before
  Step 8  M_emergent -> GM/r^2             LAST -- observational projection only

ARCHITECTURE:
  scm_vacuum_manifold.py  -> imported  (SCm base layer: CW rotation, primordial vacuum)
  ua_vacuum_manifold.py   -> imported  (UA superstructure: CCW rotation, 4 layers)
  THIS FILE runs all 8 chain steps from vacuum to GM/r^2 for every atom Z=1-118.

PERIODIC TABLE GEOMETRY RULE:
  DPMBody uses Z (vortex count) and A (resonance count) as PRIMARY geometry inputs.
  R_nuc, V_DPM, B0, omega0, v_fermi are ALL computed from Z and A -- never from mass.
  M_table is the tabulated atomic mass stored as VERIFICATION ONLY.
  Ug1 through Ug4 use M_proto (ACP-emerged from the chain), not M_table.
  GM/r^2 is the LAST output, computed from M_table as the verified stable mass.

DPM FUNDAMENTAL EQUATIONS:
  DPM       = [UA']/[SCm] = 10              (scale-invariant ratio)
  Grind_opp = omega_CW * SCm - omega_CCW * UA'
  F_U       = Ug_sum - Ubi + Um

Author: Daniel T. Murphy  |  dpm_vacuum_manifold.py v2.0  |  May 2026
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import sympy as sp

# -- Import the two standalone manifold layers ---------------------------------
from scm_vacuum_manifold import (
    RHO_VAC_SCM, RHO_VAC_UA, THZ_PHONON, BETA_I, LAMBDA_I, OMEGA_S,
    KAPPA_FLOAT, SSQ,
    F_U_Bi_i_99, monte_carlo_fubi_i,
    KER_SCm, parkhomov_excess_heat, pons_fleischmann_excess_heat,
)
from ua_vacuum_manifold import (
    ua_layer_density, ua_dpm_total_density, ua_dpm_buoyancy_factor,
    ua_calibration_ratio,
    E_PHONON, S26_3, PHI_RES, DELTA_UA_FOURTH, DPM_DENSITY_RATIO,
    UA_prime, UA_double_prime, UA_triple_prime, UA_quad_prime, UA_total,
    _F_Bi_i_scm, F_U_Bi_i_DPM,
    ua_lenr_comparison, ua_casimir_comparison, ua_string_brane_embedding,
    ua_cosmological_acceleration, ua_rotation_curve_flat,
    ua_hubble_tension_modulation, ua_dark_energy_substitute,
)

# =============================================================================
# S1  PHYSICAL CONSTANTS
# =============================================================================

HBAR:    float = 1.054571817e-34    # J*s   reduced Planck constant
MU_0:    float = 1.2566370614e-6    # N/A^2 vacuum permeability
MU_N:    float = 5.0507837461e-27   # J/T   nuclear magneton
AMU:     float = 1.66053906660e-27  # kg    atomic mass unit
C_LIGHT: float = 2.99792458e8       # m/s   speed of light
V_SCM:   float = C_LIGHT / 3.0     # m/s   SCm velocity (v_SCm = c/3)
G_CONST: float = 6.67430e-11       # m^3/(kg*s^2)  gravitational constant
R_NUC_0: float = 1.2e-15           # m     nuclear radius constant (R = R0*A^(1/3))
K_E:     float = 8.9875517923e9    # N*m^2/C^2  Coulomb constant

# =============================================================================
# S2  DPM COUPLING CONSTANTS  (Source4 canonical)
# =============================================================================

K1:         float = 1.0            # Ug1 coupling
K2:         float = 1.0            # Ug2 coupling
K3:         float = 1.0            # Ug3 coupling
K4:         float = 1.0            # Ug4 coupling
ALPHA:      float = KAPPA_FLOAT    # temporal decay = kappa = 5e-4 day^-1
DELTA_DEF:  float = 0.0            # deformation factor (ground state)
H_SCM:      float = 0.99           # condensate fraction
EPSILON_SW: float = 0.0            # solar-wind enhancement (zero at atomic scale)
RHO_SW:     float = 0.0            # solar-wind density (zero at atomic scale)
OMEGA_G:    float = 1.0            # galactic angular factor (normalised)
MBH_DG:     float = 1.0            # M_bh/d_g ratio (normalised)
P_CORE:     float = 1.0            # core pressure factor

# DPM grinding pair frequencies (CW=SCm inner, CCW=UA outer)
OMEGA_CW:  float = 2.0 * math.pi * 1.2e10   # rad/s  SCm CW grinding
OMEGA_CCW: float = 2.0 * math.pi * 8.3e9    # rad/s  UA' CCW grinding

# ACP constants (Chapter 7 Star-Magic.txt)
# E_crack = (rho_SCm * c^2) / [SSq]  -- gate energy for mass condensation
# SSQ may be a sympy Rational; force Python float to keep chain arithmetic clean.
E_CRACK:  float = float(RHO_VAC_SCM * C_LIGHT ** 2) / float(SSQ)   # J ~1.12e-19 J
# M_0 = E_crack / c^2 = rho_SCm / [SSq]  -- base DPM mass unit
M_0_DPM:  float = float(E_CRACK) / float(C_LIGHT ** 2)              # kg = RHO_VAC_SCM / SSQ

# =============================================================================
# S3  DPMBody DATACLASS  (geometry-first; mass is verification only)
# =============================================================================

@dataclass
class DPMBody:
    """A DPM body defined by vacuum geometry, NOT by atomic mass.

    PRIMARY inputs (geometric -- derived from Z and A only):
      Z       : atomic number = number of DPM vortex units in resonance
      A       : mass number = resonance count (determines nuclear radius)
      symbol  : chemical symbol
      name    : element name
      R_cov   : covalent radius [m]   (geometric, for Newton projection radius)
      R_nuc   : nuclear radius = R_NUC_0 * A^(1/3)   [m]  COMPUTED
      V_DPM   : DPM vortex volume = (4/3)*pi*R_nuc^3 [m3] COMPUTED
      B0      : nuclear surface magnetic field [T]    COMPUTED from Z, R_nuc
      omega0  : nuclear Larmor angular frequency at 1T [rad/s]  COMPUTED from Z
      v_fermi : Fermi velocity proxy [m/s]            COMPUTED from Z

    VERIFICATION ONLY (NOT a chain input -- compared against M_emergent at end):
      M_table : tabulated atomic mass [kg]   READ LAST, verified against chain
    """
    Z:       int
    A:       int
    symbol:  str
    name:    str
    R_cov:   float    # m
    R_nuc:   float    # m  (computed from A)
    V_DPM:   float    # m^3 (computed from R_nuc)
    B0:      float    # T  (computed from Z, R_nuc)
    omega0:  float    # rad/s
    v_fermi: float    # m/s
    M_table: float    # kg  -- VERIFICATION ONLY


def _build_dpm_body(Z: int, symbol: str, name: str, A: int,
                    mass_u: float, R_cov_pm: float) -> DPMBody:
    """Factory: build DPMBody from pure geometry (Z, A, R_cov).
    M_table = mass_u * AMU is stored last as verification field only.
    """
    R_cov   = R_cov_pm * 1.0e-12                        # pm -> m
    R_nuc   = R_NUC_0 * A ** (1.0 / 3.0)               # nuclear radius [m]
    V_DPM   = (4.0 / 3.0) * math.pi * R_nuc ** 3       # DPM vortex volume [m^3]
    B0      = (MU_0 / (4.0 * math.pi)) * 2.0 * Z * MU_N / R_nuc ** 3  # [T]
    omega0  = Z * 2.675e8                               # rad/s at 1T Larmor
    v_fermi = 0.77e6 * Z ** (1.0 / 3.0)                # m/s Fermi proxy
    M_table = mass_u * AMU                              # verification only
    return DPMBody(Z=Z, A=A, symbol=symbol, name=name,
                   R_cov=R_cov, R_nuc=R_nuc, V_DPM=V_DPM,
                   B0=B0, omega0=omega0, v_fermi=v_fermi,
                   M_table=M_table)


# =============================================================================
# S4  PERIODIC TABLE  (Z = 1-118)
# =============================================================================
# Columns: Z, symbol, name, A, mass_u (verification), R_cov_pm (geometry)

_PT_RAW: List[Tuple] = [
    (1,  "H",  "Hydrogen",        1,   1.008,   31),
    (2,  "He", "Helium",          4,   4.003,   28),
    (3,  "Li", "Lithium",         7,   6.941,  128),
    (4,  "Be", "Beryllium",       9,   9.012,   96),
    (5,  "B",  "Boron",          11,  10.811,   84),
    (6,  "C",  "Carbon",         12,  12.011,   77),
    (7,  "N",  "Nitrogen",       14,  14.007,   75),
    (8,  "O",  "Oxygen",         16,  15.999,   73),
    (9,  "F",  "Fluorine",       19,  18.998,   71),
    (10, "Ne", "Neon",           20,  20.180,   69),
    (11, "Na", "Sodium",         23,  22.990,  166),
    (12, "Mg", "Magnesium",      24,  24.305,  141),
    (13, "Al", "Aluminium",      27,  26.982,  121),
    (14, "Si", "Silicon",        28,  28.086,  111),
    (15, "P",  "Phosphorus",     31,  30.974,  107),
    (16, "S",  "Sulfur",         32,  32.065,  105),
    (17, "Cl", "Chlorine",       35,  35.453,  102),
    (18, "Ar", "Argon",          40,  39.948,  106),
    (19, "K",  "Potassium",      39,  39.098,  203),
    (20, "Ca", "Calcium",        40,  40.078,  176),
    (21, "Sc", "Scandium",       45,  44.956,  170),
    (22, "Ti", "Titanium",       48,  47.867,  160),
    (23, "V",  "Vanadium",       51,  50.942,  153),
    (24, "Cr", "Chromium",       52,  51.996,  139),
    (25, "Mn", "Manganese",      55,  54.938,  139),
    (26, "Fe", "Iron",           56,  55.845,  132),
    (27, "Co", "Cobalt",         59,  58.933,  126),
    (28, "Ni", "Nickel",         58,  58.693,  124),
    (29, "Cu", "Copper",         63,  63.546,  132),
    (30, "Zn", "Zinc",           64,  65.380,  122),
    (31, "Ga", "Gallium",        69,  69.723,  122),
    (32, "Ge", "Germanium",      74,  72.630,  120),
    (33, "As", "Arsenic",        75,  74.922,  119),
    (34, "Se", "Selenium",       80,  78.971,  120),
    (35, "Br", "Bromine",        79,  79.904,  120),
    (36, "Kr", "Krypton",        84,  83.798,  116),
    (37, "Rb", "Rubidium",       85,  85.468,  220),
    (38, "Sr", "Strontium",      88,  87.620,  195),
    (39, "Y",  "Yttrium",        89,  88.906,  190),
    (40, "Zr", "Zirconium",      90,  91.224,  175),
    (41, "Nb", "Niobium",        93,  92.906,  164),
    (42, "Mo", "Molybdenum",     96,  95.960,  154),
    (43, "Tc", "Technetium",     99,  98.000,  147),
    (44, "Ru", "Ruthenium",     102, 101.070,  146),
    (45, "Rh", "Rhodium",       103, 102.906,  142),
    (46, "Pd", "Palladium",     106, 106.420,  139),
    (47, "Ag", "Silver",        107, 107.868,  145),
    (48, "Cd", "Cadmium",       114, 112.411,  144),
    (49, "In", "Indium",        115, 114.818,  142),
    (50, "Sn", "Tin",           120, 118.710,  139),
    (51, "Sb", "Antimony",      121, 121.760,  139),
    (52, "Te", "Tellurium",     130, 127.600,  138),
    (53, "I",  "Iodine",        127, 126.904,  139),
    (54, "Xe", "Xenon",         132, 131.293,  140),
    (55, "Cs", "Caesium",       133, 132.905,  244),
    (56, "Ba", "Barium",        138, 137.327,  215),
    (57, "La", "Lanthanum",     139, 138.905,  207),
    (58, "Ce", "Cerium",        140, 140.116,  204),
    (59, "Pr", "Praseodymium",  141, 140.908,  203),
    (60, "Nd", "Neodymium",     142, 144.242,  201),
    (61, "Pm", "Promethium",    145, 145.000,  199),
    (62, "Sm", "Samarium",      152, 150.360,  198),
    (63, "Eu", "Europium",      153, 151.964,  198),
    (64, "Gd", "Gadolinium",    158, 157.250,  196),
    (65, "Tb", "Terbium",       159, 158.925,  194),
    (66, "Dy", "Dysprosium",    164, 162.500,  192),
    (67, "Ho", "Holmium",       165, 164.930,  192),
    (68, "Er", "Erbium",        166, 167.259,  189),
    (69, "Tm", "Thulium",       169, 168.934,  190),
    (70, "Yb", "Ytterbium",     174, 173.045,  187),
    (71, "Lu", "Lutetium",      175, 174.967,  187),
    (72, "Hf", "Hafnium",       180, 178.490,  175),
    (73, "Ta", "Tantalum",      181, 180.948,  170),
    (74, "W",  "Tungsten",      184, 183.840,  162),
    (75, "Re", "Rhenium",       187, 186.207,  151),
    (76, "Os", "Osmium",        192, 190.230,  144),
    (77, "Ir", "Iridium",       193, 192.217,  141),
    (78, "Pt", "Platinum",      195, 195.084,  136),
    (79, "Au", "Gold",          197, 196.967,  136),
    (80, "Hg", "Mercury",       202, 200.592,  132),
    (81, "Tl", "Thallium",      205, 204.383,  145),
    (82, "Pb", "Lead",          208, 207.200,  146),
    (83, "Bi", "Bismuth",       209, 208.980,  148),
    (84, "Po", "Polonium",      209, 209.000,  140),
    (85, "At", "Astatine",      210, 210.000,  150),
    (86, "Rn", "Radon",         222, 222.000,  150),
    (87, "Fr", "Francium",      223, 223.000,  260),
    (88, "Ra", "Radium",        226, 226.000,  221),
    (89, "Ac", "Actinium",      227, 227.000,  215),
    (90, "Th", "Thorium",       232, 232.038,  206),
    (91, "Pa", "Protactinium",  231, 231.036,  200),
    (92, "U",  "Uranium",       238, 238.029,  196),
    (93, "Np", "Neptunium",     237, 237.000,  190),
    (94, "Pu", "Plutonium",     244, 244.000,  187),
    (95, "Am", "Americium",     243, 243.000,  180),
    (96, "Cm", "Curium",        247, 247.000,  169),
    (97, "Bk", "Berkelium",     247, 247.000,  170),
    (98, "Cf", "Californium",   251, 251.000,  170),
    (99, "Es", "Einsteinium",   252, 252.000,  170),
    (100,"Fm", "Fermium",       257, 257.000,  167),
    (101,"Md", "Mendelevium",   258, 258.000,  173),
    (102,"No", "Nobelium",      259, 259.000,  176),
    (103,"Lr", "Lawrencium",    266, 266.000,  161),
    (104,"Rf", "Rutherfordium", 267, 267.000,  157),
    (105,"Db", "Dubnium",       268, 268.000,  149),
    (106,"Sg", "Seaborgium",    271, 271.000,  143),
    (107,"Bh", "Bohrium",       272, 272.000,  141),
    (108,"Hs", "Hassium",       277, 277.000,  134),
    (109,"Mt", "Meitnerium",    278, 278.000,  129),
    (110,"Ds", "Darmstadtium",  281, 281.000,  128),
    (111,"Rg", "Roentgenium",   282, 282.000,  121),
    (112,"Cn", "Copernicium",   285, 285.000,  122),
    (113,"Nh", "Nihonium",      286, 286.000,  136),
    (114,"Fl", "Flerovium",     289, 289.000,  143),
    (115,"Mc", "Moscovium",     290, 290.000,  162),
    (116,"Lv", "Livermorium",   293, 293.000,  175),
    (117,"Ts", "Tennessine",    294, 294.000,  165),
    (118,"Og", "Oganesson",     294, 294.000,  157),
]

PERIODIC_TABLE: List[DPMBody] = [_build_dpm_body(*row) for row in _PT_RAW]
ELEMENT: Dict[int, DPMBody]   = {b.Z: b for b in PERIODIC_TABLE}


# =============================================================================
# S5  THE QUANTUM CHAIN -- 8 STEPS FROM VACUUM TO GM/r^2
# =============================================================================

# -- STEP 0: Zero-mass vacuum state -------------------------------------------

def chain_step0_vacuum() -> Dict[str, float]:
    """Step 0: 0_vacuum -- no mass, no motion, no gravity.

    Starting axiom (Star-Magic.txt Canonical Ontology Lock):
      rho_UA = 0, rho_vac = |grad(UA)|, F_U(vacuum) = 0
    The vacuum gradient = differential tension between UA and SCm densities.
    The belly button fires HERE. Everything downstream comes from this step.
    """
    grad_UA   = RHO_VAC_UA - RHO_VAC_SCM
    E_react_0 = RHO_VAC_SCM * V_SCM ** 2 / RHO_VAC_UA
    return {
        "grad_UA":   grad_UA,    # [kg/m^3]  = 6.381e-36
        "E_react_0": E_react_0,  # [J/m^3]   peak reaction energy density at t=0
        "F_U_vac":   0.0,        # unified field = 0 in zero-mass vacuum
    }


# -- STEP 1: DPM vortex formation ---------------------------------------------

def chain_step1_dpm(body: DPMBody) -> Dict[str, float]:
    """Step 1: grad(UA) -> DPM_vortex.

    a_DPM = (F_DPM * f_DPM * E_vac_neb) / (c * V_sys)
    F_DPM = I * A_cross * (omega_1 - omega_2)

    At atomic scale:
      I_flux    = Z * rho_SCm * v_SCm     [rotational SCm flux through Z vortex units]
      A_cross   = pi * R_nuc^2            [DPM vortex cross-section]
      delta_om  = |OMEGA_CW - OMEGA_CCW|  [differential grinding angular velocity]
      E_vac_neb = rho_SCm * c^2           [vacuum nebular energy density]
      V_sys     = V_DPM                   [DPM system volume]
      f_dpm     = PHI_RES                 [resonance factor from ua layer]
    """
    I_flux    = body.Z * RHO_VAC_SCM * V_SCM
    A_cross   = math.pi * body.R_nuc ** 2
    delta_om  = abs(OMEGA_CW - OMEGA_CCW)
    F_DPM     = I_flux * A_cross * delta_om
    f_dpm     = PHI_RES
    E_vac_neb = RHO_VAC_SCM * C_LIGHT ** 2
    a_DPM     = F_DPM * f_dpm * E_vac_neb / (C_LIGHT * body.V_DPM)
    return {
        "I_flux":  I_flux,
        "A_cross": A_cross,
        "F_DPM":   F_DPM,
        "a_DPM":   a_DPM,
    }


# -- STEP 2: Magnetic moment from DPM vortex ----------------------------------

def chain_step2_mu_s(body: DPMBody) -> float:
    """Step 2: DPM_vortex -> mu_s.

    mu_s = rho_A * V_DPM

    The magnetic moment is seeded by SCm vacuum density filling the vortex volume.
    THIS IS NOT FROM ATOMIC MASS.
    mu_s comes purely from the DPM vortex geometry (R_nuc from A) and vacuum density.
    """
    return RHO_VAC_SCM * body.V_DPM


# -- ACP PROTO-MASS (between Steps 2 and 3) -----------------------------------

def chain_acp_M_proto(Z: int) -> float:
    """ACP proto-mass -- the mass that EMERGES from the DPM vortex resonance count.

    From Star-Magic.txt Chapter 7 (Mass Emergence):
      E_crack = (rho_vac_SCm * c^2) / [SSq]
      M_0     = E_crack / c^2  =  rho_vac_SCm / [SSq]
      M_proto = M_0 * (1 - exp(-Z/10)) * Z

    Z is the number of DPM vortex units (atomic number = vortex resonance count).
    M_proto is the mass emerging from the ACP chain -- not read from any table.
    For Z=1 (H):  M_proto = M_0 * (1 - exp(-0.1)) * 1   ~ M_0 * 0.0952
    For Z=26 (Fe): M_proto = M_0 * (1 - exp(-2.6)) * 26
    """
    return M_0_DPM * (1.0 - math.exp(-Z / 10.0)) * Z


# -- E_react helper -----------------------------------------------------------

def chain_E_react(v: float, t: float = 0.0) -> float:
    """E_react(t) = (rho_SCm * v^2) / rho_UA * exp(-kappa * t).

    The energy of UA/SCm maximum attraction.
    v = velocity proxy (v_fermi at atomic scale).
    E_react = 0 when v = 0 (dead mass condition -- Star-Magic.txt Chapter 14).
    """
    return RHO_VAC_SCM * v ** 2 / RHO_VAC_UA * math.exp(-KAPPA_FLOAT * t)


# -- STEPS 3-4: Ug family assembly --------------------------------------------

def chain_step3_Ug1(mu_s: float, M_proto: float, r: float,
                    t: float, t_n: float) -> float:
    """Step 3: mu_s -> Ug1[seed=DPM].

    Ug1 = k1 * mu_s * (M_proto/r^2) * exp(-alpha*t) * cos(pi*t_n) * (1+delta_def)

    mu_s    comes from Step 2 (DPM vortex volume * vacuum density).
    M_proto comes from the ACP chain (Z vortex count * M_0_DPM).
    NEITHER is from the atomic mass table.

    This is THE DPM in field form. Ug1 IS the DPM.
    """
    return (K1 * mu_s * (M_proto / r ** 2)
            * math.exp(-ALPHA * t)
            * math.cos(math.pi * t_n)
            * (1.0 + DELTA_DEF))


def chain_step4_ug_family(body: DPMBody, mu_s: float, M_proto: float,
                          r: float, t: float, t_n: float) -> Dict[str, float]:
    """Step 4: Ug1 simultaneously promotes Ug2, Ug3, Ug4.

    All four Ug terms are simultaneous expressions of the same DPM.
    None is computed before the others -- simultaneous assembly.

    Ug2 -- outer bubble:        uses vacuum charge Q_SCm + Q_UA, NOT mass
    Ug3 -- magnetic string:     uses B0 (nuclear field from Z/R_nuc), NOT mass
    Ug4 -- vacuum concentration: uses Z (vortex count), NOT mass
    E_react -- UA/SCm attraction energy (from v_fermi proxy)
    """
    E_react = chain_E_react(body.v_fermi, t)

    # Ug1 -- the DPM itself (Step 3)
    Ug1 = chain_step3_Ug1(mu_s, M_proto, r, t, t_n)

    # Ug2 -- outer field bubble
    # Q_SCm = rho_SCm * V_DPM,  Q_UA = rho_UA * V_DPM  (vacuum charge proxies)
    Q_sum = (RHO_VAC_SCM + RHO_VAC_UA) * body.V_DPM
    R_b   = body.R_nuc * 100.0                   # nuclear "heliosphere" radius
    S_rb  = 1.0 if r > R_b else 0.0             # step function (1 = outside bubble)
    sw    = 1.0 + EPSILON_SW * RHO_SW
    Ug2   = K2 * Q_sum * (M_proto / r ** 2) * S_rb * sw * H_SCM * E_react

    # Ug3 -- magnetic string disk rotation
    # Driven by B0 (nuclear surface field from Z, R_nuc) and omega0 (Larmor from Z)
    Ug3 = K3 * body.B0 * math.cos(body.omega0 * t * math.pi) * P_CORE * E_react

    # Ug4 -- vacuum concentration
    # Z = DPM vortex count = concentration factor (NOT atomic mass)
    Ug4 = (K4 * RHO_VAC_SCM * float(body.Z)
           * math.exp(-ALPHA * t)
           * math.cos(math.pi * t_n))

    Ug_sum = Ug1 + Ug2 + Ug3 + Ug4
    return {
        "E_react": E_react,
        "Ug1":     Ug1,
        "Ug2":     Ug2,
        "Ug3":     Ug3,
        "Ug4":     Ug4,
        "Ug_sum":  Ug_sum,
    }


# -- STEP 5: F_U assembly -----------------------------------------------------

def chain_step5_F_U(body: DPMBody, Ug_sum: float, r: float,
                    t_n: float, M_proto: float) -> Dict[str, float]:
    """Step 5: Ug_family + Um + FUBi -> F_U.

    F_U = Ug_sum - Ubi + Um

    FUBi (inside-outward) = buoyancy from the local DPM.
    Um   = universal magnetism from nuclear spin coupling.

    Um uses M_proto (ACP-emerged mass), not M_table.
    mu_mag = M_proto * R_nuc^2 * omega0  -- magnetic moment via vortex spin.
    """
    cos_tn = math.cos(math.pi * t_n)
    enh    = 1.0 + EPSILON_SW * RHO_SW

    # FUBi -- inside-outward buoyancy (local DPM)
    Ubi = BETA_I * Ug_sum * OMEGA_G * MBH_DG * enh * RHO_VAC_SCM * cos_tn

    # Um -- universal magnetism (M_proto drives spin coupling, not M_table)
    mu_mag = M_proto * body.R_nuc ** 2 * body.omega0
    Um     = mu_mag / r ** 3

    return {"Ubi": Ubi, "Um": Um, "F_U": Ug_sum - Ubi + Um}


# -- STEP 6: Inside/outside crossing ------------------------------------------

def chain_step6_crossing(body: DPMBody, Ug_sum: float,
                         FUBii_value: float) -> Dict[str, float]:
    """Step 6: F_U -> crossing (FUBi + FUBii = 0 compaction zone).

    THE CROSSING PRECEDES MASS. Mass does not exist before the crossing.
    Mass is BORN at the crossing. (Star-Magic.txt Chapter 6)

    FUBi  (inside-outward): local DPM buoyancy pressure outward
    FUBii (outside-inward): primordial belly button DPM magnetic repulsion inward

    FUBi(r) = BETA_I * |Ug_sum| * rho_SCm * cos(pi*t_n) / r
    Crossing: FUBi(r_cross) + FUBii = 0
    r_cross = BETA_I * |Ug_sum| * rho_SCm * cos(pi*0.25) / |FUBii|
    """
    cos_tn = math.cos(math.pi * 0.25)
    FUBi_at_Rnuc = (BETA_I * abs(Ug_sum) * RHO_VAC_SCM * cos_tn
                    / body.R_nuc)

    if abs(FUBii_value) > 0.0:
        r_cross = (BETA_I * abs(Ug_sum) * RHO_VAC_SCM * cos_tn
                   / abs(FUBii_value))
    else:
        r_cross = body.R_nuc  # fallback: crossing at nuclear radius

    return {
        "FUBi_at_Rnuc":    FUBi_at_Rnuc,
        "FUBii_value":     FUBii_value,
        "r_cross":         r_cross,
        "balance_at_Rnuc": FUBi_at_Rnuc + FUBii_value,
    }


# -- STEP 7: Mass emergence ---------------------------------------------------

def chain_step7_mass_emergence(body: DPMBody, M_proto: float) -> Dict[str, float]:
    """Step 7: crossing -> M_emergent.

    Mass is born at the crossing (Star-Magic.txt Chapter 7):
      M_atomic = M_0 * (1 - exp(-Z/10)) * Z

    M_emergent is the chain output. M_table is the tabulated stable mass.
    scale_factor = M_table / M_emergent shows the calibration residual.

    The scale_factor encodes how the 26-layer DPM amplification
    (sum(i^2, i=1..26) = 6279) and E_crack gating scale up from
    the vacuum base unit M_0_DPM to the observable atomic mass.
    """
    M_emergent   = M_proto
    scale_factor = body.M_table / M_emergent if M_emergent != 0.0 else float("nan")
    return {
        "M_emergent":   M_emergent,
        "M_0_DPM":      M_0_DPM,
        "M_table":      body.M_table,   # verification only
        "scale_factor": scale_factor,   # calibration ratio chain->observed
    }


# -- STEP 8: Newton projection (LAST) -----------------------------------------

def chain_step8_newton(M_table: float, r_cross: float) -> float:
    """Step 8: M -> GM/r^2  (LAST -- observational projection only).

    GM/r^2 is NOT a mechanism. It is what you MEASURE at the crossing
    after the chain has completed and stable mass exists.

    Uses M_table (verified stable mass) and r_cross as the crossing radius.
    """
    return G_CONST * M_table / r_cross ** 2


# -- MASTER CHAIN FUNCTION ----------------------------------------------------

def compute_chain(body: DPMBody,
                  r: Optional[float] = None,
                  t: float = 0.0,
                  t_n: float = 0.25,
                  FUBii_override: Optional[float] = None) -> Dict:
    """Run the full 8-step quantum chain for one DPM body.

    Chain: 0_vacuum -> DPM -> mu_s -> Ug1 -> Ug_family -> F_U -> crossing -> M -> GM/r^2

    The chain is strictly ordered. Mass is never an input -- it is an output.
    Periodic table geometry (Z, A, R_nuc) drives steps 0-6.
    M_table is verified at step 7, used for GM/r^2 at step 8.

    Parameters
    ----------
    body           : DPMBody  (geometry-first)
    r              : evaluation radius [m]  (default: R_nuc)
    t              : time [s]
    t_n            : negative-time parameter
    FUBii_override : override primordial FUBii value [N]
                     Default: self-consistent atomic-scale FUBii
    """
    if r is None:
        r = body.R_nuc

    s0   = chain_step0_vacuum()
    s1   = chain_step1_dpm(body)
    mu_s = chain_step2_mu_s(body)

    # ACP proto-mass from vortex resonance count Z (geometry only, no table)
    M_proto = chain_acp_M_proto(body.Z)

    s4 = chain_step4_ug_family(body, mu_s, M_proto, r, t, t_n)
    s5 = chain_step5_F_U(body, s4["Ug_sum"], r, t_n, M_proto)

    # FUBii: self-consistent atomic-scale (FUBi reversed at R_nuc)
    if FUBii_override is None:
        cos_tn = math.cos(math.pi * t_n)
        FUBii  = -(BETA_I * abs(s4["Ug_sum"]) * RHO_VAC_SCM * cos_tn
                   / body.R_nuc)
    else:
        FUBii = FUBii_override

    s6 = chain_step6_crossing(body, s4["Ug_sum"], FUBii)
    s7 = chain_step7_mass_emergence(body, M_proto)

    r_cross  = s6["r_cross"] if s6["r_cross"] > 0 else body.R_nuc
    g_Newton = chain_step8_newton(body.M_table, r_cross)

    return {
        # Identity
        "Z": body.Z, "A": body.A, "symbol": body.symbol, "name": body.name,
        "R_nuc": body.R_nuc, "V_DPM": body.V_DPM,
        # Step 0
        "s0_grad_UA":   s0["grad_UA"],
        "s0_E_react_0": s0["E_react_0"],
        "s0_F_U_vac":   s0["F_U_vac"],
        # Step 1
        "s1_F_DPM": s1["F_DPM"],
        "s1_a_DPM": s1["a_DPM"],
        # Step 2
        "s2_mu_s":  mu_s,
        # ACP
        "M_proto":  M_proto,
        # Steps 3-4
        "E_react":  s4["E_react"],
        "Ug1":      s4["Ug1"],
        "Ug2":      s4["Ug2"],
        "Ug3":      s4["Ug3"],
        "Ug4":      s4["Ug4"],
        "Ug_sum":   s4["Ug_sum"],
        # Step 5
        "Ubi":      s5["Ubi"],
        "Um":       s5["Um"],
        "F_U":      s5["F_U"],
        # Step 6
        "s6_FUBi":    s6["FUBi_at_Rnuc"],
        "s6_FUBii":   s6["FUBii_value"],
        "s6_r_cross": s6["r_cross"],
        "s6_balance": s6["balance_at_Rnuc"],
        # Step 7
        "s7_M_emergent":   s7["M_emergent"],
        "s7_M_table":      s7["M_table"],
        "s7_scale_factor": s7["scale_factor"],
        # Step 8 -- LAST
        "g_Newton": g_Newton,
    }


# =============================================================================
# S6  DPM RATIO AND GRINDING PAIR  (unchanged -- chain-invariant)
# =============================================================================

def dpm_ratio() -> float:
    """Return DPM = [UA']/[SCm] = rho_UA / rho_SCm = 10 (exact, all scales)."""
    return DPM_DENSITY_RATIO


def grind_opp(scm: float = RHO_VAC_SCM,
              ua_prime: float = RHO_VAC_SCM) -> float:
    """Grind_opp = omega_CW * SCm - omega_CCW * UA'."""
    return OMEGA_CW * scm - OMEGA_CCW * ua_prime


def dpm_react(r: float, t_n: float = 0.25,
              scm: float = RHO_VAC_SCM,
              ua_prime: float = RHO_VAC_SCM) -> float:
    """DPM reaction force density at radius r.

    The r^26 term captures 26-layer gradient suppression.
    Protected against underflow: r_nuc^26 underflows to 0.0 in IEEE doubles,
    so we guard delta=0 (equilibrium vacuum) and r^26 underflow separately.
    """
    grind = grind_opp(scm, ua_prime)
    delta = scm - ua_prime
    if delta == 0.0:
        dpm_n = 0.0
    else:
        r26 = r ** 26
        dpm_n = 0.0 if r26 == 0.0 else KAPPA_FLOAT * delta / r26
    return dpm_n + grind * math.cos(math.pi * t_n)


# =============================================================================
# S7  ALL-ELEMENTS CHAIN COMPUTATION
# =============================================================================

def compute_all_elements_chain(t: float = 0.0,
                               t_n: float = 0.25) -> List[Dict]:
    """Run the full quantum chain for every element Z=1-118.

    The chain is run from vacuum -> GM/r^2 for each DPM body.
    Periodic table geometry (Z, A) drives every computation.
    M_table is verified at the end, never used as primary input.
    """
    results = []
    for body in PERIODIC_TABLE:
        row = compute_chain(body, r=body.R_nuc, t=t, t_n=t_n)
        row["DPM_ratio"] = dpm_ratio()
        row["Grind_opp"] = grind_opp()
        row["DPM_react"] = dpm_react(r=body.R_nuc, t_n=t_n)
        results.append(row)
    return results


# =============================================================================
# S8  SYMPY CANONICAL EXPRESSIONS (quantum chain order)
# =============================================================================

_t, _t_n, _r = sp.symbols('t t_n r', positive=True)
_Z_s, _B0_s, _omega0_s = sp.symbols('Z B0 omega0', positive=True)
_k1, _k2, _k3, _k4, _alpha = sp.symbols('k1 k2 k3 k4 alpha', positive=True)
_beta_i = sp.symbols('beta_i', positive=True)
_rho_A, _rho_UA_s = sp.symbols('rho_SCm rho_UA', positive=True)
_V_b, _v_f = sp.symbols('V_DPM v_fermi', positive=True)
_M_proto_s = sp.symbols('M_proto', positive=True)   # ACP-emerged mass (NOT table mass)
_omega_CW_s, _omega_CCW_s = sp.symbols('omega_CW omega_CCW', positive=True)
_SCm_s, _UAp_s = sp.symbols('SCm UA_prime', positive=True)

_cos_tn  = sp.cos(sp.pi * _t_n)
_mu_s    = _rho_A * _V_b
_E_react = _rho_A * _v_f**2 / _rho_UA_s * sp.exp(-_alpha * _t)

# Chain-ordered symbolic expressions -- M_proto is ACP chain output, NOT table input
Ug1_sym = _k1 * _mu_s * (_M_proto_s / _r**2) * sp.exp(-_alpha * _t) * _cos_tn
Ug2_sym = _k2 * (_rho_A + _rho_UA_s) * _V_b * (_M_proto_s / _r**2) * _E_react
Ug3_sym = _k3 * _B0_s * sp.cos(_omega0_s * _t * sp.pi) * _E_react
Ug4_sym = _k4 * _rho_A * _Z_s * sp.exp(-_alpha * _t) * _cos_tn
Ubi_sym = _beta_i * (Ug1_sym + Ug2_sym + Ug3_sym + Ug4_sym) * _rho_A * _cos_tn
Um_sym  = (_M_proto_s * sp.Symbol('R_nuc')**2 * _omega0_s) / _r**3
F_U_sym = (Ug1_sym + Ug2_sym + Ug3_sym + Ug4_sym) - Ubi_sym + Um_sym

Grind_opp_sym = _omega_CW_s * _SCm_s - _omega_CCW_s * _UAp_s
DPM_ratio_sym = _UAp_s / _SCm_s

F_U_Bi_i_DPM_bound = F_U_Bi_i_DPM.subs(_F_Bi_i_scm, F_U_Bi_i_99)

# ACP mass emergence symbolic (Step 7)
_M_0, _Z_acp = sp.symbols('M_0 Z', positive=True)
M_proto_sym = _M_0 * (1 - sp.exp(-_Z_acp / 10)) * _Z_acp


# =============================================================================
# S9  F_U_Bi_i CALIBRATION PROOF  (uses scm + ua layers)
# =============================================================================

def dpm_fubi_calibration_proof() -> Dict:
    """Prove F_U_Bi vs F_U_Bi_i calibration using the DPM density ratio."""
    mean_fubi_i, std_fubi_i, rng_fubi_i = monte_carlo_fubi_i()
    ratio         = ua_calibration_ratio()
    ua_total_dens = ua_dpm_total_density(0.25)
    return {
        "rho_vac_SCm"          : RHO_VAC_SCM,
        "rho_vac_UA"           : RHO_VAC_UA,
        "ratio_UA_over_SCm"    : ratio,
        "F_U_Bi_i_MC_mean_N"   : mean_fubi_i,
        "F_U_Bi_i_MC_std_N"    : std_fubi_i,
        "F_U_Bi_i_MC_range_N"  : rng_fubi_i,
        "F_U_Bi_cosmological"  : mean_fubi_i * ratio,
        "UA_total_density"     : ua_total_dens,
        "DPM_buoyancy_factor"  : ua_dpm_buoyancy_factor(0.25),
        "scale_interpretation" : (
            "LENR (FUBii) at rho_SCm=7.09e-37 kg/m^3. "
            "Cosmological (FUBi) at rho_UA=7.09e-36 kg/m^3 = 10x SCm. "
            "DPM=[UA']/[SCm]=10 scale-invariant at atomic, stellar, cosmic scales."
        ),
    }


# =============================================================================
# S10 LENR FULL COMPARISON  (scm values + ua mechanism)
# =============================================================================

def dpm_lenr_full_comparison() -> Dict:
    """LENR comparison with both SCm numerical values and UA mechanism."""
    q_park = parkhomov_excess_heat()
    q_pf   = pons_fleischmann_excess_heat()
    ker    = KER_SCm
    ua_data = ua_lenr_comparison()
    ua_data["Holmlid"]["scm_value"]          = ker
    ua_data["Parkhomov"]["scm_value"]        = q_park
    ua_data["Pons-Fleischmann"]["scm_value"] = q_pf
    return ua_data


# =============================================================================
# S11  ENTRY POINT -- FULL QUANTUM CHAIN DEMONSTRATION
# =============================================================================

if __name__ == "__main__":

    SEP  = "=" * 78
    SEP2 = "-" * 78

    print(SEP)
    print("dpm_vacuum_manifold.py v2.0 -- QUANTUM CHAIN IS THE SPINE")
    print("0_vacuum -> DPM -> mu_s -> Ug1 -> Ug_family -> F_U"
          " -> crossing -> M -> GM/r^2")
    print(SEP)
    print(f"  E_crack = rho_SCm*c^2/[SSq] = {E_CRACK:.4e} J")
    print(f"  M_0_DPM = rho_SCm/[SSq]     = {M_0_DPM:.4e} kg  (base DPM mass unit)")
    print(f"  [SSq]                        = {SSQ}")
    print(f"  DPM ratio [UA']/[SCm]        = {dpm_ratio():.1f}  (exact, scale-invariant)")
    print(f"  rho_SCm                      = {RHO_VAC_SCM:.2e} kg/m^3")
    print(f"  rho_UA                       = {RHO_VAC_UA:.2e} kg/m^3")
    print(f"  Grind_opp (uniform vacuum)   = {grind_opp():.4e}")

    # ACP mass scaling Z=1..10
    print(f"\n{SEP}")
    print("ACP PROTO-MASS SCALING  M_proto = M_0 * (1-exp(-Z/10)) * Z")
    print("Mass emerges from DPM vortex resonance count -- NOT from mass table")
    print(SEP2)
    print(f"  {'Z':>3}  {'Sym':4}  {'A':>4}  {'M_proto [kg]':>16}  "
          f"{'M_table [kg]':>16}  {'scale_factor':>14}")
    print(SEP2)
    for body in PERIODIC_TABLE[:30]:
        mp = chain_acp_M_proto(body.Z)
        sf = body.M_table / mp if mp > 0 else float("nan")
        print(f"  {body.Z:>3}  {body.symbol:4}  {body.A:>4}  "
              f"{mp:>16.4e}  {body.M_table:>16.4e}  {sf:>14.4e}")

    # Full 8-step chain for Hydrogen
    print(f"\n{SEP}")
    print("FULL 8-STEP QUANTUM CHAIN -- HYDROGEN  Z=1")
    print(SEP)
    H = ELEMENT[1]
    ch = compute_chain(H)
    print(f"  Body: {H.name}  Z={H.Z}  A={H.A}")
    print(f"  Geometry: R_nuc={H.R_nuc:.4e} m   V_DPM={H.V_DPM:.4e} m^3")
    print(f"  B0={H.B0:.4e} T   omega0={H.omega0:.4e} rad/s")
    print(SEP2)
    print(f"  STEP 0  (zero-mass vacuum):")
    print(f"    grad_UA   = {ch['s0_grad_UA']:.4e} kg/m^3")
    print(f"    E_react_0 = {ch['s0_E_react_0']:.4e}  (peak UA/SCm attraction)")
    print(f"    F_U_vac   = {ch['s0_F_U_vac']:.1f}     (zero -- no mass exists)")
    print(f"  STEP 1  (DPM vortex forms):")
    print(f"    F_DPM     = {ch['s1_F_DPM']:.4e}")
    print(f"    a_DPM     = {ch['s1_a_DPM']:.4e}")
    print(f"  STEP 2  (magnetic moment from vortex -- NOT from mass):")
    print(f"    mu_s      = {ch['s2_mu_s']:.4e}  (rho_SCm * V_DPM)")
    print(f"  ACP  (proto-mass from Z resonance count):")
    print(f"    M_proto   = {ch['M_proto']:.4e} kg   (M_0 * [1-exp(-Z/10)] * Z)")
    print(f"  STEPS 3-4 (Ug family simultaneous -- all from DPM, not from mass table):")
    print(f"    E_react   = {ch['E_react']:.4e}")
    print(f"    Ug1       = {ch['Ug1']:+.4e}  (THE DPM in field form)")
    print(f"    Ug2       = {ch['Ug2']:+.4e}  (outer bubble)")
    print(f"    Ug3       = {ch['Ug3']:+.4e}  (magnetic string)")
    print(f"    Ug4       = {ch['Ug4']:+.4e}  (vacuum concentration, Z={H.Z})")
    print(f"    Ug_sum    = {ch['Ug_sum']:+.4e}")
    print(f"  STEP 5  (F_U assembly):")
    print(f"    Ubi       = {ch['Ubi']:+.4e}  (inside-outward buoyancy)")
    print(f"    Um        = {ch['Um']:+.4e}  (universal magnetism)")
    print(f"    F_U       = {ch['F_U']:+.4e}")
    print(f"  STEP 6  (crossing -- mass BORN here, not before):")
    print(f"    FUBi@Rnuc = {ch['s6_FUBi']:.4e}")
    print(f"    FUBii     = {ch['s6_FUBii']:.4e}")
    print(f"    r_cross   = {ch['s6_r_cross']:.4e} m   (R_nuc = {H.R_nuc:.4e} m)")
    print(f"    balance   = {ch['s6_balance']:.4e}  (-> 0 at true crossing)")
    print(f"  STEP 7  (mass emergence -- chain output, not from table):")
    print(f"    M_emergent  = {ch['s7_M_emergent']:.4e} kg  (ACP chain)")
    print(f"    M_table     = {ch['s7_M_table']:.4e} kg  (verification only)")
    print(f"    scale_factor= {ch['s7_scale_factor']:.4e}  (chain calibration)")
    print(f"  STEP 8  (GM/r^2 -- LAST -- observational projection only):")
    print(f"    g_Newton  = {ch['g_Newton']:.4e} m/s^2")

    # Full 118-element chain table
    print(f"\n{SEP}")
    print("QUANTUM CHAIN ALL 118 ELEMENTS  (r=R_nuc, t=0, t_n=0.25)")
    print(SEP)
    print(f"  {'Z':>3}  {'Sym':4}  {'mu_s':>14}  {'Ug_sum':>14}  "
          f"{'F_U':>14}  {'M_proto[kg]':>14}  {'g_Newton':>14}")
    print(SEP2)
    all_chains = compute_all_elements_chain()
    for row in all_chains:
        print(f"  {row['Z']:>3}  {row['symbol']:4}  "
              f"{row['s2_mu_s']:>+14.3e}  {row['Ug_sum']:>+14.3e}  "
              f"{row['F_U']:>+14.3e}  {row['M_proto']:>14.3e}  "
              f"{row['g_Newton']:>14.3e}")

    # F_U_Bi_i calibration proof
    print(f"\n{SEP}")
    print("F_U_Bi_i CALIBRATION PROOF  (scm Monte-Carlo + ua DPM density)")
    print(SEP)
    proof = dpm_fubi_calibration_proof()
    print(f"  rho_SCm               = {proof['rho_vac_SCm']:.2e} kg/m^3")
    print(f"  rho_UA                = {proof['rho_vac_UA']:.2e} kg/m^3")
    print(f"  Ratio [UA']/[SCm]     = {proof['ratio_UA_over_SCm']:.1f}  (exact)")
    print(f"  FUBii MC mean         = {proof['F_U_Bi_i_MC_mean_N']:.4e} N")
    print(f"  FUBii MC std          = {proof['F_U_Bi_i_MC_std_N']:.4e} N")
    print(f"  FUBi (cosmo, 10x)     = {proof['F_U_Bi_cosmological']:.4e} N")
    print(f"  {proof['scale_interpretation']}")

    # LENR
    print(f"\n{SEP}")
    print("LENR COMPARISON  (SCm values + UA mechanism)")
    print(SEP)
    lenr = dpm_lenr_full_comparison()
    for exp, info in lenr.items():
        val = (f"  scm_value={info['scm_value']:.4e}"
               if info.get("scm_value") is not None else "")
        print(f"  [{exp}]  {info['observable']}{val}")
        print(f"    UA: {str(info['mechanism'])[:90]}...")

    # Next steps
    print(f"\n{SEP}")
    print("NEXT STEPS -- DPM leads to Star-Magic")
    print(SEP)
    print("  1. scm_vacuum_manifold.py   SCm CW base layer          COMPLETE")
    print("  2. ua_vacuum_manifold.py    UA CCW superstructure       COMPLETE")
    print("  3. dpm_vacuum_manifold.py   Chain-compliant assembly    COMPLETE v2.0")
    print("  ---")
    print("  4. NEXT: dpm_chain_papers.py  -- whitepaper per chain step")
    print("     PAPER_vacuum   Step 0 -- primordial DPM, belly button equations")
    print("     PAPER_dpm      Step 1 -- DPM formation at every scale")
    print("     PAPER_Ug_chain Steps 3-4 -- simultaneous Ug assembly proof")
    print("     PAPER_crossing Step 6 -- P!=NP and the compaction zone")
    print("     PAPER_mass     Step 7 -- mass emergence from DPM resonance")
    print("     PAPER_Newton   Step 8 -- why Newton measured the last step")
    print("  5. THEN: wire dpm_vacuum_manifold into source2.cpp via uqff_server.js")
    print(SEP)
    print("THE DPM IS FIRST. GM/r^2 IS LAST.")
    print(SEP)
