# -*- coding: utf-8 -*-
"""
dpm_vacuum_manifold.py
Di-Pseudo-Monopole (DPM) Vacuum Calculator — Complete Assembly

ARCHITECTURE:
  scm_vacuum_manifold.py  →  imported here  (SCm base layer: CW rotation, primordial vacuum)
  ua_vacuum_manifold.py   →  imported here  (UA superstructure: CCW rotation, 4 layers)
  THIS FILE assembles both layers into the complete DPM and:
    1. Computes Ug1, Ug2, Ug3, Ug4, Ubi, Um for any body (stellar or atomic)
    2. Scales the DPM to every atom in the periodic table (Z = 1–118)
    3. Demonstrates Grind_opp and DPM = [UA']/[SCm] at all scales
    4. Provides the F_U_Bi_i Monte-Carlo calibration proof (requires both layers)

DPM FUNDAMENTAL EQUATIONS:
  DPM       = [UA']/[SCm]                                    (ratio, PAPER_411)
  Grind_opp = ω_CW · SCm − ω_CCW · UA'                      (grinding pair reaction)
  F_U       = (Ug1 + Ug2 + Ug3 + Ug4) − Ubi + Um            (total UQFF force)

Ugi COMPONENTS (canonical — SOURCE4 namespace):
  Ug1 = k₁ · μₛ · (M/r²) · e^{−αt} · cos(πtₙ) · (1+δ_def)
        where μₛ = ρ_A · V_body  [magnetic dipole]
  Ug2 = k₂ · (Q_SCm + Q_UA) · (M/r²) · S(r−Rb) · (1+δ_sw·v_sw) · H_SCm · E_react
        [charge-reactivity coupling]
  Ug3 = k₃ · B₀ · cos(ω₀·t·π) · P_core · E_react
        [magnetic string rotation]
  Ug4 = k₄ · ρ_vac · Z · e^{−αt} · cos(πtₙ)
        [vacuum concentration, Z = atomic number]
  Ubi = βᵢ · Σ(Ugᵢ) · Ω_g · (M_bh/d_g) · ρ_A · cos(πtₙ)
        [buoyancy force]
  Um  = μ_mag / r³   where μ_mag = M · R² · ω₀
        [universal magnetism]

PERIODIC TABLE SCALING:
  Each atom is a DPM body: M=atomic mass, R_nuc=nuclear radius, B0=nuclear field,
  omega0=nuclear Larmor frequency.  DPM ratio = 10 at ALL scales.

Author: Daniel T. Murphy  |  dpm_vacuum_manifold.py v1.0
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

import numpy as np
import sympy as sp

# ── Import the two standalone manifold layers ─────────────────────────────────
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

# ─────────────────────────────────────────────────────────────────────────────
# §1  PHYSICAL CONSTANTS
# ─────────────────────────────────────────────────────────────────────────────

HBAR:    float = 1.054571817e-34    # J·s   reduced Planck constant
MU_0:    float = 1.2566370614e-6    # N/A²  vacuum permeability
MU_N:    float = 5.0507837461e-27   # J/T   nuclear magneton
AMU:     float = 1.66053906660e-27  # kg    atomic mass unit
C_LIGHT: float = 2.99792458e8       # m/s   speed of light
R_NUC_0: float = 1.2e-15            # m     nuclear radius constant (R = R0·A^1/3)
K_E:     float = 8.9875517923e9     # N·m²/C²  Coulomb constant

# ─────────────────────────────────────────────────────────────────────────────
# §2  DPM COUPLING CONSTANTS  (Source4 canonical)
# ─────────────────────────────────────────────────────────────────────────────

K1:         float = 1.0     # Ug1 coupling (magnetic dipole)
K2:         float = 1.0     # Ug2 coupling (charge-reactivity)
K3:         float = 1.0     # Ug3 coupling (magnetic string rotation)
K4:         float = 1.0     # Ug4 coupling (vacuum concentration)
ALPHA:      float = KAPPA_FLOAT     # decay rate (= κ = 5e-4)
DELTA_DEF:  float = 0.0     # deformation factor (ground state = 0)
H_SCM:      float = 0.99    # H_SCm coupling constant
EPSILON_SW: float = 0.0     # solar-wind enhancement (atomic context = 0)
RHO_SW:     float = 0.0     # solar-wind density (atomic context = 0)
OMEGA_G:    float = 1.0     # galactic angular factor (normalised for atomic calc)
MBH_DG:     float = 1.0     # M_bh/d_g galactic ratio (normalised)
P_CORE:     float = 1.0     # core pressure factor (dimensionless)

# ── DPM grinding frequencies (from CondensedPhysics4.py OMEGA_CW/CCW) ────────
OMEGA_CW:  float = 2.0 * math.pi * 1.2e10   # rad/s  CW  SCm grinding frequency
OMEGA_CCW: float = 2.0 * math.pi * 8.3e9    # rad/s  CCW UA' grinding frequency

# ─────────────────────────────────────────────────────────────────────────────
# §3  AtomicElement DATACLASS
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class AtomicElement:
    """A periodic-table element treated as a DPM body.

    Parameters
    ----------
    Z       : atomic number
    symbol  : chemical symbol
    name    : element name
    A       : mass number (most abundant / most stable isotope)
    M       : atomic mass [kg]
    R_cov   : covalent radius [m]
    R_nuc   : nuclear radius [m]   = R_NUC_0 · A^(1/3)
    V_nuc   : nuclear volume [m³]
    B0      : nuclear surface magnetic field [T]
    omega0  : nuclear Larmor angular frequency [rad/s]  at B_ref = 1 T
    v_fermi : Fermi velocity proxy [m/s]
    """
    Z:       int
    symbol:  str
    name:    str
    A:       int
    M:       float    # kg
    R_cov:   float    # m
    R_nuc:   float    # m
    V_nuc:   float    # m³
    B0:      float    # T
    omega0:  float    # rad/s
    v_fermi: float    # m/s


def _build_element(Z: int, symbol: str, name: str, A: int,
                   mass_u: float, R_cov_pm: float) -> AtomicElement:
    """Factory: compute derived nuclear parameters and return AtomicElement."""
    M       = mass_u * AMU
    R_cov   = R_cov_pm * 1.0e-12          # pm → m
    R_nuc   = R_NUC_0 * A ** (1.0 / 3.0) # m
    V_nuc   = (4.0 / 3.0) * math.pi * R_nuc ** 3
    # Nuclear surface field: B0 = (μ₀/4π) · 2·Z·μ_N / R_nuc³
    B0      = (MU_0 / (4.0 * math.pi)) * 2.0 * Z * MU_N / R_nuc ** 3
    # Larmor frequency at 1 T: ω₀ = Z · γ_proton  (crude but consistent scaling)
    # γ_proton = 2.675e8 rad/s/T
    omega0  = Z * 2.675e8                  # rad/s  at B=1T reference
    # Fermi velocity: v_F = v_F0 · Z^(1/3)  (free-electron model proxy)
    v_fermi = 0.77e6 * Z ** (1.0 / 3.0)   # m/s
    return AtomicElement(Z=Z, symbol=symbol, name=name, A=A,
                         M=M, R_cov=R_cov, R_nuc=R_nuc, V_nuc=V_nuc,
                         B0=B0, omega0=omega0, v_fermi=v_fermi)


# ─────────────────────────────────────────────────────────────────────────────
# §4  PERIODIC TABLE  (Z = 1–118)
# ─────────────────────────────────────────────────────────────────────────────
# Columns: Z, symbol, name, A (most abundant/stable), mass_u, R_cov_pm

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

# Build all 118 AtomicElement instances
PERIODIC_TABLE: List[AtomicElement] = [
    _build_element(*row) for row in _PT_RAW
]

# Convenience lookup by Z
ELEMENT: Dict[int, AtomicElement] = {el.Z: el for el in PERIODIC_TABLE}


# ─────────────────────────────────────────────────────────────────────────────
# §5  Ugi NUMERICAL COMPUTE FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────
# Evaluation defaults: t=0 (instantaneous), t_n=0.25 (cos(π/4)=√2/2 ≈ 0.707)

def compute_Ug1(el: AtomicElement,
                r: float | None = None,
                t: float = 0.0,
                t_n: float = 0.25) -> float:
    """Ug1 — Magnetic Dipole Force (canonical SOURCE4 equation).

    Ug1 = k₁ · μₛ · (M/r²) · e^{−αt} · cos(πtₙ) · (1+δ_def)
    where μₛ = ρ_A · V_nuc  [vacuum magnetic moment proxy, J/T]

    Parameters
    ----------
    el  : AtomicElement
    r   : evaluation radius [m]  (default: covalent radius)
    t   : time [s]
    t_n : negative-time parameter
    """
    if r is None:
        r = el.R_cov
    mu_s    = RHO_VAC_SCM * el.V_nuc           # [J/T] vacuum magnetic moment
    grad_M  = el.M / (r * r)                   # [kg/m²]
    exp_t   = math.exp(-ALPHA * t)
    cos_tn  = math.cos(math.pi * t_n)
    deform  = 1.0 + DELTA_DEF
    return K1 * mu_s * grad_M * exp_t * cos_tn * deform


def compute_Ug2(el: AtomicElement,
                r: float | None = None,
                t: float = 0.0,
                t_n: float = 0.25) -> float:
    """Ug2 — Charge-Reactivity Coupling (canonical SOURCE4 equation).

    Ug2 = k₂ · (Q_SCm + Q_UA) · (M/r²) · S(r−Rb) · (1+δ_sw·v_sw) · H_SCm · E_react
    where:
      Q_SCm  = ρ_SCm · V_nuc
      Q_UA   = ρ_UA  · V_nuc
      R_b    = R_nuc × 100  (nuclear "heliosphere")
      E_react = ρ_SCm · v² / ρ_UA · exp(−κt)
    """
    if r is None:
        r = el.R_cov
    Q_SCm   = RHO_VAC_SCM * el.V_nuc
    Q_UA    = RHO_VAC_UA  * el.V_nuc
    R_b     = el.R_nuc * 100.0              # nuclear bubble radius
    S_rb    = 1.0 if r > R_b else 0.0      # step function
    sw_fac  = 1.0 + EPSILON_SW * RHO_SW
    E_react = RHO_VAC_SCM * el.v_fermi ** 2 / RHO_VAC_UA * math.exp(-ALPHA * t)
    return K2 * (Q_SCm + Q_UA) * (el.M / (r * r)) * S_rb * sw_fac * H_SCM * E_react


def compute_Ug3(el: AtomicElement,
                r: float | None = None,
                t: float = 0.0,
                t_n: float = 0.25) -> float:
    """Ug3 — Magnetic String Rotation 90° (canonical SOURCE4 equation).

    Ug3 = k₃ · B₀ · cos(ω₀·t·π) · P_core · E_react
    where E_react = ρ_SCm · v² / ρ_UA · exp(−κt)
    """
    if r is None:
        r = el.R_cov
    rot_term = math.cos(el.omega0 * t * math.pi)
    E_react  = RHO_VAC_SCM * el.v_fermi ** 2 / RHO_VAC_UA * math.exp(-ALPHA * t)
    return K3 * el.B0 * rot_term * P_CORE * E_react


def compute_Ug4(el: AtomicElement,
                t: float = 0.0,
                t_n: float = 0.25) -> float:
    """Ug4 — Vacuum Concentration (canonical SOURCE4 equation).

    Ug4 = k₄ · ρ_vac · Z · e^{−αt} · cos(πtₙ)
    At atomic scale, the concentration factor C_concentration = Z (atomic number):
    higher Z → denser vacuum concentration in the nuclear region.
    """
    exp_t  = math.exp(-ALPHA * t)
    cos_tn = math.cos(math.pi * t_n)
    return K4 * RHO_VAC_SCM * float(el.Z) * exp_t * cos_tn


def compute_Ubi(el: AtomicElement,
                Ug_sum: float,
                r: float | None = None,
                t_n: float = 0.25) -> float:
    """Ubi — Buoyancy Force (canonical SOURCE4 equation).

    Ubi = βᵢ · Ug_sum · Ω_g · (M_bh/d_g) · (1 + ε_sw·ρ_sw) · ρ_A · cos(πtₙ)
    For atomic context: Ω_g = 1, M_bh/d_g = 1, ε_sw = 0.
    """
    enhancement = 1.0 + EPSILON_SW * RHO_SW
    cos_tn      = math.cos(math.pi * t_n)
    return BETA_I * Ug_sum * OMEGA_G * MBH_DG * enhancement * RHO_VAC_SCM * cos_tn


def compute_Um(el: AtomicElement,
               r: float | None = None) -> float:
    """Um — Universal Magnetism (canonical SOURCE4 equation).

    Um = μ_mag / r³   where μ_mag = M · R_nuc² · ω₀

    The nuclear magnetic moment proxy μ_mag captures the spin-mass coupling.
    """
    if r is None:
        r = el.R_cov
    mu_mag = el.M * el.R_nuc ** 2 * el.omega0
    return mu_mag / (r ** 3)


def compute_F_U(el: AtomicElement,
                r: float | None = None,
                t: float = 0.0,
                t_n: float = 0.25) -> Dict[str, float]:
    """Compute complete F_U = (Ug1+Ug2+Ug3+Ug4) − Ubi + Um for one element.

    Returns a dict with all components plus the total.
    """
    if r is None:
        r = el.R_cov
    ug1 = compute_Ug1(el, r, t, t_n)
    ug2 = compute_Ug2(el, r, t, t_n)
    ug3 = compute_Ug3(el, r, t, t_n)
    ug4 = compute_Ug4(el, t, t_n)
    ug_sum = ug1 + ug2 + ug3 + ug4
    ubi = compute_Ubi(el, ug_sum, r, t_n)
    um  = compute_Um(el, r)
    F_U = ug_sum - ubi + um
    return {
        "Ug1": ug1,
        "Ug2": ug2,
        "Ug3": ug3,
        "Ug4": ug4,
        "Ug_sum": ug_sum,
        "Ubi": ubi,
        "Um":  um,
        "F_U": F_U,
    }


# ─────────────────────────────────────────────────────────────────────────────
# §6  DPM RATIO AND GRINDING PAIR
# ─────────────────────────────────────────────────────────────────────────────

def dpm_ratio() -> float:
    """Return DPM = [UA']/[SCm] = ρ_vac_UA / ρ_vac_SCm = 10 (exact, all scales)."""
    return DPM_DENSITY_RATIO


def grind_opp(scm: float = RHO_VAC_SCM,
              ua_prime: float = RHO_VAC_SCM) -> float:
    """Grind_opp = ω_CW · SCm − ω_CCW · UA'

    The CW grinding component (SCm) minus CCW grinding component (UA').
    Net positive → CW dominates; net negative → CCW dominates.

    Parameters
    ----------
    scm      : SCm vacuum density value [kg/m³]  (default: RHO_VAC_SCM)
    ua_prime : UA' vacuum density value [kg/m³]  (default: UA' = RHO_VAC_SCM)
    """
    return OMEGA_CW * scm - OMEGA_CCW * ua_prime


def dpm_react(r: float, t_n: float = 0.25,
              scm: float = RHO_VAC_SCM,
              ua_prime: float = RHO_VAC_SCM) -> float:
    """DPM reaction force density.

    DPM_react = κ · (DPM_n_SCm − DPM_s_UA') / r^{26}
              + Grind_opp · cos(π·t_n)   [26D oscillatory term]

    The r^{26} denominator encodes the 26-dimensional projection.
    cos(π·t_n) is the negative-time modulation.
    """
    grind = grind_opp(scm, ua_prime)
    dpm_n = KAPPA_FLOAT * (scm - ua_prime) / (r ** 26)
    osc   = grind * math.cos(math.pi * t_n)
    return dpm_n + osc


# ─────────────────────────────────────────────────────────────────────────────
# §7  SYMPY CANONICAL EXPRESSIONS
# ─────────────────────────────────────────────────────────────────────────────

_t, _t_n, _r, _M_s, _B0_s, _omega0_s, _Z_s = sp.symbols(
    't t_n r M B0 omega0 Z', positive=True)
_k1, _k2, _k3, _k4, _alpha = sp.symbols('k1 k2 k3 k4 alpha', positive=True)
_beta_i = sp.symbols('beta_i', positive=True)
_rho_A  = sp.symbols('rho_A', positive=True)
_V_b    = sp.symbols('V_b', positive=True)
_v_f    = sp.symbols('v_f', positive=True)
_rho_UA_s = sp.symbols('rho_UA', positive=True)
_omega_CW_s, _omega_CCW_s = sp.symbols('omega_CW omega_CCW', positive=True)
_SCm_s, _UAp_s = sp.symbols('SCm UA_prime', positive=True)

_cos_tn  = sp.cos(sp.pi * _t_n)
_mu_s    = _rho_A * _V_b
_E_react = _rho_A * _v_f**2 / _rho_UA_s * sp.exp(-_alpha * _t)

# Canonical Ugi symbolic forms
Ug1_sym = _k1 * _mu_s * (_M_s / _r**2) * sp.exp(-_alpha * _t) * _cos_tn
Ug2_sym = _k2 * (_rho_A + _rho_UA_s) * _V_b * (_M_s / _r**2) * _E_react
Ug3_sym = _k3 * _B0_s * sp.cos(_omega0_s * _t * sp.pi) * _E_react
Ug4_sym = _k4 * _rho_A * _Z_s * sp.exp(-_alpha * _t) * _cos_tn
Ubi_sym = _beta_i * (Ug1_sym + Ug2_sym + Ug3_sym + Ug4_sym) * _rho_A * _cos_tn
Um_sym  = (_M_s * sp.Symbol('R_nuc')**2 * _omega0_s) / _r**3

F_U_sym = (Ug1_sym + Ug2_sym + Ug3_sym + Ug4_sym) - Ubi_sym + Um_sym

# DPM grinding pair
Grind_opp_sym = _omega_CW_s * _SCm_s - _omega_CCW_s * _UAp_s
DPM_ratio_sym = _UAp_s / _SCm_s

# Full DPM buoyancy (binding the placeholder from ua_vacuum_manifold)
# F_U_Bi_i_DPM (from ua) with F_Bi_i_scm → actual F_U_Bi_i_99 sum from scm
F_U_Bi_i_DPM_bound = F_U_Bi_i_DPM.subs(_F_Bi_i_scm, F_U_Bi_i_99)


# ─────────────────────────────────────────────────────────────────────────────
# §8  PERIODIC TABLE DPM SCALING  (all 118 elements)
# ─────────────────────────────────────────────────────────────────────────────

def compute_all_elements(t: float = 0.0,
                         t_n: float = 0.25) -> List[Dict]:
    """Compute DPM Ugi, Ubi, Um, F_U for every element Z=1–118.

    Returns
    -------
    list of dicts, one per element, containing all force components + metadata.
    """
    results = []
    for el in PERIODIC_TABLE:
        fu = compute_F_U(el, r=el.R_cov, t=t, t_n=t_n)
        grind = grind_opp()                          # universal — does not vary per element
        dpm_r = dpm_react(r=el.R_cov, t_n=t_n)
        results.append({
            "Z":      el.Z,
            "symbol": el.symbol,
            "name":   el.name,
            "A":      el.A,
            "M_kg":   el.M,
            "R_cov_m":el.R_cov,
            "R_nuc_m":el.R_nuc,
            "B0_T":   el.B0,
            "omega0": el.omega0,
            **fu,
            "Grind_opp": grind,
            "DPM_react": dpm_r,
            "DPM_ratio": dpm_ratio(),
        })
    return results


# ─────────────────────────────────────────────────────────────────────────────
# §9  F_U_Bi_i CALIBRATION PROOF  (requires both scm + ua layers)
# ─────────────────────────────────────────────────────────────────────────────

def dpm_fubi_calibration_proof() -> Dict[str, object]:
    """Prove F_U_Bi vs F_U_Bi_i calibration using the DPM density ratio.

    F_U_Bi   (inside→outside, cosmological)  ∝ ρ_vac_UA
    F_U_Bi_i (outside→inside, LENR/atomic)   ∝ ρ_vac_SCm
    Ratio = ρ_vac_UA / ρ_vac_SCm = 10  (exact, scale-invariant)

    Returns
    -------
    dict with Monte-Carlo statistics and interpretation.
    """
    mean_fubi_i, std_fubi_i, rng_fubi_i = monte_carlo_fubi_i()
    ratio = ua_calibration_ratio()
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
            "LENR (F_U_Bi_i) at ρ_vac_SCm = 7.09e-37 kg/m³. "
            "Cosmology (F_U_Bi) at ρ_vac_UA = 7.09e-36 kg/m³ = 10× SCm. "
            "DPM ratio [UA']/[SCm] = 10 is scale-invariant: "
            "identical at atomic, stellar, and cosmological scales."
        ),
    }


# ─────────────────────────────────────────────────────────────────────────────
# §10 LENR FULL COMPARISON  (scm values + ua mechanism — both layers)
# ─────────────────────────────────────────────────────────────────────────────

def dpm_lenr_full_comparison() -> Dict[str, Dict]:
    """Full LENR comparison with both SCm numerical values and UA mechanism.

    Merges scm_vacuum_manifold numerical data with ua_vacuum_manifold
    UA-layer mechanism explanations.
    """
    q_park = parkhomov_excess_heat()
    q_pf   = pons_fleischmann_excess_heat()
    ker    = KER_SCm
    ua_data = ua_lenr_comparison()

    # Inject scm values into ua entries
    ua_data["Holmlid"]["scm_value"]         = ker
    ua_data["Parkhomov"]["scm_value"]       = q_park
    ua_data["Pons-Fleischmann"]["scm_value"] = q_pf
    return ua_data


# ─────────────────────────────────────────────────────────────────────────────
# §11  ENTRY POINT — FULL DPM DEMONSTRATION
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":

    SEP  = "=" * 78
    SEP2 = "-" * 78

    # ── §1  DPM fundamental constants ────────────────────────────────────────
    print(SEP)
    print("§1  DPM FUNDAMENTAL CONSTANTS")
    print(SEP)
    print(f"  DPM ratio [UA']/[SCm]    = {dpm_ratio():.1f}  (exact, scale-invariant)")
    print(f"  ρ_vac_SCm (CW  base)     = {RHO_VAC_SCM:.2e} kg/m³")
    print(f"  ρ_vac_UA  (CCW super)    = {RHO_VAC_UA:.2e} kg/m³")
    print(f"  ω_CW  grinding freq      = {OMEGA_CW:.4e} rad/s")
    print(f"  ω_CCW grinding freq      = {OMEGA_CCW:.4e} rad/s")
    grind = grind_opp()
    print(f"  Grind_opp = ω_CW·SCm − ω_CCW·UA' = {grind:.6e}")

    # ── §2  Sympy canonical Ugi expressions ──────────────────────────────────
    print(f"\n{SEP}")
    print("§2  CANONICAL Ugi SYMBOLIC EXPRESSIONS")
    print(SEP)
    print(f"  Ug1 = {sp.latex(Ug1_sym)}")
    print(f"  Ug2 = {sp.latex(Ug2_sym)}")
    print(f"  Ug3 = {sp.latex(Ug3_sym)}")
    print(f"  Ug4 = {sp.latex(Ug4_sym)}")
    print(f"  Ubi = β_i · (Ug1+Ug2+Ug3+Ug4) · ρ_A · cos(πt_n)")
    print(f"  Um  = (M·R_nuc²·ω₀) / r³")
    print(f"  F_U = Ug1+Ug2+Ug3+Ug4 − Ubi + Um")
    print(f"\n  DPM grinding pair (sympy):")
    print(f"  Grind_opp = {sp.latex(Grind_opp_sym)}")
    print(f"  DPM_ratio = {sp.latex(DPM_ratio_sym)}")
    print(f"\n  F_U_Bi_i_DPM (bound):")
    print(f"  = F_U_Bi_i_99 × UA_total  (scm sum × ua superstructure density)")

    # ── §3  Ugi for hydrogen (Z=1) ────────────────────────────────────────────
    print(f"\n{SEP}")
    print("§3  DPM Ugi DEMONSTRATION — HYDROGEN  (Z=1)")
    print(SEP)
    H = ELEMENT[1]
    fu_H = compute_F_U(H)
    print(f"  Element : {H.name} (Z={H.Z}, A={H.A})")
    print(f"  M       = {H.M:.4e} kg")
    print(f"  R_cov   = {H.R_cov:.4e} m   (31 pm)")
    print(f"  R_nuc   = {H.R_nuc:.4e} m   (proton radius ~1.2 fm)")
    print(f"  B0      = {H.B0:.4e} T    (nuclear surface field)")
    print(f"  ω₀      = {H.omega0:.4e} rad/s  (Larmor at 1T)")
    print(SEP2)
    for k, v in fu_H.items():
        print(f"  {k:8s} = {v:+.6e}")

    # ── §4  Periodic table DPM scaling (all 118) ─────────────────────────────
    print(f"\n{SEP}")
    print("§4  PERIODIC TABLE DPM SCALING  (Z=1–118,  t_n=0.25)")
    print(SEP)
    print(f"  {'Z':>3}  {'Sym':4}  {'Ug1':>14}  {'Ug2':>14}  {'Ug3':>14}  "
          f"{'Ug4':>14}  {'Ubi':>14}  {'Um':>14}  {'F_U':>14}")
    print(SEP2)
    all_el = compute_all_elements()
    for row in all_el:
        print(f"  {row['Z']:>3}  {row['symbol']:4}  "
              f"{row['Ug1']:>+14.4e}  {row['Ug2']:>+14.4e}  "
              f"{row['Ug3']:>+14.4e}  {row['Ug4']:>+14.4e}  "
              f"{row['Ubi']:>+14.4e}  {row['Um']:>+14.4e}  "
              f"{row['F_U']:>+14.4e}")

    # ── §5  F_U_Bi_i Monte-Carlo + calibration proof ─────────────────────────
    print(f"\n{SEP}")
    print("§5  F_U_Bi_i MONTE-CARLO CALIBRATION PROOF  (scm + ua combined)")
    print(SEP)
    proof = dpm_fubi_calibration_proof()
    print(f"  ρ_vac_SCm              = {proof['rho_vac_SCm']:.2e} kg/m³")
    print(f"  ρ_vac_UA               = {proof['rho_vac_UA']:.2e} kg/m³")
    print(f"  Ratio [UA']/[SCm]      = {proof['ratio_UA_over_SCm']:.1f}  (exact)")
    print(f"  F_U_Bi_i MC mean       = {proof['F_U_Bi_i_MC_mean_N']:.4e} N")
    print(f"  F_U_Bi_i MC std        = {proof['F_U_Bi_i_MC_std_N']:.4e} N")
    print(f"  F_U_Bi (cosmo, 10×)    = {proof['F_U_Bi_cosmological']:.4e} N")
    print(f"  UA total density       = {proof['UA_total_density']:.4e} kg/m³")
    print(f"  DPM buoyancy factor    = {proof['DPM_buoyancy_factor']:.6f}")
    print(f"  Interpretation:")
    print(f"    {proof['scale_interpretation']}")

    # ── §6  LENR full comparison ──────────────────────────────────────────────
    print(f"\n{SEP}")
    print("§6  LENR FULL COMPARISON  (SCm + UA dual-layer)")
    print(SEP)
    lenr = dpm_lenr_full_comparison()
    for exp, info in lenr.items():
        val = (f"  scm_value={info['scm_value']:.4e}"
               if info["scm_value"] is not None else "")
        print(f"  [{exp}]  {info['observable']}{val}")
        print(f"    UA mechanism: {info['mechanism'][:90]}...")

    # ── §7  DPM react at selected radii ──────────────────────────────────────
    print(f"\n{SEP}")
    print("§7  DPM_react  (grinding pair reaction at selected radii)")
    print(SEP)
    for r_val in [1e-15, 1e-12, 1e-10, 1e-6, 1.0, 1e6, 1e11]:
        try:
            dr = dpm_react(r=r_val, t_n=0.25)
            print(f"  r = {r_val:.1e} m  →  DPM_react = {dr:.6e}")
        except (OverflowError, ZeroDivisionError):
            print(f"  r = {r_val:.1e} m  →  DPM_react = (overflow at r^26)")

    print(f"\n{SEP}")
    print("✅  dpm_vacuum_manifold.py COMPLETE")
    print("    scm_vacuum_manifold (SCm CW base) + ua_vacuum_manifold (UA CCW)")
    print("    Ugi computed for all 118 elements.  DPM ratio = 10 at all scales.")
    print(SEP)
