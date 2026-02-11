"""
CondensedPhysics_InputData.py
=============================

INPUT DATA REPOSITORY for CondensedPhysics.py Calculator

This file stores OBSERVATIONAL and EMPIRICAL parameters used as inputs
to the UQFF equations in CondensedPhysics.py.

ARCHITECTURE:
    source2.cpp (HEAD) → API Fetch → bodies_YYYYMMDD_HHMMSS.csv
                                          ↓
    CondensedPhysics_InputData.py   →   CondensedPhysics.py (CALCULATOR)
                                          ↓
                                   CondensedPhysics_OutputData.py

RULES:
    - This file contains INPUT PARAMETERS ONLY (no equations)
    - Parameters are organized by observational source/event
    - All values must have units documented
    - Data can be updated as new observations become available

Copyright © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import numpy as np
from typing import Dict, Any

# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICAL CONSTANTS
# ═══════════════════════════════════════════════════════════════════════════════

SOLAR_MASS_KG = 1.989e30  # kg
SPEED_OF_LIGHT = 2.998e8  # m/s
G = 6.674e-11  # m³/kg/s²


# ═══════════════════════════════════════════════════════════════════════════════
# GW170817 BINARY NEUTRON STAR MERGER
# Source: LIGO/Virgo detection Aug 17, 2017
# References: arXiv:2503.17445, A&A aa52290-24.pdf, MNRAS 2025 analyses
# ═══════════════════════════════════════════════════════════════════════════════

GW170817_PARAMS = {
    # Event metadata
    'event_name': 'GW170817',
    'event_date': '2017-08-17',
    'detection': 'LIGO/Virgo gravitational waves + EM counterpart AT2017gfo',
    
    # Binary neutron star parameters
    'M_ns_total': 2.8,          # M_☉ - total mass of binary
    'M_ns_total_kg': 2.8 * SOLAR_MASS_KG,  # kg
    'M_ns_1': 1.36,             # M_☉ - primary NS mass estimate
    'M_ns_2': 1.17,             # M_☉ - secondary NS mass estimate (chirp mass derived)
    
    # Ejecta parameters (2025 NR simulations)
    'M_ej_total': 0.05,         # M_☉ - total ejecta mass
    'M_ej_total_kg': 0.05 * SOLAR_MASS_KG,  # kg
    'M_ej_dynamical': 0.02,     # M_☉ - dynamical ejecta (~40%)
    'M_ej_wind': 0.04,          # M_☉ - wind ejecta (~80% of wind+dyn total)
    'M_ej_fraction_dynamical': 0.40,  # dimensionless - 40% dynamical fraction
    
    # Velocity parameters
    'v_ej': 0.1,                # c - ejecta velocity (0.1c to 0.3c range)
    'v_ej_ms': 0.1 * SPEED_OF_LIGHT,  # m/s
    'v_ej_min': 0.1,            # c - minimum observed
    'v_ej_max': 0.3,            # c - maximum observed (fast tail)
    
    # Neutron star density
    'rho_ns': 1e15,             # kg/m³ - NS core density
    'rho_ns_surface': 1e14,     # kg/m³ - NS surface density
    
    # Electron fraction (r-process critical)
    'Ye_midplane': 0.1,         # dimensionless - neutron-rich midplane
    'Ye_outflow': 0.2,          # dimensionless - less neutron-rich outflow
    'Ye_effective': 0.15,       # dimensionless - effective average
    'Ye_rprocess_threshold': 0.25,  # Ye < 0.25 required for r-process
    
    # R-process nucleosynthesis
    'r_process_solar': 0.95,    # dimensionless - 95% solar abundance match
    'A_min': 140,               # mass number - r-process lower bound
    'A_max': 254,               # mass number - predicted from exponential term
    'A_predicted': 254,         # mass number - UQFF prediction
    
    # Neutrino parameters
    'E_nu_total_erg': 1e53,     # erg - total neutrino energy
    'E_nu_total_J': 1e53 * 1e-7,  # J
    'neutrino_outflow': 0.70,   # dimensionless - 70% outflow
    'neutrino_inflow': 0.30,    # dimensionless - 30% inflow
    
    # UQFF calibrated parameters
    'beta_i': 0.61,             # dimensionless - buoyancy opposition strength
    'Ub_i_magnitude': 1e27,     # J/m³ - buoyancy energy density (scaled from Sun)
    'omega': np.pi,             # rad/s - modulation frequency
    
    # Distance and redshift
    'distance_Mpc': 40,         # Mpc - distance to NGC 4993 host
    'distance_m': 40 * 3.086e22,  # m
    'redshift': 0.0099,         # z - cosmological redshift
    
    # Kilonova AT2017gfo
    'kilonova_peak_mag': -16,   # mag - peak absolute magnitude
    'kilonova_decay_days': 7,   # days - characteristic decay time
    
    # Verification sources
    'sources': [
        'LIGO/Virgo Collaboration 2017',
        'arXiv:2503.17445 (2025 NR simulations)',
        'A&A aa52290-24.pdf',
        'MNRAS 510 L7 (2022)',
        'Frontiers Physics 8:355 (2020)'
    ]
}


# ═══════════════════════════════════════════════════════════════════════════════
# R-PROCESS OUTFLOW GENERIC PARAMETERS
# For use with any NS merger or collapsar event
# ═══════════════════════════════════════════════════════════════════════════════

RPROCESS_DEFAULTS = {
    # UQFF calibrated constants
    'beta_i': 0.61,             # dimensionless - opposition strength
    'kappa': 0.0005,            # 1/day - decay constant
    'SSq': 0.57,                # dimensionless - calibrated [SSq]
    
    # Vacuum densities
    'rho_vac_SCm': 7.09e-37,    # J/m³
    'rho_vac_UA': 7.09e-36,     # J/m³
    'lambda_vac_sw': 1e-10,     # dimensionless - solar wind vacuum coupling
    
    # Default NS parameters
    'M_ns_default': 1.4,        # M_☉ - canonical NS mass
    'R_ns_default': 12e3,       # m - typical NS radius 12 km
    'rho_ns_default': 1e15,     # kg/m³
    
    # R-process thresholds
    'Ye_rprocess_max': 0.25,    # Ye < 0.25 for r-process
    'A_rprocess_min': 56,       # Fe-56 is starting point
    'A_rprocess_heavy': 140,    # Heavy r-process nuclei
}


# ═══════════════════════════════════════════════════════════════════════════════
# JCAP COSMOLOGICAL PARAMETERS
# Document: UQFF proof set verification of ρ_vac ratios for JCAP DM density_
#           λ_vac alignment_28Sept2025.docx
# References:
#   - JCAP01(2025)021: Solar DM density ~0.47 GeV/cm³
#   - JCAP07(2025)033: Primordial DM ~10^{-26} J/m³
#   - arXiv:2505.17828, 2408.00822
# ═══════════════════════════════════════════════════════════════════════════════

JCAP_COSMOLOGY_PARAMS = {
    # Cosmological constant
    'Lambda': 1.1e-52,              # m⁻² - cosmological constant
    'Lambda_energy': 5.35e-10,      # J/m³ - dark energy density ρ_Λ c²
    
    # Dark energy density
    'rho_Lambda': 5.96e-27,         # kg/m³ - dark energy mass density
    'rho_Lambda_J': 5.35e-10,       # J/m³ - dark energy ~10^{-9} J/m³
    
    # Dark matter densities (JCAP papers)
    'rho_DM_local': 5.35e-25,       # J/m³ - local DM halo ~0.3 GeV/cm³
    'rho_DM_solar': 8.4e-25,        # J/m³ - Solar DM ~0.47 GeV/cm³ (JCAP01(2025)021)
    'rho_DM_primordial': 1e-26,     # J/m³ - primordial DM (JCAP07(2025)033)
    'rho_DM_GeV_cm3': 0.3,          # GeV/cm³ - standard local DM
    'rho_DM_range_min': 0.3,        # GeV/cm³
    'rho_DM_range_max': 0.5,        # GeV/cm³
    
    # Vacuum energy density components (UQFF)
    'rho_vac_SCm': 7.09e-37,        # J/m³ - [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,         # J/m³ - [UA] vacuum density
    'rho_vac_A': 1e-23,             # J/m³ - Universal Aether vacuum density
    'rho_vac_Ui': 2.84e-36,         # J/m³ - Ui vacuum component
    
    # λ_vac sum parameters
    'E_0': 1e-20,                   # J - base energy for E_i = E_0 × 10^i
    'n_layers': 26,                 # number of vacuum layers
    'lambda_vac_cosmic': 1e-9,      # J/m³ - cosmic scale λ_vac ≈ ρ_Λ
    
    # Ratio verification targets
    'ratio_SCm_to_lambda_vac': 7.09e-28,  # ρ_vac,[SCm] / λ_vac
    'ratio_UA_to_lambda_vac': 7.09e-27,   # ρ_vac,[UA] / λ_vac
    'ratio_target_order': -28,            # log10 order of ratios (~10^{-28})
    'ratio_document_order': -38,          # document states ~10^{-38} (may be specific sub-ratio)
    
    # Physical constants for DE calculation
    'c': SPEED_OF_LIGHT,            # m/s
    'G': 6.674e-11,                 # m³/kg/s²
    
    # Inertia fraction bounds
    'f_i_sum': 1.0,                 # Σ f_i = 1
    'f_i_min': 0.0,                 # minimum fraction
    'f_i_max': 1.0,                 # maximum fraction
}


# ═══════════════════════════════════════════════════════════════════════════════
# RACS J0320-35 QUASAR JET PARAMETERS
# Document: UQFF's Navier-Stokes demonstrates RACS J0320-35's fluid jets_28Sept2025.docx
# References:
#   - Chandra X-ray Observatory 2025 observations
#   - NASA release: chandra.si.edu/photo/2025/red6/
# ═══════════════════════════════════════════════════════════════════════════════

RACS_J0320_35_PARAMS = {
    # Source identification
    'name': 'RACS J0320-35',
    'type': 'Quasar / Super-Eddington SMBH',
    'discovery': 'Chandra 2025',
    
    # Distance and redshift
    'redshift': 6.5,                        # z ~ 6.5
    'distance_ly': 12.8e9,                  # 12.8 billion light-years
    'distance_m': 12.8e9 * 9.461e15,        # m
    
    # Black hole parameters
    'M_bh_Msun': 4e8,                       # ~400 million M_sun (estimate)
    'M_bh_kg': 4e8 * SOLAR_MASS_KG,         # kg
    'eddington_ratio': 2.4,                 # 2.4× Eddington limit
    'accretion_rate_Msun_yr': 3000,         # ~3000 M_sun/year
    
    # Jet parameters (Navier-Stokes verification)
    'v_jet_c': 0.99,                        # ~0.99c relativistic jets
    'v_jet_ms': 0.99 * SPEED_OF_LIGHT,      # m/s
    'v_SCm': 1e8,                           # m/s - [SCm] expulsion velocity (~c/3)
    'jet_asymmetry_ratio': 1.5,             # Observed asymmetry (1.5-2)
    'jet_length_kpc': 100,                  # ~100 kpc typical jet length
    'jet_length_m': 100 * 3.086e19,         # m
    
    # Fluid dynamics parameters
    'rho_jet_kgm3': 1e-21,                  # Jet density (relativistic plasma)
    'mu_viscosity': 1e-11,                  # Dynamic viscosity (Pa·s)
    'nu_kinematic': 1e10,                   # Kinematic viscosity (m²/s)
    'reynolds_number': 1e15,                # Re >> 1 (turbulent)
    
    # UQFF buoyancy parameters
    'beta_i': 0.61,                         # Buoyancy opposition strength
    'omega': np.pi,                         # rad/s - modulation frequency
    'gamma_decay': 0.0005,                  # day⁻¹ - decay rate
    'delta_sw': 1.0,                        # Solar wind modulation
    'lambda_vac_sw': 1e-10,                 # Vacuum wavelength
    
    # Time-reversal zone (TRZ) parameters
    # t_n2=0.25 gives cos(π×0.25)=cos(π/4)=0.707, asymmetry ratio ~1.414
    't_n1': 0.0,                            # Jet 1 normalized time
    't_n2': 0.25,                           # Jet 2 normalized time (π/4 phase shift)
    'delta_t': 0.25,                        # Phase difference causing asymmetry
    
    # Energy parameters
    'E_react_Wm3': 1e46,                    # Reactive energy density W/m³
    'f_quasi': 0.1,                         # Quasi-equilibrium factor
    'jet_luminosity_W': 1e40,               # Typical jet luminosity
    
    # CRP parameters (cosmic ray propagation)
    'p_max_eV': 1e16,                       # CRP momentum cutoff
    'spectral_index': 2.2,                  # Power law n(p) ∝ p^{-2.2}
    'D_E_exponent': 0.5,                    # Diffusion D ∝ E^0.5
    
    # Physical constants
    'c': SPEED_OF_LIGHT,
    'G': 6.674e-11,
    
    # UQFF vacuum components
    'rho_vac_UA': 7.09e-36,                 # [UA] vacuum density
    'rho_vac_SCm': 7.09e-37,                # [SCm] vacuum density
}


# ═══════════════════════════════════════════════════════════════════════════════
# HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def get_event_params(event_name: str) -> Dict[str, Any]:
    """
    Retrieve input parameters for a specific astrophysical event.
    
    Parameters
    ----------
    event_name : str
        Name of event (e.g., 'GW170817', 'JCAP_COSMOLOGY', 'RACS_J0320_35')
    
    Returns
    -------
    dict : Event parameters
    """
    events = {
        'GW170817': GW170817_PARAMS,
        'JCAP_COSMOLOGY': JCAP_COSMOLOGY_PARAMS,
        'RACS_J0320_35': RACS_J0320_35_PARAMS,
        'YANG_MILLS': YANG_MILLS_PARAMS,
        'RIEMANN_HYPOTHESIS': RIEMANN_HYPOTHESIS_PARAMS,
        'P_VS_NP': P_VS_NP_PARAMS,
        'DUST_YIELD': DUST_YIELD_PARAMS,
        'SGRA_STAR': SGRA_STAR_PARAMS,
        'SHEAR_CHI_SQUARED': SHEAR_CHI_SQUARED_PARAMS,
        'TIME_ASYMMETRY': TIME_ASYMMETRY_PARAMS,
        'NUCLEAR_BINDING_SHELL': NUCLEAR_BINDING_SHELL_PARAMS,
        'SOLAR_WIND_PARKER_PROBE': SOLAR_WIND_PARKER_PROBE_PARAMS,
        'ALPHA_BEC_LENR': ALPHA_BEC_LENR_PARAMS,
        'UQFF_MASTER_FRAMEWORK': UQFF_MASTER_FRAMEWORK_PARAMS,
    }
    
    if event_name in events:
        return events[event_name]
    else:
        raise ValueError(f"Unknown event: {event_name}. Available: {list(events.keys())}")


def create_rprocess_params(
    M_ns: float,
    M_ej: float,
    v_ej: float,
    Ye: float,
    rho_ns: float = None,
    beta_i: float = None
) -> Dict[str, Any]:
    """
    Create parameter dictionary for generic r-process outflow calculation.
    
    Parameters
    ----------
    M_ns : float
        Neutron star mass in solar masses
    M_ej : float
        Ejecta mass in solar masses
    v_ej : float
        Ejecta velocity in units of c
    Ye : float
        Electron fraction (dimensionless)
    rho_ns : float, optional
        NS density in kg/m³ (default: 1e15)
    beta_i : float, optional
        Buoyancy opposition strength (default: 0.61)
    
    Returns
    -------
    dict : Complete parameter set for RProcessOutflowModel
    """
    defaults = RPROCESS_DEFAULTS
    
    return {
        'M_ns': M_ns,
        'M_ns_kg': M_ns * SOLAR_MASS_KG,
        'M_ej': M_ej,
        'M_ej_kg': M_ej * SOLAR_MASS_KG,
        'v_ej': v_ej,
        'v_ej_ms': v_ej * SPEED_OF_LIGHT,
        'Ye': Ye,
        'rho_ns': rho_ns if rho_ns else defaults['rho_ns_default'],
        'beta_i': beta_i if beta_i else defaults['beta_i'],
        'kappa': defaults['kappa'],
        'SSq': defaults['SSq'],
        'rho_vac_SCm': defaults['rho_vac_SCm'],
        'rho_vac_UA': defaults['rho_vac_UA'],
    }


def get_gw170817_rprocess_params() -> Dict[str, Any]:
    """
    Get parameters formatted for RProcessOutflowModel.compute_complete_rprocess_outflow().
    
    Returns dictionary compatible with:
        model.compute_complete_rprocess_outflow(**get_gw170817_rprocess_params())
    """
    p = GW170817_PARAMS
    return {
        'M_ns': p['M_ns_total_kg'],
        'M_ej': p['M_ej_total_kg'],
        'd': 1e7,  # Characteristic merger separation ~10 km
        'rho': p['rho_ns'],
        'beta_i': p['beta_i'],
        't_n': 0.0
    }


def get_vacuum_alignment_params() -> Dict[str, Any]:
    """
    Get parameters formatted for VacuumDensityAlignmentModel.
    
    Returns dictionary compatible with:
        model.compute_complete_alignment(**get_vacuum_alignment_params())
    """
    p = JCAP_COSMOLOGY_PARAMS
    return {
        'rho_vac_components': {
            'SCm': p['rho_vac_SCm'],
            'UA': p['rho_vac_UA'],
            'A': p['rho_vac_A'],
            'Ui': p['rho_vac_Ui'],
        },
        'lambda_vac_target': p['lambda_vac_cosmic'],
        'E_0': p['E_0'],
        'n_layers': p['n_layers'],
        'rho_DM_local': p['rho_DM_local'],
        'rho_Lambda': p['rho_Lambda_J'],
        'Lambda': p['Lambda'],
    }


def get_jet_asymmetry_params() -> Dict[str, Any]:
    """
    Get parameters formatted for RelativisticJetAsymmetryModel.
    
    Returns dictionary compatible with:
        model.compute_complete_jet_analysis(**get_jet_asymmetry_params())
    """
    p = RACS_J0320_35_PARAMS
    return {
        'M_bh': p['M_bh_kg'],
        'v_jet': p['v_jet_c'],
        'rho_jet': p['rho_jet_kgm3'],
        'mu': p['mu_viscosity'],
        'L_jet': p['jet_length_m'],
        'beta_i': p['beta_i'],
        'omega': p['omega'],
        't_n1': p['t_n1'],
        't_n2': p['t_n2'],
        'gamma_decay': p['gamma_decay'],
        'E_react': p['E_react_Wm3'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# YANG-MILLS MASS GAP PARAMETERS
# Document: Yang-Mills Mass Gap Proof_20April2025
# Clay Millennium Prize Problem - Existence and Mass Gap proof via UQFF
# ═══════════════════════════════════════════════════════════════════════════════

YANG_MILLS_PARAMS = {
    # UQFF-Yang-Mills coupling
    'lambda_H': 1.0,                        # Higgs coupling constant
    'g_UQFF': 0.1,                          # UQFF coupling = ρ_vac,[SCm]/ρ_vac,[UA]
    'f_quasi': 0.01,                        # Quasi-equilibrium factor
    
    # Vacuum energy densities (J/m³)
    'rho_vac_UA_prime': 1e-23,              # [UA'] vacuum density
    'rho_vac_SCm': 7.09e-37,                # [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,                 # [UA] vacuum density
    'rho_vac_ratio': 0.1,                   # ρ_vac,[SCm]/ρ_vac,[UA]
    
    # Higgs field parameters
    'omega_H': 1.585e-8,                    # s⁻¹ (2π / 3.96×10⁸)
    'omega_c': 2 * np.pi / 3.96e8,          # Characteristic frequency
    'v_H_GeV': 246.0,                       # Higgs VEV (GeV)
    'mu_H_squared': -8.9e4,                 # GeV² (negative for SSB)
    'lambda_quartic': 0.13,                 # Higgs quartic self-coupling
    
    # Pseudo-monopole states
    'phi': 1.0,                             # Phase factor
    'n_states': 26,                         # Quantized states (1-26)
    'SSq': 1.0,                             # [SSq] suppression factor
    
    # Magnetic dipole (for vector potential)
    'mu_j_0': 3.38e23,                      # T·m³ base magnetic moment
    'mu_j_amplitude': 0.4,                  # Oscillation amplitude
    
    # Gauge group parameters (SU(2) or SU(3))
    'gauge_group': 'SU(3)',                 # Default QCD gauge group
    'N_generators_SU2': 3,                  # Lie algebra generators for SU(2)
    'N_generators_SU3': 8,                  # Lie algebra generators for SU(3)
    
    # Physical constants
    'c': 2.998e8,                           # m/s
    'hbar': 1.055e-34,                      # J·s
    'hbar_c_MeV_fm': 197.3,                 # MeV·fm (natural units)
    
    # Mass gap targets (corrected calculation)
    # Note: 69.8 MeV → m = E/c² = 69.8×1.602e-13 / (3e8)² = 1.24e-28 kg
    'm1_kg': 1.24e-28,                      # Mass gap m_1 in kg (corresponds to 70 MeV)
    'E1_J': 1.118e-11,                      # E_1 = m_1 c² in J
    'E1_MeV': 69.8,                         # E_1 in MeV
    
    # Correlation length
    'xi': 1e-15,                            # Correlation length ~1 fm
    
    # Structure constants (for non-abelian)
    'epsilon_abc': 1.0,                     # Levi-Civita normalization
}


def get_yang_mills_params() -> Dict[str, Any]:
    """
    Get parameters formatted for YangMillsMassGapModel.
    
    Returns dictionary compatible with:
        model.compute_complete_analysis(**get_yang_mills_params())
    """
    p = YANG_MILLS_PARAMS
    return {
        'lambda_H': p['lambda_H'],
        'rho_vac_UA_prime': p['rho_vac_UA_prime'],
        'rho_vac_SCm': p['rho_vac_SCm'],
        'rho_vac_UA': p['rho_vac_UA'],
        'omega_H': p['omega_H'],
        'f_quasi': p['f_quasi'],
        'SSq': p['SSq'],
        'n_states': p['n_states'],
        'phi': p['phi'],
        'gauge_group': p['gauge_group'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# RIEMANN HYPOTHESIS PARAMETERS
# Document: Riemann Hypothesis_20April2025
# Clay Millennium Prize Problem - Critical Line proof via UQFF
# ═══════════════════════════════════════════════════════════════════════════════

RIEMANN_HYPOTHESIS_PARAMS = {
    # Pseudo-monopole quantization
    'phi': 1.0,                             # Normalized Higgs field strength
    'n_states': 26,                         # Number of UQFF quantum states
    'SSq': 1.0,                             # [SSq] suppression factor
    
    # Vacuum energy densities (J/m³)
    'rho_vac_UA_prime': 1e-23,              # [UA'] vacuum density
    'rho_vac_SCm': 7.09e-37,                # [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,                 # [UA] vacuum density
    'vac_ratio': 0.1,                       # ρ_vac,[SCm]/ρ_vac,[UA]
    
    # Universal Magnetism parameters
    'mu_j_0': 3.38e20,                      # T·m³ base magnetic moment
    'mu_j_amplitude': 0.4,                  # Oscillation amplitude
    'omega_c': 2 * np.pi / 3.96e8,          # Characteristic frequency
    'gamma': 0.0005,                        # Decay rate
    
    # Known zeta zeros (for validation)
    't_1': 14.134725,                       # First non-trivial zero imaginary part
    't_2': 21.022040,                       # Second zero
    't_3': 25.010858,                       # Third zero
    't_4': 30.424876,                       # Fourth zero
    't_5': 32.935062,                       # Fifth zero
    
    # Critical line
    'sigma_critical': 0.5,                  # σ = 1/2 (Riemann Hypothesis)
    
    # UQFF-zeta function analog terms
    'f_quasi': 0.01,                        # Quasi-equilibrium factor
    'pi': np.pi,                            # π constant
    
    # Resonance conditions
    'delta_6': 2 * np.pi,                   # δ_6 = 2π (full rotation)
    'omega_n_base': 1.0,                    # Base frequency for ω_n
    
    # Physical constants
    'hbar': 1.055e-34,                      # J·s
    'c': 2.998e8,                           # m/s
}


def get_riemann_params() -> Dict[str, Any]:
    """
    Get parameters formatted for RiemannHypothesisModel.
    
    Returns dictionary compatible with:
        model.compute_complete_analysis(**get_riemann_params())
    """
    p = RIEMANN_HYPOTHESIS_PARAMS
    return {
        'phi': p['phi'],
        'n_states': p['n_states'],
        'SSq': p['SSq'],
        'rho_vac_UA_prime': p['rho_vac_UA_prime'],
        'rho_vac_SCm': p['rho_vac_SCm'],
        'rho_vac_UA': p['rho_vac_UA'],
        'sigma_critical': p['sigma_critical'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# P VS NP PARAMETERS
# Document: P vs. NP Proof_20April2025
# Clay Millennium Prize Problem - Computational complexity via UQFF non-locality
# ═══════════════════════════════════════════════════════════════════════════════

P_VS_NP_PARAMS = {
    # Universal Magnetism components
    'mu_j_0': 3.38e20,                      # T·m³ base magnetic moment
    'mu_j_amplitude': 0.4,                  # Oscillation amplitude
    'mu_j_base': 1e3,                       # Base oscillation (10³)
    'omega_c': 2 * np.pi / 3.96e8,          # Characteristic frequency
    'gamma': 0.0005,                        # Decay rate κ
    
    # Vacuum energy densities (J/m³)
    'rho_vac_UA_prime': 1e-23,              # [UA'] vacuum density
    'rho_vac_SCm': 7.09e-37,                # [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,                 # [UA] vacuum density
    
    # UQFF quantum states
    'n_states': 26,                         # Number of quantum states
    'SSq': 1.0,                             # [SSq] suppression factor (normalized)
    
    # Reactivity and quasi-equilibrium
    'E_react_0': 1e46,                      # E_react(0) base energy
    'kappa': 0.0005,                        # Decay rate /day
    'f_Heaviside': 0.01,                    # Heaviside step factor
    'f_quasi': 0.01,                        # Quasi-equilibrium factor
    
    # SCm pressure term
    'P_SCm': 1.0,                           # Normalized SCm pressure
    
    # Computational complexity parameters
    'm_exponent': 2,                        # Polynomial exponent for P problems
    'k_input_size': 100,                    # Reference input size
    
    # Physical constants
    'pi': np.pi,
    'e': np.e,
    
    # Test problems (SAT reference)
    'sat_verify_exponent': 2,               # T_verify ∝ k²
    'sat_solve_base': 'exponential',        # T_solve ∝ e^([SSq]n/26 · k)
}


def get_p_vs_np_params() -> Dict[str, Any]:
    """
    Get parameters formatted for PvsNPComplexityModel.
    
    Returns dictionary compatible with:
        model.compute_complete_analysis(**get_p_vs_np_params())
    """
    p = P_VS_NP_PARAMS
    return {
        'mu_j_0': p['mu_j_0'],
        'mu_j_amplitude': p['mu_j_amplitude'],
        'mu_j_base': p['mu_j_base'],
        'omega_c': p['omega_c'],
        'gamma': p['gamma'],
        'rho_vac_UA_prime': p['rho_vac_UA_prime'],
        'rho_vac_SCm': p['rho_vac_SCm'],
        'E_react_0': p['E_react_0'],
        'kappa': p['kappa'],
        'f_Heaviside': p['f_Heaviside'],
        'f_quasi': p['f_quasi'],
        'n_states': p['n_states'],
        'SSq': p['SSq'],
        'm_exponent': p['m_exponent'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# DUST YIELD AND EXTINCTION PARAMETERS
# 71-Equation Catalog: A_V extinction and y_dust yield equations
# References: Milky Way dust models, star formation regions
# ═══════════════════════════════════════════════════════════════════════════════

DUST_YIELD_PARAMS = {
    # Dust-to-gas ratios
    'M_dust': 1e-2,                         # M_☉ - dust mass (typical for 1 M_☉ gas)
    'M_gas': 1.0,                           # M_☉ - gas mass
    'dust_to_gas_ratio': 0.01,              # dimensionless - MW average
    
    # Dust opacity
    'kappa_dust': 1e4,                      # m²/kg - dust opacity at optical
    'kappa_V': 1e4,                         # m²/kg - V-band opacity
    
    # Metallicity
    'Z': 0.02,                              # Solar metallicity
    'Z_solar': 0.02,                        # Reference solar metallicity
    
    # Timescales
    'tau': 1e7,                             # Gyr → years (10 Myr typical)
    'tau_SF': 1e6,                          # years - star formation timescale
    'tau_depletion': 2e9,                   # years - gas depletion time
    
    # UQFF fundamental parameters
    'nu_fund': 0.618,                       # Golden ratio exponent
    'kappa': 0.0005,                        # 1/day - decay constant
    'SSq': 0.57,                            # [SSq] calibrated value
    
    # Extinction coefficients
    'R_V': 3.1,                             # Total-to-selective extinction ratio
    'A_V_MW_mean': 1.0,                     # mag - typical MW extinction
    
    # Dust yield targets
    'y_dust_SN': 0.1,                       # M_☉ - SN dust yield per event
    'y_dust_AGB': 0.01,                     # M_☉ - AGB dust yield
    
    # Physical constants
    'mag_factor': 1.086,                    # -2.5 / ln(10) ≈ 1.086
}


def get_dust_yield_params() -> Dict[str, Any]:
    """
    Get parameters formatted for DustYieldExtinctionModel.
    
    Returns dictionary compatible with:
        model.compute_complete_analysis(**get_dust_yield_params())
    """
    p = DUST_YIELD_PARAMS
    return {
        'M_dust': p['M_dust'],
        'M_gas': p['M_gas'],
        'kappa_dust': p['kappa_dust'],
        'Z': p['Z'],
        'tau': p['tau'],
        'tau_SF': p['tau_SF'],
        'nu_fund': p['nu_fund'],
        'R_V': p['R_V'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SAGITTARIUS A* GRAVITY PARAMETERS
# 71-Equation Catalog: g_SgrA*(r,t) with M(t), B(t) time-dependence
# References: GRAVITY Collaboration, EHT 2022, ALMA observations
# ═══════════════════════════════════════════════════════════════════════════════

SGRA_STAR_PARAMS = {
    # Sagittarius A* properties
    'M_SgrA': 4.0e6,                        # M_☉ - SMBH mass
    'M_SgrA_kg': 4.0e6 * SOLAR_MASS_KG,     # kg
    'R_s': 1.2e10,                          # m - Schwarzschild radius (2GM/c²)
    'distance_pc': 8000,                    # pc - distance to galactic center
    'distance_m': 8000 * 3.086e16,          # m
    
    # Accretion parameters
    'M_dot_Edd': 1e-3,                      # M_☉/yr - Eddington accretion rate
    'M_dot_SgrA': 1e-8,                     # M_☉/yr - actual accretion rate (~10⁻⁵ Eddington)
    'M_dot_SgrA_kg_s': 1e-8 * SOLAR_MASS_KG / (365.25 * 24 * 3600),  # kg/s
    'accretion_efficiency': 0.1,            # η - radiative efficiency
    
    # Magnetic field parameters
    'B_horizon': 1e2,                       # T - field at event horizon
    'B_decay_rate': 1e-10,                  # 1/s - slow decay timescale
    'mu_0': 1e20,                           # T·m³ - magnetic moment at t=0
    
    # UQFF Ug terms
    'Ug1_amplitude': 1e-8,                  # m/s² - magnetic dipole contribution
    'Ug2_amplitude': 1e-12,                 # m/s² - charge-reactivity
    'Ug3_amplitude': 1e-10,                 # m/s² - string rotation
    'Ug4_amplitude': 1e-15,                 # m/s² - vacuum concentration
    
    # Cosmological correction
    'Lambda': 1.1e-52,                      # m⁻² - cosmological constant
    
    # UQFF modulation
    'kappa': 0.0005,                        # 1/day - decay constant
    'SSq': 0.57,                            # [SSq] calibrated
    'beta_i': 0.603,                        # Buoyancy parameter
    
    # Physical constants
    'G': G,
    'c': SPEED_OF_LIGHT,
}


def get_sgra_star_params() -> Dict[str, Any]:
    """
    Get parameters formatted for SgrAStarGravityModel.
    
    Returns dictionary compatible with:
        model.compute_g_SgrA(**get_sgra_star_params())
    """
    p = SGRA_STAR_PARAMS
    return {
        'M_0': p['M_SgrA_kg'],
        'r': p['R_s'] * 10,                 # 10 R_s reference
        't': 0.0,
        'M_dot': p['M_dot_SgrA_kg_s'],
        'B_0': p['B_horizon'],
        'tau_B': 1.0 / p['B_decay_rate'] if p['B_decay_rate'] > 0 else 1e15,
        'Lambda': p['Lambda'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SHEAR CHI-SQUARED PARAMETERS
# 71-Equation Catalog: χ² = Σ(P_obs - P_ucf(δ_τ))²/σ_P² for shear maps
# References: JWST WL surveys, DES Y3, KiDS-1000
# ═══════════════════════════════════════════════════════════════════════════════

SHEAR_CHI_SQUARED_PARAMS = {
    # JWST survey parameters
    'survey_name': 'JWST_WL',
    'survey_area_deg2': 0.6,                # deg² - typical JWST pointing
    'n_galaxies': 1000,                     # Number of background galaxies
    
    # Shear calibration
    'delta_tau_best': 0.05,                 # Optimal time calibration offset
    'delta_tau_range_min': 0.01,            # Minimum δ_τ
    'delta_tau_range_max': 0.10,            # Maximum δ_τ
    'delta_tau_search_points': 50,          # Grid search resolution
    
    # Observed shear power spectrum (mock data P_ell for ℓ = 100,200,...1000)
    'ell_values': [100, 200, 300, 500, 700, 1000],
    'P_obs_values': [1.2e-4, 0.9e-4, 0.7e-4, 0.5e-4, 0.35e-4, 0.22e-4],
    'sigma_P_values': [0.1e-4, 0.08e-4, 0.07e-4, 0.06e-4, 0.05e-4, 0.04e-4],
    
    # UQFF model parameters
    'A_amp': 1e-3,                          # Shear amplitude normalization
    'n_spectral': -2.0,                     # Spectral index
    'z_source': 1.0,                        # Source redshift
    
    # Chi-squared thresholds
    'chi2_reduced_good': 1.5,               # χ²_red < 1.5 is good fit
    'chi2_reduced_acceptable': 2.0,         # χ²_red < 2.0 is acceptable
    'dof_default': 5,                       # Degrees of freedom (N_points - 1)
    
    # UQFF calibration
    'kappa': 0.0005,                        # 1/day
    'SSq': 0.57,                            # [SSq]
}


def get_shear_chi_squared_params() -> Dict[str, Any]:
    """
    Get parameters formatted for ShearChiSquaredModel.
    
    Returns dictionary compatible with:
        model.compute_chi_squared(**get_shear_chi_squared_params())
    """
    p = SHEAR_CHI_SQUARED_PARAMS
    return {
        'P_obs': p['P_obs_values'],
        'sigma_P': p['sigma_P_values'],
        'ell_values': p['ell_values'],
        'delta_tau': p['delta_tau_best'],
        'A_amp': p['A_amp'],
        'n_spectral': p['n_spectral'],
        'z': p['z_source'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# TIME ASYMMETRY PARAMETERS
# 71-Equation Catalog: ∂t ≠ 0 irreversibility (dE/dt < 0, dS/dt > 0)
# References: UQFF Second Law, arrow of time framework
# ═══════════════════════════════════════════════════════════════════════════════

TIME_ASYMMETRY_PARAMS = {
    # Reactivity energy
    'E_react_0': 1e46,                      # J - initial reactivity energy
    'E_react_unit': 'J',
    
    # Decay parameters
    'kappa': 0.0005,                        # 1/day - UQFF calibrated decay constant
    'kappa_unit': '1/day',
    
    # Entropy production
    'S_0': 1e80,                            # J/K - cosmic entropy scale
    'T_ref': 2.725,                         # K - CMB temperature
    
    # Quasi-equilibrium evolution
    'f_quasi_0': 0.01,                      # Initial quasi-equilibrium factor
    'f_quasi_final': 1.0,                   # Final value as t → ∞
    
    # Time scales for verification
    't_test_min': 1.0,                      # days - minimum test time
    't_test_max': 1e5,                      # days - maximum test time (avoiding underflow)
    't_test_points': 50,                    # Number of test points
    
    # Arrow of time conditions
    'dE_dt_sign': 'negative',               # dE/dt < 0 (dissipation)
    'dS_dt_sign': 'positive',               # dS/dt > 0 (entropy production)
    
    # UQFF parameters
    'SSq': 0.57,                            # [SSq] calibrated
    'beta_i': 0.603,                        # Buoyancy opposition
    
    # Physical interpretation
    'description': 'Time irreversibility encoded in UQFF via κ decay and f_quasi → 1',
}


def get_time_asymmetry_params() -> Dict[str, Any]:
    """
    Get parameters formatted for TimeAsymmetryModel.
    
    Returns dictionary compatible with:
        model.verify_time_asymmetry(**get_time_asymmetry_params())
    """
    p = TIME_ASYMMETRY_PARAMS
    return {
        'E_react_0': p['E_react_0'],
        'kappa': p['kappa'],
        'S_0': p['S_0'],
        'T_ref': p['T_ref'],
        'f_quasi_0': p['f_quasi_0'],
        't_min': p['t_test_min'],
        't_max': p['t_test_max'],
        'n_points': p['t_test_points'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# NUCLEAR BINDING SHELL LEVELS PARAMETERS
# Document: UQFF proof set for Nuclear Binding Shell Levels_28Sept2025.docx
# PDG 2025 verification of E_n = E_0 × 10^n (n=8 nuclear binding, n=12 Higgs)
# References:
#   - PDG 2025: https://pdg.lbl.gov/2025/reviews/rpp2024-rev-passage-particles-matter.pdf
#   - ResearchGate: Shell Model Calculations 133/135Sn, 133/135Sb
#   - ENSDF: https://www.nndc.bnl.gov/ensdf/
# ═══════════════════════════════════════════════════════════════════════════════

NUCLEAR_BINDING_SHELL_PARAMS = {
    # ─────────────────────────────────────────────────────────────────────────
    # UQFF 26-LEVEL POLYNOMIAL BASE
    # E_n = E_0 × 10^n where E_0 = 10^{-20} J
    # ─────────────────────────────────────────────────────────────────────────
    'E_0': 1e-20,                           # J - vacuum fluctuation minimum
    'E_0_unit': 'J',
    'n_levels': 26,                         # Total levels in hierarchy
    'n_nuclear': 8,                         # Level for nuclear binding (~10^-12 J)
    'n_proton': 10,                         # Level for proton mass (~10^-10 J)
    'n_higgs': 12,                          # Level for Higgs boson (~10^-8 J)
    
    # ─────────────────────────────────────────────────────────────────────────
    # PDG 2025 REFERENCE VALUES
    # ─────────────────────────────────────────────────────────────────────────
    'B_per_nucleon_MeV': 8.0,               # MeV - average binding per nucleon
    'B_per_nucleon_J': 1.28e-12,            # J - converted
    'm_Higgs_GeV': 125.18,                  # GeV - PDG 2025
    'm_Higgs_J': 2.005e-8,                  # J - E_H = m_H c²
    
    # ─────────────────────────────────────────────────────────────────────────
    # SEMI-EMPIRICAL MASS FORMULA (BETHE-WEIZSÄCKER)
    # B(A,Z) = a_V×A - a_S×A^(2/3) - a_C×Z²/A^(1/3) - a_A×(N-Z)²/A + δ
    # ─────────────────────────────────────────────────────────────────────────
    'a_V': 15.8,                            # MeV - Volume term coefficient
    'a_S': 18.3,                            # MeV - Surface term coefficient
    'a_C': 0.714,                           # MeV - Coulomb term coefficient
    'a_A': 23.2,                            # MeV - Asymmetry term coefficient
    'a_P': 12.0,                            # MeV - Pairing term coefficient
    
    # ─────────────────────────────────────────────────────────────────────────
    # MAGIC NUMBERS (SHELL CLOSURES)
    # ─────────────────────────────────────────────────────────────────────────
    'magic_numbers': [2, 8, 20, 28, 50, 82, 126],  # Nuclear magic numbers
    
    # ─────────────────────────────────────────────────────────────────────────
    # Pb-206 ENSDF REFERENCE DATA
    # Source: NNDC ENSDF 2025, Nuclear Data Sheets 201, 346
    # ─────────────────────────────────────────────────────────────────────────
    'Pb206_A': 206,                         # Mass number
    'Pb206_Z': 82,                          # Atomic number (magic)
    'Pb206_N': 124,                         # Neutron number (near magic 126)
    'Pb206_B_per_A': 7.88,                  # MeV/nucleon (ENSDF)
    'Pb206_levels_MeV': [                   # Sample ENSDF levels (17 points)
        0.0, 0.044, 0.137, 0.334, 0.583, 0.802, 1.028,
        1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0
    ],
    
    # ─────────────────────────────────────────────────────────────────────────
    # Sn/Sb ISOTOPE REFERENCE DATA (from ResearchGate paper)
    # Shell Model Calculations of Energy Levels and Binding Energy of
    # 133/135Sn and 133/135Sb Nuclei
    # ─────────────────────────────────────────────────────────────────────────
    'Sn133_A': 133,
    'Sn133_Z': 50,                          # Magic Z=50
    'Sn133_N': 83,                          # Near magic N=82
    'Sn135_A': 135,
    'Sn135_Z': 50,
    'Sn135_N': 85,
    'Sb133_A': 133,
    'Sb133_Z': 51,                          # One proton above magic
    'Sb133_N': 82,                          # Magic N=82
    'Sb135_A': 135,
    'Sb135_Z': 51,
    'Sb135_N': 84,
    
    # ─────────────────────────────────────────────────────────────────────────
    # POLYNOMIAL FIT PARAMETERS
    # V(r) ≈ Σ_{n=1}^{26} a_n r^n
    # R² ≈ 0.95 for low degrees, ~1 for deg=26 (overfit)
    # ─────────────────────────────────────────────────────────────────────────
    'poly_deg_low': 5,                      # Low degree for physical fit
    'poly_deg_high': 16,                    # Max for 17 data points
    'R_squared_target': 0.95,               # Physical fit quality target
    
    # ─────────────────────────────────────────────────────────────────────────
    # CONVERSION FACTORS
    # ─────────────────────────────────────────────────────────────────────────
    'MeV_to_J': 1.602e-13,                  # MeV → J
    'GeV_to_J': 1.602e-10,                  # GeV → J
    'keV_to_J': 1.602e-16,                  # keV → J
    'fm_to_m': 1e-15,                       # fm → m
    
    # ─────────────────────────────────────────────────────────────────────────
    # QCD CORNELL POTENTIAL PARAMETERS
    # V(r) = σr - α_s/r (confinement + Coulomb)
    # ─────────────────────────────────────────────────────────────────────────
    'sigma_GeV2': 0.18,                     # GeV² - string tension
    'alpha_s': 0.39,                        # Strong coupling at charm scale
    'GeV_fm': 0.197,                        # ℏc in GeV·fm (natural units)
    
    # ─────────────────────────────────────────────────────────────────────────
    # LOW-n LEVEL PARAMETERS (n=1-5)
    # ─────────────────────────────────────────────────────────────────────────
    'n_low_max': 5,                         # Maximum low-n level
    'E_4_keV': 0.624,                       # E_4 ≈ 1 keV (LHC virtual quarks)
    'ATLAS_CONF': 'ATLAS-CONF-2025-007',    # LHC verification source
    
    # Description
    'description': 'PDG 2025 verified UQFF 26-level polynomial structure',
    'document': 'UQFF proof set for Nuclear Binding Shell Levels_28Sept2025.docx',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
}


def get_nuclear_binding_shell_params() -> Dict[str, Any]:
    """
    Get parameters formatted for NuclearBindingShellLevelsModel.
    
    Returns dictionary compatible with:
        model = NuclearBindingShellLevelsModel()
        model.compute_semi_empirical_binding(A, Z)
        model.polynomial_fit_shell_levels(deg)
    """
    p = NUCLEAR_BINDING_SHELL_PARAMS
    return {
        'E_0': p['E_0'],
        'B_per_nucleon_MeV': p['B_per_nucleon_MeV'],
        'm_Higgs_GeV': p['m_Higgs_GeV'],
        'bethe_weiszacker': {
            'a_V': p['a_V'],
            'a_S': p['a_S'],
            'a_C': p['a_C'],
            'a_A': p['a_A'],
            'a_P': p['a_P'],
        },
        'magic_numbers': p['magic_numbers'],
        'Pb206_levels': p['Pb206_levels_MeV'],
        'poly_deg': p['poly_deg_low'],
        'R_squared_target': p['R_squared_target'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UQFF MASTER FRAMEWORK PARAMETERS
# Document: UQFF Compressed_Resonant_Buoyancy_Superconductive_Triadic_Quadratic_
#           Master Buoyancy_29Sept2025.docx
# Complete 7-mode UQFF operational framework
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_MASTER_FRAMEWORK_PARAMS = {
    # ─────────────────────────────────────────────────────────────────────────
    # COUPLING CONSTANTS (k_i for Ug_i gravity terms)
    # ─────────────────────────────────────────────────────────────────────────
    'k_1': 1.5,                             # Ug1 (magnetic dipole) coupling
    'k_2': 1.2,                             # Ug2 (bubble) coupling
    'k_3': 1.8,                             # Ug3 (disk) coupling
    'k_4': 1.0,                             # Ug4 (BH) coupling
    
    # ─────────────────────────────────────────────────────────────────────────
    # BUOYANCY PARAMETERS (Mode 3)
    # ─────────────────────────────────────────────────────────────────────────
    'beta_i': 0.61,                         # Buoyancy opposition scale
    'delta_sw': 0.01,                       # Solar wind modulation
    'lambda_vac_sw': 7.2e-4,                # J/m³ - solar wind vacuum density
    'UA_charge': 1e-11,                     # C - Universal Aether charge
    
    # ─────────────────────────────────────────────────────────────────────────
    # GALACTIC PARAMETERS
    # ─────────────────────────────────────────────────────────────────────────
    'omega_g': 7.3e-16,                     # rad/s - galactic spin (220 km/s ÷ 8 kpc)
    'M_bh': 8.15e36,                        # kg - Sgr A* mass (4.1×10⁶ M_☉)
    'd_g': 2.55e20,                         # m - galactic distance (27,000 ly)
    
    # ─────────────────────────────────────────────────────────────────────────
    # SUPERCONDUCTIVE PARAMETERS (Mode 4)
    # E_react = (ρ_vac,[SCm] × v_SCm²) / ρ_A × e^(-κt) = 10^46 × e^(-0.0005t)
    # ─────────────────────────────────────────────────────────────────────────
    'rho_vac_SCm': 7.09e-37,                # J/m³ - [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,                 # J/m³ - [UA] vacuum density
    'rho_A': 1e-23,                         # J/m³ - Aether vacuum density
    'v_SCm': 1e8,                           # m/s - SCm velocity (c/3)
    'kappa': 0.0005,                        # 1/day - reactivity decay constant
    'E_react_base': 1e46,                   # Base reactivity energy
    
    # ─────────────────────────────────────────────────────────────────────────
    # RESONANT PARAMETERS (Mode 2)
    # cos(πt_n) oscillation with period 2 days
    # ─────────────────────────────────────────────────────────────────────────
    'gamma': 5e-5,                          # 1/day - decay rate (τ ~55 yr)
    't_0': 0.0,                             # Reference time for t_n = t - t_0
    'f_TRZ': 0.1,                           # TRZ factor for negentropy
    
    # ─────────────────────────────────────────────────────────────────────────
    # STRING PARAMETERS (μ_j contribution)
    # ─────────────────────────────────────────────────────────────────────────
    'mu_j_base': 3.38e20,                   # T·m³ - base magnetic moment
    'mu_j_oscillation': 1e3,                # Base oscillation amplitude (10³)
    'mu_j_amplitude': 0.4,                  # Oscillation amplitude factor
    'r_j': 1.496e13,                        # m - string distance (AU scale)
    'phi_j': 1.0,                           # Unit vector (dimensionless)
    
    # ─────────────────────────────────────────────────────────────────────────
    # METRIC AND STRESS-ENERGY (Mode 1 contribution)
    # ─────────────────────────────────────────────────────────────────────────
    'g_munu_trace': -2,                     # Minkowski metric trace [1,-1,-1,-1]
    'eta_coupling': 1e-22,                  # Metric coupling constant
    'T_s_munu': 1.123e7,                    # J/m³ - stress-energy (1.27e3 + 1.11e7)
    
    # ─────────────────────────────────────────────────────────────────────────
    # INERTIA PARAMETERS
    # ─────────────────────────────────────────────────────────────────────────
    'delta_i': 1.0,                         # Inertial coupling
    'lambda_i': 1.0,                        # Inertia scaling factor
    
    # ─────────────────────────────────────────────────────────────────────────
    # CRP (COSMIC RAY PROPAGATION) PARAMETERS
    # n = p^(-2.2) × exp(-p/p_max), D_E = D_0 × E^0.5
    # ─────────────────────────────────────────────────────────────────────────
    'D_0': 1e28,                            # m²/s - diffusion coefficient base
    'gamma_CRP': 2.2,                       # Spectral index for n(p)
    'p_max': 1e6,                           # GeV/c - maximum momentum cutoff
    
    # ─────────────────────────────────────────────────────────────────────────
    # TRIADIC PARAMETERS (Mode 5)
    # F_U_tri = (Ug3 × Ub_i × Um)^(1/3) × exp(-[SSq]n/26)
    # ─────────────────────────────────────────────────────────────────────────
    'n_triadic': 13,                        # Triadic plasma layer
    'SSq_triadic': 38,                      # Self-similar quotient (log ratio ~10^-38)
    
    # ─────────────────────────────────────────────────────────────────────────
    # QUADRATIC APPROXIMATION (Mode 6)
    # V(r) ≈ a_0 + a_1 × r + a_2 × r²
    # ─────────────────────────────────────────────────────────────────────────
    'a_0': 0.0,                             # Constant term
    'a_1': 0.0,                             # Linear term
    'a_2': 1e-12,                           # Quadratic term (for n=8 bindings)
    'R_squared': 0.95,                      # Fit quality for low-degree approx
    
    # ─────────────────────────────────────────────────────────────────────────
    # MASTER BUOYANCY (Mode 7)
    # Master Ub_i = Ub_i + exp(-(π - t)) × Um / ρ_vac,[UA]
    # ─────────────────────────────────────────────────────────────────────────
    'Um_reference': 2.28e65,                # J/m³ - magnetism reference (Sun)
    'mayan_alignment': True,                # Use exp(-(π - t)) term
    
    # ─────────────────────────────────────────────────────────────────────────
    # 7 OPERATIONAL MODES
    # ─────────────────────────────────────────────────────────────────────────
    'modes': {
        1: 'UQFF_Compressed',               # Compact F_U form
        2: 'UQFF_Resonant',                 # cos(πt_n) oscillations
        3: 'UQFF_Buoyancy',                 # Ub_i opposition
        4: 'UQFF_Superconductive',          # E_react [SCm] modulation
        5: 'UQFF_Triadic',                  # Geometric mean
        6: 'UQFF_Quadratic',                # V(r) approximation
        7: 'UQFF_Master_Buoyancy',          # Extended Ub_i
    },
    
    # Document reference
    'document': 'UQFF Compressed_Resonant_Buoyancy_Superconductive_Triadic_Quadratic_Master Buoyancy_29Sept2025.docx',
}


def get_uqff_master_framework_params() -> Dict[str, Any]:
    """
    Get parameters formatted for UQFFMasterFramework.
    
    Returns dictionary compatible with:
        framework = UQFFMasterFramework()
        result = framework.compute_complete_F_U(**get_uqff_master_framework_params())
    """
    p = UQFF_MASTER_FRAMEWORK_PARAMS
    return {
        'k': [p['k_1'], p['k_2'], p['k_3'], p['k_4']],
        'beta_i': p['beta_i'],
        'omega_g': p['omega_g'],
        'M_bh': p['M_bh'],
        'd_g': p['d_g'],
        'kappa': p['kappa'],
        'gamma': p['gamma'],
        'rho_vac_SCm': p['rho_vac_SCm'],
        'rho_vac_UA': p['rho_vac_UA'],
        'v_SCm': p['v_SCm'],
        'Um_reference': p['Um_reference'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# SOLAR WIND PARKER PROBE PARAMETERS
# Document: UQFF proof set for Solar Wind Density for Parker Solar Probe (CDAWeb 2025)_28Sept2025.docx
# Verification: δ_sw=0.01, v_sw=5e5 m/s, ρ_sw~8e-21 kg/m³
# References:
#   - NASA CDAWeb: https://cdaweb.gsfc.nasa.gov/
#   - Parker Solar Probe SWEAP: https://sweap.cfa.harvard.edu/
#   - Grok conversation: https://x.com/i/grok?conversation=1972174124559557054
# ═══════════════════════════════════════════════════════════════════════════════

SOLAR_WIND_PARKER_PROBE_PARAMS = {
    # ─────────────────────────────────────────────────────────────────────────
    # UQFF SOLAR WIND MODULATION PARAMETERS
    # Term: (1 + δ_sw × v_sw) in Ug2
    # ─────────────────────────────────────────────────────────────────────────
    'delta_sw': 0.01,                       # Solar wind modulation factor (unitless)
    'delta_sw_unit': 'dimensionless',
    
    'v_sw': 5e5,                            # m/s - model solar wind velocity (500 km/s)
    'v_sw_km_s': 500,                       # km/s - same in conventional units
    'v_sw_unit': 'm/s',
    
    # ─────────────────────────────────────────────────────────────────────────
    # SOLAR WIND MODULATION FACTOR
    # 1 + δ_sw × v_sw = 1 + 0.01 × 5×10⁵ = 5001
    # ─────────────────────────────────────────────────────────────────────────
    'sw_factor': 5001,                      # 1 + δ_sw × v_sw = 5001× enhancement
    'sw_factor_formula': '1 + 0.01 × 5e5 = 5001',
    
    # ─────────────────────────────────────────────────────────────────────────
    # PARKER SOLAR PROBE CDAWeb 2025 OBSERVED VALUES
    # Source: SWEAP instrument, Encounters 20-25 at 1 AU
    # ─────────────────────────────────────────────────────────────────────────
    'n_p_min': 4e6,                         # m⁻³ - minimum proton density (4 cm⁻³)
    'n_p_max': 10e6,                        # m⁻³ - maximum proton density (10 cm⁻³)
    'n_p_avg': 7e6,                         # m⁻³ - average proton density (7 cm⁻³)
    'n_p_typical': 5e6,                     # m⁻³ - typical density for ρ_sw~8e-21
    
    'v_sw_slow': 3e5,                       # m/s - slow solar wind (300 km/s)
    'v_sw_fast': 8e5,                       # m/s - fast solar wind (800 km/s)
    'v_sw_observed_avg': 5.5e5,             # m/s - observed average (550 km/s)
    
    # ─────────────────────────────────────────────────────────────────────────
    # SOLAR WIND MASS DENSITY
    # ρ_sw = m_p × n_p
    # ─────────────────────────────────────────────────────────────────────────
    'm_p': 1.672e-27,                       # kg - proton mass
    'rho_sw_expected': 8e-21,               # kg/m³ - expected solar wind density
    'rho_sw_min': 6.7e-21,                  # kg/m³ - at n_p=4 cm⁻³
    'rho_sw_max': 1.67e-20,                 # kg/m³ - at n_p=10 cm⁻³
    'rho_sw_formula': 'ρ_sw = m_p × n_p = 1.672e-27 × 5e6 = 8.36e-21 kg/m³',
    
    # ─────────────────────────────────────────────────────────────────────────
    # HELIOSPHERE PARAMETERS
    # ─────────────────────────────────────────────────────────────────────────
    'R_b': 1.496e13,                        # m - 100 AU boundary for Ug2
    'heliopause': 1.824e13,                 # m - ~122 AU heliopause
    'termination_shock': 1.2e13,            # m - ~80-100 AU
    
    # ─────────────────────────────────────────────────────────────────────────
    # VERIFICATION METRICS
    # ─────────────────────────────────────────────────────────────────────────
    'v_error_percent_max': 20.0,            # Acceptable velocity error (%)
    'rho_error_percent_max': 50.0,          # Acceptable density error (%)
    
    # Description and references
    'description': 'Parker Solar Probe CDAWeb 2025 solar wind verification',
    'document': 'UQFF proof set for Solar Wind Density for Parker Solar Probe (CDAWeb 2025)_28Sept2025.docx',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
    'data_source': 'NASA CDAWeb, Parker Solar Probe SWEAP',
}


# ═══════════════════════════════════════════════════════════════════════════════
# ALPHA BEC LENR PARAMETERS
# Tohsaki et al. AMD verification of Bose term N_B and T_c shifts
# ═══════════════════════════════════════════════════════════════════════════════

ALPHA_BEC_LENR_PARAMS = {
    # Physical constants
    'hbar': 1.055e-34,                      # J·s (reduced Planck constant)
    'k_B': 1.381e-23,                       # J/K (Boltzmann constant)
    'm_proton': 1.673e-27,                  # kg (proton mass)
    'm_alpha': 6.646e-27,                   # kg (alpha particle mass = 4 × m_p)
    
    # Alpha particle properties
    'alpha_binding_MeV': 28.3,              # MeV (alpha binding energy)
    'alpha_radius_fm': 1.67,                # fm (RMS charge radius)
    
    # Nuclear densities
    'rho_0_fm3': 0.17,                      # fm⁻³ (nuclear saturation density)
    'rho_BEC_fm3': 0.03,                    # fm⁻³ (dilute limit for alpha BEC)
    
    # BEC parameters
    'zeta_3_2': 2.612,                      # ζ(3/2) Riemann zeta function
    
    # ¹²C Hoyle state (7.65 MeV)
    'E_Hoyle_MeV': 7.65,                    # MeV above ground state
    'N_B_C12': 3,                           # 3-alpha cluster
    'condensate_fraction_Hoyle': 0.70,      # 70% in 0S state (arXiv:1103.3940)
    'Hoyle_lifetime_s': 2.4e-16,            # s (Hoyle state lifetime)
    
    # ¹⁶O alpha states
    'E_O16_alpha_MeV': 14.4,                # MeV (0⁺ state)
    'N_B_O16': 4,                           # 4-alpha cluster
    'condensate_fraction_O16': 0.60,        # 60% condensate
    
    # ⁸Be (unstable 2-alpha)
    'N_B_Be8': 2,                           # 2-alpha cluster
    'condensate_fraction_Be8': 0.95,        # 95% condensate (ground state)
    'Be8_lifetime_s': 1e-16,                # s (extremely short)
    
    # ²⁰Ne 5-alpha
    'N_B_Ne20': 5,                          # 5-alpha cluster
    'condensate_fraction_Ne20': 0.50,       # 50% condensate
    
    # LENR parameters
    'E_react_Wm3': 1e46,                    # W/m³ (UQFF vacuum reactor energy density)
    'SCm_UA_ratio': 1e-38,                  # [SCm]/[UA] typical ratio
    'E_nuclear_J': 1e-12,                   # J (nuclear binding energy scale)
    
    # T_c calculation parameters
    'T_c_base_K': 1.2e6,                    # K (base critical temperature)
    'delta_T_c_LENR_K': 300,                # K (LENR temperature shift)
    
    # AMD wave function parameters (Tohsaki)
    'b_gaussian_fm': 1.52,                  # fm (Gaussian width parameter)
    
    # Conversion factors
    'MeV_to_J': 1.602e-13,
    'fm_to_m': 1e-15,
    
    # Verification targets
    'N_B_AMD_C12': 3,                       # AMD reference for ¹²C
    'N_B_AMD_O16': 4,                       # AMD reference for ¹⁶O
    'condensate_fraction_AMD_C12': 0.70,    # AMD reference 70% ± 10%
    'condensate_fraction_AMD_error': 0.10,  # ±10% uncertainty
    
    # Description and references
    'description': 'Alpha BEC LENR verification via Tohsaki et al. AMD',
    'document': 'UQFF proof set verification of Bose term N_B, T_c shifts for alpha BEC_29Sept2025.docx',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
    'arxiv_reference': 'arXiv:1103.3940',
    'inis_reference': 'https://inis.iaea.org/records/3164a-q0271',
}


# ═══════════════════════════════════════════════════════════════════════════════
# 26D COSMIC EGG HYPERGRAPH PARAMETERS
# Document: BigBangHypergraphTheory_12Dec2025.docx
# Theory: 26D Egg Universe with SCm-UA encapsulation layers, Higgs as shift marker
# ═══════════════════════════════════════════════════════════════════════════════

COSMIC_EGG_HYPERGRAPH_PARAMS = {
    # Document metadata
    'document': 'BigBangHypergraphTheory_12Dec2025.docx',
    'grok_conversation': 'https://x.com/i/grok?conversation=1999530577406439555',
    'date': '12 December 2025',
    
    # 26D Egg Universe Structure
    'n_dimensions': 26,                     # 26-dimensional egg structure
    'n_cosmic_eggs': 26,                    # Hypergraph of nested cosmic eggs
    'n_encapsulation_layers': 5,            # [UA'] to [UA''''']
    
    # UA Encapsulation Layer Names
    'UA_layers': ['[UA]', "[UA']", "[UA'']", "[UA''']", "[UA'''']", "[UA''''']"],
    
    # UA Layer Density Fractions (halving progression)
    'UA_layer_fractions': [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125],
    
    # Speed parameters
    'v_BigBang_initial': 2.998e8,           # m/s - initial Big Bang speed (~c)
    'v_current_expansion': 0.0695,          # c - current Hubble flow ~6.95% c
    
    # Higgs as Inertial Gradient Shift Marker
    'VEV_Higgs_GeV': 246.0,                 # GeV - Higgs VEV
    'VEV_Higgs_J': 3.94e-8,                 # J - VEV in Joules
    'm_Higgs_GeV': 125.09,                  # GeV - Higgs mass
    'kappa_Higgs_shift': 0.01,              # Higgs contribution factor
    
    # DPM Grinding Parameters
    'omega_CW_SCm_north': 1.0e-10,          # rad/s - CW rotation (SCm, north pole)
    'omega_CCW_UA_south': -1.0e-10,         # rad/s - CCW rotation ([UA'], south pole)
    'grind_decay_factor': 0.1,              # Exponential decay per layer
    
    # Proto-Hydrogen 26D Shell
    'n_proto_H_shells': 26,                 # Empty shell alignments
    'E_proto_H_symbolic': 'c^26',           # Symbolic energy (infinity-like)
    'shell_fill_rate': 1e-15,               # s^-1 - Shell filling rate
    
    # Metallicity at [UA''''']
    'metallicity_peak_layer': 5,            # [UA'''''] = densest metallicity
    'superconductive_metal_formation': True,
    'globular_cluster_epoch': 1,            # 1st epoch black holes
    
    # Vacuum and Buoyancy Standards
    'rho_vac_UA_base': 7.09e-36,            # J/m³ - base UA vacuum density
    'rho_vac_SCm_base': 7.09e-37,           # J/m³ - base SCm vacuum density
    'buoyancy_pre_mass': True,              # Buoyancy established before mass
    
    # Big Bang Initiation Parameters
    'BBDT_factor': 1.0,                     # Big Bang Dynamo Term
    'f_TRZ_max_BigBang': 1.0,               # Maximum TRZ at Big Bang
    'Planck_density_J_m3': 5.16e96,         # J/m³ - Planck density
    
    # Time Relevancy Adjustments
    't_adj_formula': 't_obs / (1 + delta_rel)',
    'delta_rel_typical': 0.1,               # Typical relativistic dilation
    
    # Destruction Limit Philosophy
    'destruction_yields': 'one piece per complete action',
    'full_blueprint_possible': False,       # Insufficient matter in universe
    
    # Probability and Partition
    'partition_9D': 9,                      # 9D partition function base
    'entropy_26D_initial': 1e120,           # Initial 26D egg entropy (Planck units)
    
    # Reverse Analogy (Higgs flavors)
    'higgs_flavors_meaning': 'inside-out non-reversible representations',
    '2D_wreckage_planes': True,             # Collision debris as 2D projections
    
    # Cross-references to existing models
    'related_models': ['DPMModel', 'BigBangOriginModel', 'HiggsSCmIntegrationModel'],
}


def get_cosmic_egg_hypergraph_params() -> Dict[str, Any]:
    """
    Get parameters formatted for CosmicEggHypergraphModel.
    
    Returns dictionary compatible with:
        model = CosmicEggHypergraphModel()
        result = model.compute_26D_egg_energy()
    """
    p = COSMIC_EGG_HYPERGRAPH_PARAMS
    return {
        'n_dimensions': p['n_dimensions'],
        'n_encapsulation_layers': p['n_encapsulation_layers'],
        'UA_layers': p['UA_layers'],
        'UA_layer_fractions': p['UA_layer_fractions'],
        'v_init': p['v_BigBang_initial'],
        'v_current': p['v_current_expansion'] * p['v_BigBang_initial'],
        'VEV_Higgs': p['VEV_Higgs_J'],
        'omega_CW': p['omega_CW_SCm_north'],
        'omega_CCW': p['omega_CCW_UA_south'],
        'n_proto_H_shells': p['n_proto_H_shells'],
        'rho_vac_UA': p['rho_vac_UA_base'],
        'rho_vac_SCm': p['rho_vac_SCm_base'],
    }


def get_alpha_bec_lenr_params() -> Dict[str, Any]:
    """
    Get parameters formatted for AlphaBECModel.
    
    Returns dictionary compatible with:
        model = AlphaBECModel()
        result = model.verify_Tohsaki_AMD()
    """
    p = ALPHA_BEC_LENR_PARAMS
    return {
        'm_alpha': p['m_alpha'],
        'N_B_C12': p['N_B_C12'],
        'N_B_O16': p['N_B_O16'],
        'condensate_fraction_Hoyle': p['condensate_fraction_Hoyle'],
        'T_c_base': p['T_c_base_K'],
        'delta_T_c': p['delta_T_c_LENR_K'],
        'rho_BEC': p['rho_BEC_fm3'],
        'b_gaussian': p['b_gaussian_fm'],
    }


def get_solar_wind_parker_probe_params() -> Dict[str, Any]:
    """
    Get parameters formatted for SolarWindModel.
    
    Returns dictionary compatible with:
        model = SolarWindModel()
        model.verify_PSP_CDAWeb_2025()
    """
    p = SOLAR_WIND_PARKER_PROBE_PARAMS
    return {
        'delta_sw': p['delta_sw'],
        'v_sw': p['v_sw'],
        'n_p_range': (p['n_p_min'], p['n_p_max']),
        'v_sw_range': (p['v_sw_slow'], p['v_sw_fast']),
        'rho_expected': p['rho_sw_expected'],
        'sw_factor': p['sw_factor'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UQFF COMPRESSED SUMMARY PARAMETERS
# Document: Compressed Summary of Your Unified Quantum Field Equation System
# Complete 26-level polynomial structure, F_U equation, component equations
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_COMPRESSED_SUMMARY_PARAMS = {
    # Document metadata
    'document': 'Compressed Summary of Your Unified Quantum Field Equation System',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
    'date': '2025-09-29',
    
    # 26-Level Polynomial Nuclear Structure
    'E_0': 1e-20,                           # J - Base energy (smallest quantum energy)
    'n_levels': 26,                         # Total polynomial levels
    'E_n_formula': 'E_n = E_0 × 10^n',      # Energy scaling formula
    'total_span': 1e25,                     # Orders of magnitude span
    
    # 26-Level Energy Table
    '26_level_energies': {
        1: {'E_n': 1e-19, 'scale': 'Sub-quantum fluctuations'},
        2: {'E_n': 1e-18, 'scale': 'Planck-like vacuum'},
        3: {'E_n': 1e-17, 'scale': 'Weak interactions'},
        4: {'E_n': 1e-16, 'scale': 'Electron bindings'},
        5: {'E_n': 1e-15, 'scale': 'Atomic excitations'},
        6: {'E_n': 1e-14, 'scale': 'Nuclear gamma rays'},
        7: {'E_n': 1e-13, 'scale': 'Neutron bindings'},
        8: {'E_n': 1e-12, 'scale': 'Proton-neutron pairs'},
        9: {'E_n': 1e-11, 'scale': 'Alpha clusters'},
        10: {'E_n': 1e-10, 'scale': 'Atomic solids'},
        11: {'E_n': 1e-9, 'scale': 'Molecular'},
        12: {'E_n': 1e-8, 'scale': 'Macroscopic'},
        13: {'E_n': 1e-7, 'scale': 'Cosmic plasma'},
        14: {'E_n': 1e-6, 'scale': 'Low-energy astrophysics'},
        15: {'E_n': 1e-5, 'scale': 'Stellar winds'},
        16: {'E_n': 1e-4, 'scale': 'Planetary cores'},
        17: {'E_n': 1e-3, 'scale': 'Solar flares'},
        18: {'E_n': 1e-2, 'scale': 'Higgs boson'},
        19: {'E_n': 1e-1, 'scale': 'High-energy particles'},
        20: {'E_n': 1e0, 'scale': 'Galactic vacuum (Ug4)'},
        21: {'E_n': 1e1, 'scale': 'Black hole influences'},
        22: {'E_n': 1e2, 'scale': 'Quasar jets'},
        23: {'E_n': 1e3, 'scale': 'Galactic spins'},
        24: {'E_n': 1e4, 'scale': 'Intergalactic'},
        25: {'E_n': 1e5, 'scale': 'Cosmic rays'},
        26: {'E_n': 1e6, 'scale': 'Universal scales'},
    },
    
    # Coupling Constants (refined from solar data)
    'k_1': 1.5,                             # Ug1 coupling (dipole)
    'k_2': 1.2,                             # Ug2 coupling (heliosphere)
    'k_3': 1.8,                             # Ug3 coupling (strings disk)
    'k_4': 1.0,                             # Ug4 coupling (star-black hole)
    'beta_i': 0.6,                          # Buoyancy coupling
    'lambda_weights': [1.0, 0.8, 1.2, 0.5], # Inertia couplings
    
    # Galactic Parameters
    'omega_g': 7.3e-16,                     # rad/s - Galactic spin
    'M_bh': 8.15e36,                        # kg - Sgr A* black hole mass
    'd_g': 2.44e20,                         # m - Sun-Sgr A* distance (verified VERA/GAIA)
    'd_g_uqff': 2.55e20,                    # m - UQFF approximation (5% error)
    
    # Reactor Efficiency
    'E_react_0': 1e46,                      # W/m³ - Initial reactor efficiency
    'kappa': 0.0005,                        # day^-1 - Decay rate
    'E_react_formula': 'E_react = 10^46 × e^(-0.0005t)',
    'tau_days': 2000,                       # days - Characteristic time
    
    # Magnetic Parameters
    'mu_j_base': 1e3,                       # T·m³ - Base magnetic moment
    'mu_j_osc': 0.4,                        # Oscillation amplitude
    'mu_j_scale': 3.38e20,                  # T·m³ - Scaling factor
    'r_j': 1.496e13,                        # m - String distance (1 AU base)
    'gamma_decay': 5e-5,                    # day^-1 - Decay rate
    
    # Solar Wind Parameters
    'delta_sw': 0.01,                       # Wind modulation factor
    'v_sw': 5e5,                            # m/s - Wind velocity
    'H_SCm': 1.0,                           # Heliosphere factor
    
    # Aether Parameters
    'g_mu_nu': [1, -1, -1, -1],             # Minkowski metric signature
    'eta': 1e-22,                           # Aether coupling
    'T_s_mu_nu_total': 1.112e7,             # J/m³ - Total stress-energy
    'T_s_mu_nu_UA': 1.27e3,                 # J/m³ - UA contribution
    'T_s_mu_nu_SCm': 1.11e7,                # J/m³ - SCm contribution
    
    # Vacuum Energy
    'lambda_vac_formula': 'λ_vac = Σ(f_i E_i)/V',
    'lambda_vac_cosmic': 1e-9,              # J/m³ - Cosmic scale
    'lambda_vac_SCm': 7.09e-37,             # J/m³ - SCm density
    'lambda_vac_UA': 7.09e-36,              # J/m³ - UA density
    
    # SCm Properties
    'SCm_density': 1e15,                    # kg/m³ - No quantum signature
    'SCm_Qs': None,                         # No quantum signature
    'UA_charge': 1e-11,                     # C - Trapped aether
    'rho_A': 1e-23,                         # kg/m³ - Aether density
    
    # Core Penetration Factors
    'P_core_Sun': 1.0,                      # Full penetration for Sun
    'P_core_planet': 1e-3,                  # Reduced for planets
    'P_SCm_Sun': 1.0,                       # Full SCm penetration for Sun
    'P_SCm_planet': 1e-3,                   # Reduced for planets
    
    # Time Parameters
    't_n_definition': 't - t_0 (can be negative for reversals)',
    'omega_cycle': np.pi,                   # rad/s - Cycle constant
    'alpha_decay': 0.001,                   # day^-1 - Time decay
    
    # Solar Reference Values (t=0, t_n=0)
    'Ug1_Sun': 1.39e26,                     # J/m³
    'Ug2_Sun': 1.18e53,                     # J/m³
    'Ug3_Sun': 1.8e49,                      # J/m³
    'Ug4_Sun': 2.50e-20,                    # J/m³
    'Um_Sun': 2.28e65,                      # J/m³ (DOMINANT)
    'UI_Sun': 1.38e-47,                     # J/m³
    'A_mu_nu_Sun': 1.12e-15,                # J/m³
    'Ub1_Sun': -1.94e27,                    # J/m³ (opposes Ug)
    'F_U_Sun': 2.28e65,                     # J/m³ (dominated by Um)
    
    # Verification Data
    'verification': {
        'Sun_SgrA_distance_verified': {'uqff': 2.55e20, 'vera_gaia': 2.44e20, 'error_pct': 5},
        'quasar_luminosity_range': (1e39, 1e47),  # W
        'E_react_matches': True,
        'vacuum_cosmological_constant': 1e-9,  # J/m³
        'nuclear_binding_n8': 1e-12,  # J matches
        'F_U_normalized': 1e27,  # N/m²
    },
    
    # Key Concepts Summary
    'key_concepts': {
        'SCm': 'Dense, undetectable (no Qs), drives Ug3, quasar jets, planetary cores',
        'Ug_ranges': 'Ug1 (dipole), Ug2 (heliosphere), Ug3 (disk strings), Ug4 (star-BH)',
        'Um': 'Near-lossless magnetic strings in 90° disk',
        'Ub': 'Opposes Ug, proportional to galactic spin/BH strength',
        'UA': 'Medium for interactions; negative time for reversals',
        'Quasars': 'Ug failure to trap SCm (reactive/fastest substance)',
        '26_levels': 'E_n = E_0 × 10^n, polynomial hierarchy n=1-26',
    },
    
    # Related Existing Models
    'related_models': [
        'UnifiedFieldEquation',
        'CompleteUnifiedFieldModel',
        'ReactorEfficiencyModel',
        'FermiLAT4LACBlazarModel',
        'NuclearBindingShellLevelsModel',
    ],
}


def get_uqff_compressed_summary_params() -> Dict[str, Any]:
    """
    Get parameters for UQFF Compressed Summary calculations.
    
    Returns dictionary compatible with:
        model = UnifiedFieldEquation()
        result = model.compute_F_U(params)
    """
    p = UQFF_COMPRESSED_SUMMARY_PARAMS
    return {
        'E_0': p['E_0'],
        'n_levels': p['n_levels'],
        'k_i': [p['k_1'], p['k_2'], p['k_3'], p['k_4']],
        'beta_i': p['beta_i'],
        'omega_g': p['omega_g'],
        'M_bh': p['M_bh'],
        'd_g': p['d_g'],
        'E_react': p['E_react_0'],
        'kappa': p['kappa'],
        'eta': p['eta'],
        'delta_sw': p['delta_sw'],
        'v_sw': p['v_sw'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# POLYNOMIAL NUCLEAR STRUCTURE 2025 DATASET VERIFICATION
# Document: 26-Level Polynomial Verification with High-Energy Datasets (2025)
# Verification against: Fermi LAT, Chandra, Parker, Voyager, Gaia, ENSDF, LHC
# References:
#   - HEASARC: https://heasarc.gsfc.nasa.gov/
#   - Fermi LAT 4LAC-DR4: https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR4/
#   - Chandra: https://cxc.harvard.edu/
#   - Parker Solar Probe CDAWeb: https://cdaweb.gsfc.nasa.gov/
#   - Voyager: https://voyager.jpl.nasa.gov/
#   - Gaia DR3/DR4: https://gea.esac.esa.int/archive/
#   - ENSDF (NNDC): https://www.nndc.bnl.gov/ensdf/
#   - ATLAS-CONF-2025-007: https://cds.cern.ch/
# ═══════════════════════════════════════════════════════════════════════════════

DATASET_VERIFICATION_2025_PARAMS = {
    # Document metadata
    'document': '26-Level Polynomial Verification with High-Energy Datasets (2025)',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
    'date': '2025-09-29',
    'verification_method': 'Multi-dataset cross-validation',
    
    # ─────────────────────────────────────────────────────────────────────────
    # POLYNOMIAL FIT PARAMETERS (Sample Pb-206 ENSDF levels)
    # Code: import pandas as pd; levels = [0, 0.044, 0.137, 0.334, 0.583, 0.802, 1.028] * 1.6e-13
    #       poly = np.polyfit(range(len(levels)), levels, 26); print(poly)
    # ─────────────────────────────────────────────────────────────────────────
    'Pb206_sample_levels_MeV': [0.0, 0.044, 0.137, 0.334, 0.583, 0.802, 1.028],  # Sample levels
    'Pb206_sample_levels_J': [
        0.0, 0.044 * 1.6e-13, 0.137 * 1.6e-13, 0.334 * 1.6e-13,
        0.583 * 1.6e-13, 0.802 * 1.6e-13, 1.028 * 1.6e-13
    ],
    'poly_deg_test': 26,                    # Polynomial degree tested
    'R_squared_deg26': 1.0,                 # Perfect fit (overfit)
    'R_squared_low_deg': 0.95,              # Physical fit quality
    'poly_fit_note': 'R²~0.95 for low deg, overfit for deg=26; no standard 26-deg polynomial in shell models per NNDC/IAEA',
    
    # ─────────────────────────────────────────────────────────────────────────
    # 26-LEVEL ENERGY SCALE VERIFICATION TABLE
    # E_n = E_0 × 10^n (E_0 = 10^{-20} J)
    # ─────────────────────────────────────────────────────────────────────────
    'energy_level_verification': {
        # Sub-nuclear (n=1-5): Verified vs LHC quark energies
        'n1_5': {
            'E_range_J': (1e-19, 1e-15),
            'scale': 'Sub-nuclear',
            'verification': 'Aligns with LHC quark energies ~10^-16 J from ATLAS-CONF-2025-007',
            'verified': True,
        },
        # Nuclear bindings (n=6-10): Verified vs ENSDF A=206
        'n6_10': {
            'E_range_J': (1e-14, 1e-10),
            'scale': 'Nuclear bindings',
            'verification': 'Fits ENSDF A=206 levels: ground to ~10^-12 J; polynomial fit R²~0.95',
            'ENSDF_A': 206,
            'ENSDF_max_MeV': 10,            # ~10 MeV = 10^-12 J matches n=8
            'verified': True,
        },
        # Excitations/molecular (n=11-15): Verified vs LHC ion collisions
        'n11_15': {
            'E_range_J': (1e-9, 1e-5),
            'scale': 'Excitations/molecular',
            'verification': 'Matches LHC ion collisions, arXiv:2504.00790',
            'arxiv': '2504.00790',
            'verified': True,
        },
        # Stellar/plasma (n=16-20): Verified vs Parker solar wind
        'n16_20': {
            'E_range_J': (1e-4, 1e0),
            'scale': 'Stellar/plasma',
            'verification': 'Aligns Parker solar wind energies ~10^-6 J/proton',
            'parker_E_per_proton_J': 1e-6,
            'verified': True,
        },
        # Galactic (n=21-26): Verified vs Fermi quasar jets
        'n21_26': {
            'E_range_J': (1e1, 1e6),
            'scale': 'Galactic',
            'verification': 'Fits Fermi quasar jets ~10^6 J events; no direct 26-level evidence',
            'fermi_quasar_E_J': 1e6,
            'verified': False,              # Speculative at high-n
        },
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # COMPUTED Ug VALUES FOR SUN/SGR A* (2025 DATA)
    # ─────────────────────────────────────────────────────────────────────────
    'Ug1_Sun_SgrA': 9.26e22,                # e^{-0.001 t} (Fermi solar flare alignments)
    'Ug1_alpha': 0.001,                     # day^-1 decay rate
    'Ug1_verification': 'Fermi solar flare energies align',
    
    'Ug2_Sun_SgrA': 8.91e6,                 # Parker wind v_sw=5×10^5 m/s
    'Ug2_v_sw': 5e5,                        # m/s solar wind velocity
    'Ug2_delta_sw': 0.01,                   # Wind modulation factor
    'Ug2_verification': 'Parker CDAWeb wind density ~8×10^-21 kg/m³',
    
    'Ug3_Sun_SgrA': 1e3,                    # cos(2.5×10^-6 t) (Chandra magnetic fields)
    'Ug3_omega_s': 2.5e-6,                  # rad/s spin frequency
    'Ug3_verification': 'Chandra magnetic field measurements',
    
    'Ug4_Sun_SgrA': 3.19e16,                # Gaia Sgr A* M_bh=4.1×10^6 M_☉
    'Ug4_M_bh_Msun': 4.1e6,                 # M_☉ (Gaia DR3)
    'Ug4_verification': 'Gaia DR3 Sgr A* mass measurement',
    
    'Ubi_Sun_SgrA': -1.08e23,               # e^{-0.001 t} (JCAP DM spike data)
    'Ubi_verification': 'JCAP dark matter spike profile at MW center',
    
    'Um_Sun_SgrA': 2.26e16,                 # (1 - e^{-0.0001 t}) (Fermi blazar jets)
    'Um_gamma': 0.0001,                     # day^-1 gamma decay
    'Um_verification': 'Fermi blazar jets magnetic energy transfer',
    
    'UA_mu_nu': [1, -1, -1, -1],            # Minkowski metric
    'UA_eta_correction': 1.27e-20,          # Cosmological λ_vac ~10^-9 J/m³ from JCAP
    
    # ─────────────────────────────────────────────────────────────────────────
    # VARIABLE TABLE WITH DATA VERIFICATIONS
    # ─────────────────────────────────────────────────────────────────────────
    'variable_table': {
        'k_i': {
            'value': (1.2, 1.8),
            'tuning': 'LHC/ATLAS Higgs m_H=125 GeV ~10^-8 J (n=12)',
            'source': 'ATLAS-CONF-2025-007',
        },
        'lambda_vac': {
            'value': 1e-9,
            'description': 'Σ f_i E_i / V',
            'verification': 'Matches JCAP MW DM spike; no SCm evidence',
            'unit': 'J/m³',
        },
        'SCm_density': {
            'value': 1e15,
            'unit': 'kg/m³',
            'note': 'Speculative; fastest v_SCm~10^8 m/s fits Fermi jets',
            'v_SCm': 1e8,                   # m/s
        },
        't_n': {
            'description': 't - t_0 (can be <0 for reversals)',
            'asymmetry_source': 'Chandra 3C 273 quasar jet asymmetry',
        },
        'E_react': {
            'formula': '10^46 × e^{-0.0005 t}',
            'observed_range_W': (1e39, 1e47),
            'source': 'Fermi 4LAC blazars',
        },
        'R_b': {
            'value': 1.496e13,              # m (not AU, this is ~100 AU)
            'description': 'Heliosphere boundary',
            'source': 'Parker/Voyager',
        },
        'd_g': {
            'value_gaia': 2.47e20,          # m
            'value_uqff': 2.55e20,          # m (5% error)
            'description': 'Sun to Sgr A* distance',
            'source': 'Gaia DR3 (2025 update)',
        },
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # DATASET VERIFICATION SUMMARY
    # ─────────────────────────────────────────────────────────────────────────
    'verification_summary': {
        'nuclear_polynomial': {
            'source': 'ENSDF (NNDC)',
            'status': 'PARTIAL',
            'note': 'A=206 ~20-30 excitations up to 10 MeV; polynomial fit overfits (R²=1 for deg=26, but unphysical vs shell model ~10 levels)',
            'lhc_support': 'ATLAS-CONF-2025-007 Higgs/off-shell data - no 26-level support (speculative)',
        },
        'quasars_jets': {
            'source': 'Fermi LAT 4LAC (HEASARC)',
            'status': 'VERIFIED',
            'observed': '90% gamma sources are blazars',
            'luminosity_match': 'E_react matches 4LAC luminosities',
            'chandra_source': 'RACS J0320-35 (2025) jet growth fluid/unequal as predicted',
            'scm_detection': 'No SCm detection',
        },
        'heliosphere_solar_wind': {
            'source': 'Parker CDAWeb + Voyager',
            'status': 'VERIFIED',
            'parker_density': 8e-21,        # kg/m³ (aligns Ug2)
            'voyager_boundary_AU': 122,     # AU (interstellar boundary)
            'age_correlation': 'Speculative liquid link unverified',
        },
        'sgr_a_vacuum': {
            'source': 'Gaia DR3/DR4 + JCAP',
            'status': 'PARTIAL',
            'gaia_note': 'DR4 (mid-2026 preview 2025) star orbits near Sgr A show no Ug4 signature',
            'jcap_dm_density': 1e-9,        # J/m³ (fits λ_vac)
        },
        'overall': 'Model internally consistent; partial empirical alignments (energies, jets)',
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # MILLENNIUM PROBLEM TIES
    # Key concepts referenced: Navier-Stokes, Yang-Mills, Riemann
    # ─────────────────────────────────────────────────────────────────────────
    'millennium_ties': {
        'navier_stokes': 'Quasar jets modeled as fluid dynamics (unequal jets, SCm expulsion)',
        'yang_mills': 'SCm mass gap (no quantum signature Qs, explains confinement)',
        'riemann': 'π cycles in resonance frequencies (26-level periodicity)',
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # KEY CONCEPTS SUMMARY FROM DOCUMENT
    # ─────────────────────────────────────────────────────────────────────────
    'key_concepts': {
        '26_level_polynomial': 'V(r) ≈ Σ_{n=1}^{26} a_n r^n, exponentially scaled E_n = E_0 × 10^n',
        'SCm_cosmic_glue': 'Dense, no Qs, Ug3 driver, quasar/core interactions',
        'UA_medium': 'Universal Aether with negative time derivations (UA prime, UA double-prime) for reversals',
        'Ug_failure': 'Quasars form when Ug fails to trap SCm, igniting UA in fluid unequal jets',
        'heliosphere_Ug2': 'Transmutates winds to liquids, thickness correlates planetary volumes to stellar age',
        'planetary_cores': 'SCm + UA exclusively interact with Ug3 for orbits/spins',
        'vacuum_energy': 'SCm/UA inertia influences λ_vac',
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # PHENOMENA EXPLAINED
    # ─────────────────────────────────────────────────────────────────────────
    'phenomena_explained': {
        'heliosphere': 'Ug2 transmutates winds to liquids; thickness + planetary volumes correlate to stellar age',
        'quasars': 'Ug failure expels SCm; ignites UA in fluid unequal jets',
        'planetary_cores': 'SCm + UA interact exclusively with Ug3 for orbits/spins',
        'vacuum_energy': 'SCm/UA inertia influences λ_vac (cosmological constant)',
    },
    
    # Related models in CondensedPhysics.py
    'related_models': [
        'NuclearBindingShellLevelsModel',
        'FermiLAT4LACBlazarModel',
        'HelioUg2Model',
        'SagittariusAStarModel',
        'QuasarJetFluidModel',
    ],
}


def get_dataset_verification_2025_params() -> Dict[str, Any]:
    """
    Get parameters formatted for 2025 Dataset Verification.
    
    Returns dictionary compatible with cross-validation against:
        - ENSDF/NNDC nuclear data
        - Fermi LAT 4LAC blazars
        - Chandra X-ray observations
        - Parker Solar Probe wind data
        - Voyager interstellar boundary
        - Gaia DR3/DR4 Sgr A* distance
    """
    p = DATASET_VERIFICATION_2025_PARAMS
    return {
        # Polynomial fit data
        'Pb206_levels': p['Pb206_sample_levels_MeV'],
        'poly_R_squared': p['R_squared_low_deg'],
        
        # Computed Ug values
        'Ug1': p['Ug1_Sun_SgrA'],
        'Ug2': p['Ug2_Sun_SgrA'],
        'Ug3': p['Ug3_Sun_SgrA'],
        'Ug4': p['Ug4_Sun_SgrA'],
        'Ubi': p['Ubi_Sun_SgrA'],
        'Um': p['Um_Sun_SgrA'],
        
        # Key physical parameters
        'v_sw': p['Ug2_v_sw'],
        'M_bh': p['Ug4_M_bh_Msun'] * SOLAR_MASS_KG,
        'd_g': p['variable_table']['d_g']['value_gaia'],
        'lambda_vac': p['variable_table']['lambda_vac']['value'],
        'E_react_range': p['variable_table']['E_react']['observed_range_W'],
        
        # Verification status
        'verification': p['verification_summary'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# TSP Q-SCOPE SUPERCONDUCTIVE FRAMEWORK PARAMETERS
# Document: Universal Quantum Field Superconductive Framework (UQFF/TSP)
# Theory of Superconductive Permanence with Q-Scope Data from Groups #1-12
# References:
#   - HTSC-2025: https://www.htsc-2025.info/
#   - THz Superconductors: arXiv:2403.xxxxx (quantum echo 2025)
#   - Ultrafast Oscilloscope: https://www.keysight.com/oscilloscopes
#   - THz Astronomy: https://science.nasa.gov/astrophysics/
# ═══════════════════════════════════════════════════════════════════════════════

TSP_QSCOPE_SUPERCONDUCTIVE_PARAMS = {
    # Document metadata
    'document': 'Universal Quantum Field Superconductive Framework (UQFF/TSP)',
    'grok_conversation': 'https://x.com/i/grok?conversation=1972174124559557054',
    'date': '2025-05-15',
    'corpus_size': '~7000 pages (compressed to ~15%)',
    
    # ─────────────────────────────────────────────────────────────────────────
    # Q-SCOPE DATA PARAMETERS (Groups #1-12, Images 1-73)
    # Two-channel measurements over ~38.4 seconds
    # ─────────────────────────────────────────────────────────────────────────
    
    # Channel 1: Smooth Q-Wave (Variable Frequency/Amplitude)
    'A1_initial': 0.491,                    # V - Initial amplitude (Group #1)
    'A1_DPM': 0.4604,                       # V - Average amplitude (DPM calibration)
    'f_primary_initial': 5455,              # Hz - Initial frequency (5.455 kHz)
    'f_primary_final': 976.68,              # Hz - Final frequency (Group #12)
    'dT_initial_ms': 8,                     # ms - Initial period
    'dT_final_ms': 25,                      # ms - Final period
    'f_dT_initial': 125,                    # Hz (1/dT = 1/0.008)
    'f_dT_final': 40,                       # Hz (1/dT = 1/0.025)
    
    # Channel 2: Eccentric (Stable Amplitude)
    'A2_channel': 3.102,                    # V - Constant amplitude (flux-pinned)
    'A2_stability': 'constant',             # Stable across Groups #1-12
    
    # Total evolution time
    't_evolution_sec': 38.4,                # s - Total measurement time
    
    # ─────────────────────────────────────────────────────────────────────────
    # 1.2 THz HOLE PARAMETERS
    # Low-energy anomaly facilitating signal reversal and stabilization
    # ─────────────────────────────────────────────────────────────────────────
    'f_THz_hole': 1.2e12,                   # Hz - THz hole frequency
    'f_THz_hole_min': 1.2e12,               # Hz - Minimum
    'f_THz_hole_max': 1.3e12,               # Hz - Maximum
    'THz_hole_function': 'Low-energy reversal and stabilization',
    'THz_hole_verified': False,             # Not empirically verified as of 2025
    
    # ─────────────────────────────────────────────────────────────────────────
    # GINZBURG-LANDAU PARAMETERS (Ug - Order Parameter ψ)
    # ∇²ψ + αψ + β|ψ|²ψ = 0
    # ─────────────────────────────────────────────────────────────────────────
    'alpha_GL': -1e6,                       # α ∝ (T - T_c) < 0 below T_c
    'beta_GL': 1e12,                        # β > 0 (nonlinear coefficient)
    'psi_stable': 1.0,                      # |ψ|² ≈ 1 from A₂ = 3.102 V
    'coherence_length_m': 1e-6,             # ξ - Coherence length (m)
    
    # ─────────────────────────────────────────────────────────────────────────
    # BOGOLIUBOV-DE GENNES PARAMETERS (Ub - Quasiparticle)
    # ─────────────────────────────────────────────────────────────────────────
    'hbar': 1.055e-34,                      # ℏ (J·s)
    'm_electron': 9.109e-31,                # kg
    'Delta_gap_initial': 5e-22,             # J - Gap energy (initial)
    'Delta_gap_final': 1.6e-22,             # J - Gap energy (from k_Δ × f_dT=40)
    'k_Delta': 4e-24,                       # Gap-frequency coupling (J/Hz)
    'E_quasiparticle': 'constant',          # Stable from A₂ constant
    
    # ─────────────────────────────────────────────────────────────────────────
    # FLUX PINNING PARAMETERS (Um - Magnetic Flux)
    # Um = Φ₀ Σ_i δ(r - r_i)
    # ─────────────────────────────────────────────────────────────────────────
    'Phi_0': 2.067833848e-15,               # Wb - Flux quantum
    'vortex_positions_fixed': True,         # Fixed r_i from A₂ stability
    'pinning_energy': 1e-17,                # J - Typical pinning energy
    
    # ─────────────────────────────────────────────────────────────────────────
    # Q-WAVE RESONANCE PARAMETERS (Ur)
    # Ur = A sin(2πft) + A₂ sin(2πft + ϕ)
    # ─────────────────────────────────────────────────────────────────────────
    'A_Ur': 0.491,                          # V - Primary amplitude
    'A2_Ur': 3.102,                         # V - Secondary amplitude
    'f_Ur': 976.68,                         # Hz - Resonance frequency
    'phi_Ur': 0,                            # rad - Phase difference
    'dA_UA': 2.611,                         # V - Amplitude difference (A₂ - A₁)
    
    # ─────────────────────────────────────────────────────────────────────────
    # DI-PSEUDO-MONOPOLE (U_dp) PARAMETERS
    # U_dp = k × (A₁ × A₂) / f_dp² × cos(ϕ_dp)
    # ─────────────────────────────────────────────────────────────────────────
    'k_DPM': 6.674e-11,                     # m³/kg/s² (= G)
    'A1_DPM_final': 0.491,                  # V
    'A2_DPM_final': 3.102,                  # V
    'f_dp': 40,                             # Hz (= f_dT final)
    'phi_dp': 0,                            # rad - Phase
    'U_dp_computed': 0.0008926,             # (A₁ × A₂) / f_dp² at Group #12
    
    # ─────────────────────────────────────────────────────────────────────────
    # TEMPORAL EVOLUTION (Ut)
    # Ut = 1/dT
    # ─────────────────────────────────────────────────────────────────────────
    'Ut_initial': 125,                      # Hz (1/8ms)
    'Ut_final': 40,                         # Hz (1/25ms)
    'Ut_trend': 'slowing',                  # Indicates cooling/stabilization
    
    # ─────────────────────────────────────────────────────────────────────────
    # SC_m COHERENCE METRIC
    # SC_m = |ψ|² / ∫|ψ|² dV
    # ─────────────────────────────────────────────────────────────────────────
    'SC_m_stable': 1.0,                     # ≈1 from A₂ stability
    
    # ─────────────────────────────────────────────────────────────────────────
    # BRAIN WAVE SUBHARMONIC MAPPING
    # f_sub = f/n → brain wave bands
    # ─────────────────────────────────────────────────────────────────────────
    'brainwave_bands': {
        'delta': (0.5, 4),                  # Hz - Deep sleep
        'theta': (4, 8),                    # Hz - Meditation/relaxation
        'alpha': (8, 13),                   # Hz - Relaxation/calm focus
        'beta': (13, 30),                   # Hz - Active thinking
        'gamma': (30, 100),                 # Hz - High activity/cognition
    },
    'f_sub_example': 976.68 / 20,           # = 48.834 Hz (gamma range)
    'subharmonic_factor_range': (10, 120),  # n values for brain wave mapping
    
    # Emotional state mapping
    'emotion_mapping': {
        'gamma': 'High activity, peak cognition, problem-solving',
        'beta': 'Active thinking, alertness, anxiety',
        'alpha': 'Relaxation, calm focus, meditation entry',
        'theta': 'Deep relaxation, creativity, dreaming',
        'delta': 'Deep sleep, regeneration, unconscious',
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # NAVIER-STOKES SUPPORT (Vortex Dynamics)
    # ρ(v·∇v) = -∇p + μ∇²v
    # ─────────────────────────────────────────────────────────────────────────
    'rho_fluid': 6e3,                       # kg/m³ - Typical superconductor density
    'mu_viscosity': 1e-6,                   # Pa·s - Effective viscosity
    'flow_transition': 'turbulent → laminar',
    'dT_slowing_mechanism': 'Vortex pinning, magnetic stabilization',
    
    # ─────────────────────────────────────────────────────────────────────────
    # 26-LEVEL POLYNOMIAL TIES (Root Structure)
    # E_n = E_0 × 10^n (E_0 = 10^{-20} J)
    # ─────────────────────────────────────────────────────────────────────────
    '26_level_ties': {
        'n1_5': 'Sub-quantum (LHC quark ~10^{-16} J)',
        'n6_10': 'Nuclear (ENSDF Pb-206 ~10^{-12} J)',
        'n11_15': 'Molecular/THz (q-scope kHz-THz translations)',
        'n16_20': 'Stellar/plasma (Parker wind ~10^{-6} J)',
        'n21_26': 'Cosmic (Fermi jets ~10^6 J; THz emissions)',
    },
    'ramanujan_extension': 'Extends 6-10th level polynomials to unify quantum/cosmic',
    
    # ─────────────────────────────────────────────────────────────────────────
    # MILLENNIUM PROBLEM TIES
    # ─────────────────────────────────────────────────────────────────────────
    'millennium_ties': {
        'navier_stokes': 'Smoothness in laminar regimes from slowing dT',
        'yang_mills': 'Mass gap via 1.2 THz hole (low-energy states)',
        'riemann': 'π cycles in resonance, encoding primes in q-waves',
    },
    
    # ─────────────────────────────────────────────────────────────────────────
    # VERIFICATION STATUS
    # ─────────────────────────────────────────────────────────────────────────
    'verification': {
        'q_scope_internal': True,           # Internal consistency high
        'thz_superconductors': 'Partial',   # THz experiments exist, not exact "hole"
        'brain_wave_links': False,          # Speculative, unverified
        '1.2_THz_hole': False,              # No empirical verification
        'htsc_2025_alignment': True,        # Catalog aligns Tc models
    },
    
    # Related models in CondensedPhysics.py
    'related_models': [
        'GinzburgLandauModel',
        'BogoliubovDeGennesModel',
        'BrainWaveSubharmonicModel',
        'DPMAttractionModel',
        'FluxPinningModel',
        'QWaveResonanceModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# FINAL PARSEC PROBLEM PARAMETERS (Drawing 3)
# ═══════════════════════════════════════════════════════════════════════════════
# SMBH binary merger dynamics at ~1 parsec separation
# UQFF [SCm]-[UA] mechanism resolves stalling
# ═══════════════════════════════════════════════════════════════════════════════

FINAL_PARSEC_PROBLEM_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # DOCUMENT METADATA
    # ───────────────────────────────────────────────────────────────────────────
    'document': 'Final Parsec Problem - SMBH Binary Merger Dynamics',
    'source': 'Star Magic / UQFF Framework (Drawing 3)',
    'date': '2025-02-09',
    'grok_session': 'https://x.com/i/grok',
    
    # ───────────────────────────────────────────────────────────────────────────
    # PHYSICAL SCALES
    # ───────────────────────────────────────────────────────────────────────────
    'parsec_m': 3.086e16,                    # 1 parsec in meters
    'final_parsec_separation': 3.086e16,     # ~1 pc stalling distance (m)
    'gw_regime_separation': 3.086e15,        # 0.1 pc - GW dominant (m)
    'dynamical_friction_cutoff': 3.086e19,   # ~1 kpc (m)
    't_Hubble_s': 4.35e17,                   # Hubble time (s) = 13.8 Gyr
    't_Hubble_yr': 1.38e10,                  # Hubble time (years)
    
    # ───────────────────────────────────────────────────────────────────────────
    # SMBH PARAMETERS (typical merger scenario)
    # ───────────────────────────────────────────────────────────────────────────
    'M_SMBH_range_kg': (1.989e36, 1.989e39), # 10^6 - 10^9 M_☉
    'M_SMBH_typical_kg': 1.989e38,           # 10^8 M_☉
    'M_SMBH_typical_Msun': 1e8,
    'mass_ratio_range': (0.01, 1.0),         # q = M_2 / M_1
    
    # ───────────────────────────────────────────────────────────────────────────
    # MERGER PHASES (Classical Model)
    # ───────────────────────────────────────────────────────────────────────────
    'phases': {
        'dynamical_friction': {
            'separation_range': (3.086e16, 3.086e19),  # 1 pc to 1 kpc
            'mechanism': 'Dynamical friction against stars/gas/DM',
            'timescale': '~1 Gyr',
            'efficiency': 'High when stellar density sufficient',
        },
        'binary_hardening': {
            'separation_range': (3.086e16, 3.086e17),  # 1 pc to 10 pc
            'mechanism': 'Three-body stellar ejection (slingshots)',
            'timescale': '~100 Myr',
            'efficiency': 'Decreases as loss cone depletes',
        },
        'final_parsec_stall': {
            'separation_range': (3.086e15, 3.086e16),  # 0.1 pc to 1 pc
            'mechanism': 'Loss cone depleted; GW too weak',
            'timescale': '> τ_Hubble (classical)',
            'efficiency': 'STALLED - THE PROBLEM',
        },
        'gw_inspiral': {
            'separation_range': (0, 3.086e15),  # 0 to 0.1 pc
            'mechanism': 'Gravitational wave emission',
            'timescale': 'Years to Myr',
            'efficiency': 'High once entered',
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # CLASSICAL PROPOSED SOLUTIONS
    # ───────────────────────────────────────────────────────────────────────────
    'classical_solutions': {
        'gas_rich_environment': {
            'mechanism': 'Circumbinary gas disk - viscous drag',
            'effectiveness': 'Works in wet mergers',
            'issue': 'Not all mergers are gas-rich',
        },
        'stellar_refilling': {
            'mechanism': 'Triaxial galaxy - chaotic orbits refill loss cone',
            'effectiveness': 'Partial resolution',
            'issue': 'Requires specific galaxy morphology',
        },
        'triple_smbh': {
            'mechanism': 'Kozai-Lidov oscillations from third BH',
            'effectiveness': 'Induces eccentricity → faster inspiral',
            'issue': 'Requires third BH presence',
        },
        'dark_matter': {
            'mechanism': 'Self-interacting DM provides friction',
            'effectiveness': 'Under active research (2024-2025)',
            'issue': 'Requires exotic DM properties',
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # UQFF [SCm]-[UA] SOLUTION PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    'uqff_solution': {
        # [SCm] expulsion mechanism
        'scm_expulsion': {
            'description': '[SCm] from SMBHs ignites against [UA]',
            'outcome': 'Quasar jet-like energy release',
            'effect': 'Angular momentum extraction from binary',
        },
        # [UA] string interactions
        'ua_string_interactions': {
            'description': 'Aether strings mediate energy transfer',
            'mechanism': 'Bidirectional longitudinal wavepairs (Whittaker)',
            'symmetry': '4-symmetry breaks 3-symmetry stalling',
        },
        # Time-Reversal Zones
        'trz_enhancement': {
            'description': 'Bearden TRZs enable negentropic reordering',
            'f_TRZ': 0.1,  # Enhancement factor
            'effect': 'Vacuum energy extraction accelerates merger',
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # VACUUM ENERGY DENSITIES (Level 13 - Cosmic Scale)
    # ───────────────────────────────────────────────────────────────────────────
    # ρ_vac,X = Σf_i E_i / V_object (J/m³)
    'rho_vac_SCm_BH': 2.39e-22,              # Black hole scale
    'rho_vac_UA_BH': 7.09e-36,               # [UA] background at BH
    'rho_vac_Ug4_BH': 1.19e-24,              # BH interaction energy
    'quantum_level': 13,                      # Cosmic, plasma-dominated
    
    # ───────────────────────────────────────────────────────────────────────────
    # Ug4 EQUATION PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    # Ug4 = k_4 ρ_vac,[SCm] (M_bh / d_g) e^(-αt) cos(ωt_n) (1 + f_feedback)
    'k_4': 1.0,
    'alpha': 0.0005,                         # Decay constant (1/day)
    'omega': np.pi,                          # Angular frequency
    'f_feedback_base': 0.1,                  # Feedback enhancement
    'f_TRZ': 0.1,                            # TRZ enhancement
    
    # ───────────────────────────────────────────────────────────────────────────
    # GW TIMESCALE PARAMETERS (Peters formula)
    # ───────────────────────────────────────────────────────────────────────────
    # τ_GW = (5/256) × (c⁵a⁴) / (G³M₁M₂(M₁+M₂))
    'G': 6.674e-11,
    'c': 2.998e8,
    'M_sun': 1.989e30,
    
    # ───────────────────────────────────────────────────────────────────────────
    # CGM METAL RETENTION (Drawing 31 - Sanchez et al.)
    # ───────────────────────────────────────────────────────────────────────────
    'cgm_metal_retention': {
        'f_Z_halo_overmassive': 0.89,        # Metal fraction - over-massive SMBH
        'f_Z_halo_undermassive': 0.85,       # Metal fraction - under-massive
        'f_Z_disk_starforming': 0.73,        # Disk - star-forming
        'f_Z_disk_quenched': 0.51,           # Disk - quenched
        'f_CGM_typical': 0.15,               # CGM baryon fraction
        # M-σ relation deviation drives feedback
        # δM_BH > 0 → metal ejection; δM_BH < 0 → retention
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # OBSERVATIONAL REFERENCES
    # ───────────────────────────────────────────────────────────────────────────
    'observational': {
        'lisa': 'Future: LISA for SMBH mergers (mHz)',
        'ligo_virgo': 'Current: stellar-mass BH mergers',
        'pta': 'NANOGrav 15-yr: stochastic GW background hints',
        'm87': 'M87 SMBH imaging (EHT)',
        'sgr_a': 'Sgr A* SMBH imaging (EHT)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # KEY EQUATIONS
    # ───────────────────────────────────────────────────────────────────────────
    'equations': {
        'gw_timescale': 'τ_GW = (5/256) × (c⁵a⁴) / (G³M₁M₂(M₁+M₂))',
        'scm_ua_extraction': 'dE/dt = Ug4 × V_interaction × (1 + f_TRZ) / τ_cross',
        'vacuum_energy_density': 'ρ_vac,X = Σf_i E_i / V_object',
        'ug4_full': 'Ug4 = k_4 ρ_vac,[SCm] (M_bh/d_g) e^(-αt) cos(ωt_n) (1 + f_feedback)',
        'modified_timescale': '1/τ_mod = 1/τ_GW + 1/τ_[SCm]-[UA]',
        'metal_retention': 'f_Z = (M_Z,disk_gas + M_Z,disk_stars) / M_Z,formed',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS IN CondensedPhysics.py
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'FinalParsecProblemModel',
        'CGMMetalRetentionModel',
        'SMBHInteractionModel',
        'GravitationalWaveModel',
        'MergerRateModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# EQUATIONS OF THE ATOM (Document 16) - UQFF Standard Model Framework
# ═══════════════════════════════════════════════════════════════════════════════
EQUATIONS_OF_ATOM_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # DOCUMENT METADATA
    # ───────────────────────────────────────────────────────────────────────────
    'document': 'Compressed Summary of Equations of The Atom in UQFF Framework',
    'document_id': 16,
    'source': 'Star Magic / UQFF Framework',
    'date': '2025-02-09',
    'grok_session': 'https://x.com/i/grok',
    
    # ───────────────────────────────────────────────────────────────────────────
    # FUNDAMENTAL CONSTANTS
    # ───────────────────────────────────────────────────────────────────────────
    'c_light': 2.998e8,                # m/s
    'eV_to_J': 1.602e-19,              # J/eV
    'MeV_to_J': 1.602e-13,             # J/MeV
    'GeV_to_J': 1.602e-10,             # J/GeV
    'h_bar': 1.055e-34,                # J·s
    'alpha_em': 1.0/137.036,           # Fine structure constant
    'alpha_s_MZ': 0.118,               # Strong coupling at M_Z
    
    # ───────────────────────────────────────────────────────────────────────────
    # UQFF BASE PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    'E_0': 1e-20,                      # Base energy (J) for 26-level
    'SSq': 0.57,                       # [SSq] calibrated value
    'gamma_decay': 0.0005,             # κ decay constant (per day)
    'rho_vac_SCm': 7.09e-37,           # J/m³ - [SCm] vacuum density
    'rho_vac_UA': 7.09e-36,            # J/m³ - [UA] vacuum density
    'f_TRZ': 0.1,                      # TRZ boost factor
    
    # ───────────────────────────────────────────────────────────────────────────
    # 26-LEVEL ENERGY STRUCTURE
    # ───────────────────────────────────────────────────────────────────────────
    '26_level_formula': 'E_n = E_0 × 10^n (E_0 = 10^{-20} J)',
    '26_level_mapping': {
        'n_1_5': 'Sub-quantum vortices (quarks/gluons)',
        'n_6_10': 'Leptons/quarks',
        'n_11_15': 'Kaons/mesons (plasma decay)',
        'n_16_20': 'Bosons/Higgs (mass generation)',
        'n_21_26': 'Cosmic (TRZs, non-local jumps)',
    },
    'level_energies_J': {
        1: 1e-19, 2: 1e-18, 3: 1e-17, 4: 1e-16, 5: 1e-15,
        6: 1e-14, 7: 1e-13, 8: 1e-12, 9: 1e-11, 10: 1e-10,
        11: 1e-9, 12: 1e-8, 13: 1e-7, 14: 1e-6, 15: 1e-5,
        16: 1e-4, 17: 1e-3, 18: 1e-2, 19: 1e-1, 20: 1.0,
        21: 1e1, 22: 1e2, 23: 1e3, 24: 1e4, 25: 1e5, 26: 1e6,
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # PDG 2025 PARTICLE DATA (verified values)
    # ───────────────────────────────────────────────────────────────────────────
    'pdg_2025': {
        'u_quark': {'mass_MeV': 2.16, 'n': 1, 'charge': 2/3, 'spin': 0.5},
        'd_quark': {'mass_MeV': 4.67, 'n': 2, 'charge': -1/3, 'spin': 0.5},
        's_quark': {'mass_MeV': 93.0, 'n': 5, 'charge': -1/3, 'spin': 0.5},
        'c_quark': {'mass_GeV': 1.27, 'n': 9, 'charge': 2/3, 'spin': 0.5},
        'b_quark': {'mass_GeV': 4.18, 'n': 10, 'charge': -1/3, 'spin': 0.5},
        't_quark': {'mass_GeV': 172.76, 'n': 14, 'charge': 2/3, 'spin': 0.5},
        'electron': {'mass_MeV': 0.511, 'n': 6, 'charge': -1, 'spin': 0.5},
        'muon': {'mass_MeV': 105.66, 'n': 8, 'charge': -1, 'spin': 0.5},
        'tau': {'mass_MeV': 1776.9, 'n': 9, 'charge': -1, 'spin': 0.5},
        'nu_e': {'mass_eV': 0.0, 'n': 1, 'charge': 0, 'spin': 0.5},
        'nu_mu': {'mass_eV': 0.0, 'n': 1, 'charge': 0, 'spin': 0.5},
        'nu_tau': {'mass_eV': 0.0, 'n': 1, 'charge': 0, 'spin': 0.5},
        'photon': {'mass_eV': 0.0, 'n': 0, 'charge': 0, 'spin': 1},
        'gluon': {'mass_eV': 0.0, 'n': 0, 'charge': 0, 'spin': 1},
        'W_boson': {'mass_GeV': 80.377, 'n': 16, 'charge': 1, 'spin': 1},
        'Z_boson': {'mass_GeV': 91.1876, 'n': 17, 'charge': 0, 'spin': 1},
        'Higgs': {'mass_GeV': 125.25, 'n': 18, 'charge': 0, 'spin': 0},
        'proton': {'mass_MeV': 938.272, 'n': 7, 'charge': 1, 'spin': 0.5},
        'neutron': {'mass_MeV': 939.565, 'n': 7, 'charge': 0, 'spin': 0.5},
        'K_plus': {'mass_MeV': 493.677, 'n': 12, 'charge': 1, 'spin': 0},
        'K_minus': {'mass_MeV': 493.677, 'n': 12, 'charge': -1, 'spin': 0},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # DOCUMENT 16 UQFF ENERGY SOLUTIONS (Table 1)
    # ───────────────────────────────────────────────────────────────────────────
    'uqff_energies': {
        'u_quark': {'E_J': 3.68e-13, 'n': 1, 'mass_MeV': 2.16},
        'd_quark': {'E_J': 7.94e-13, 'n': 2, 'mass_MeV': 4.67},
        'electron': {'E_J': 8.19e-14, 'n': 6, 'mass_MeV': 0.511},
        'muon': {'E_J': 1.69e-11, 'n': 8, 'mass_MeV': 105.66},
        'tau': {'E_J': 2.85e-10, 'n': 9, 'mass_MeV': 1776.9},
        'K_plus': {'E_J': 7.91e-11, 'n': 12, 'mass_MeV': 493.7},
        'proton': {'E_J': 1.50e-10, 'n': 7, 'mass_MeV': 938.3},
        'neutron': {'E_J': 1.51e-10, 'n': 7, 'mass_MeV': 939.6},
        'W_boson': {'E_J': 1.29e-8, 'n': 16, 'mass_GeV': 80.4},
        'Z_boson': {'E_J': 1.46e-8, 'n': 17, 'mass_GeV': 91.2},
        'Higgs': {'E_J': 2.00e-8, 'n': 18, 'mass_GeV': 125.25},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # KEY EQUATIONS (Document 16)
    # ───────────────────────────────────────────────────────────────────────────
    'equations': {
        'particle_energy_density': 'ρ_particle = λ × ρ_vac,[SCm] × ρ_vac,[UA] × ω(t) × cos(πt_n) × (1+f_TRZ) × exp(-[SSq]^{n/26} × exp(-π-t))',
        'gluon_field_tensor': 'G_μν = α_s × (ρ_vac,[UA] / r) × exp(-γt)',
        'jump_probability': 'P_jump = 1 - exp(-λ_g × r)',
        'ponderomotive_force': 'F_p = -(e² / 4mω²) × ∇(E²)',
        'negative_time_operator': 't⁻ = -t_n × exp(π - t_n)',
        'nonlocal_jump': '[SSq]^{n/26} × exp(-π - t)',
        '26_level_energy': 'E_n = E_0 × 10^n',
        'uqff_energy': 'E_UQFF = m × c² × exp(n/26)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # GLUON FIELD PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    'gluon_params': {
        'lambda_g': 1e15,              # Jump coupling constant (1/m)
        'r_nuclear': 1e-15,            # Nuclear scale (1 fm)
        'color_charges': ['r', 'g', 'b', 'r̄', 'ḡ', 'b̄'],
        'num_gluons': 8,
        'confinement_scale_fm': 1.0,
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # PARTICLE CLASSIFICATIONS AS [UA] VORTICES
    # ───────────────────────────────────────────────────────────────────────────
    'particle_classifications': {
        'quarks': 'Sub-quantum vortices in [UA] field',
        'gluons': 'Color flux tubes mediating strong force through [UA]',
        'leptons': 'Stable vortices (electron) or metastable (muon, tau)',
        'kaons': 'Meson vortex pairs with s-quark plasma decay modes',
        'bosons': 'Mass generation through [UA]-Higgs coupling',
        'proton_neutron': 'Bound quark vortex triplets (uud, udd)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS IN CondensedPhysics.py
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'EquationsOfTheAtomModel',
        'StandardModelUQFFModel',
        'AtomicModelUQFF',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# UPDATED COMPRESSED SUMMARY - UQFF/STAR MAGIC (Document 17)
# ═══════════════════════════════════════════════════════════════════════════════
UQFF_UPDATED_SUMMARY_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # DOCUMENT METADATA
    # ───────────────────────────────────────────────────────────────────────────
    'document': 'Updated Compressed Summary of UQFF/Star Magic Framework',
    'document_id': 17,
    'source': 'Star Magic Framework (~7000 pages compressed)',
    'date': '2025-02-09',
    'grok_session': 'https://x.com/i/grok',
    
    # ───────────────────────────────────────────────────────────────────────────
    # FUNDAMENTAL CONSTANTS
    # ───────────────────────────────────────────────────────────────────────────
    'c_light': 2.998e8,                # m/s
    'G': 6.674e-11,                    # m³/(kg·s²)
    'h_bar': 1.055e-34,                # J·s
    'M_sun': 1.989e30,                 # kg
    'r_sun': 6.957e8,                  # m
    
    # ───────────────────────────────────────────────────────────────────────────
    # UQFF COUPLING CONSTANTS
    # ───────────────────────────────────────────────────────────────────────────
    'k_1': 1.5,                        # Ug1 coupling (internal dipole)
    'k_2': 1.2,                        # Ug2 coupling (outer field)
    'k_3': 1.8,                        # Ug3 coupling (magnetic disk)
    'k_4': 1.0,                        # Ug4 coupling (galactic SMBH)
    'beta_i': 0.6,                     # Buoyancy coupling (60% opposition)
    'epsilon_sw': 0.001,               # Solar wind modulation
    'delta_sw': 0.01,                  # Solar wind velocity modulation
    'eta': 1e-22,                      # Aether coupling (fine-structure-like)
    'gamma': 0.0001,                   # String decay constant (1/day)
    'alpha_kappa': 0.0005,             # κ - E_react decay rate (1/day)
    
    # ───────────────────────────────────────────────────────────────────────────
    # GALACTIC PARAMETERS (UQFF PREDICTIONS)
    # ───────────────────────────────────────────────────────────────────────────
    'M_bh_UQFF': 8.15e36,              # kg (UQFF prediction ~4.1e6 M_☉)
    'd_g_UQFF': 2.55e20,               # m (UQFF prediction ~27,000 ly)
    'omega_g_UQFF': 7.3e-16,           # rad/s (UQFF prediction)
    
    # ───────────────────────────────────────────────────────────────────────────
    # 2025 OBSERVED VALUES
    # ───────────────────────────────────────────────────────────────────────────
    'M_bh_2025': 8.55e36,              # kg (4.3e6 M_☉, GRAVITY/Keck)
    'd_g_2025': 2.44e20,               # m (25,800 ly, Gaia DR4)
    'omega_g_2025': 9.5e-16,           # rad/s (233 km/s at 25.8 kly)
    'rho_sw_2025': 8e-21,              # kg/m³ (Parker Solar Probe)
    
    # ───────────────────────────────────────────────────────────────────────────
    # SUN REFERENCE FIELD VALUES (t=0, t_n=0)
    # ───────────────────────────────────────────────────────────────────────────
    'Ug1_Sun': 1.39e26,                # J/m³ (internal dipole)
    'Ug2_Sun': 1.18e53,                # J/m³ (outer field)
    'Ug3_Sun': 1.8e49,                 # J/m³ (magnetic disk)
    'Ug4_Sun': 2.5e-20,                # J/m³ (galactic influence)
    'Ubi_Sun': -1.94e27,               # J/m³ (buoyancy opposition)
    'Um_Sun': 2.26e16,                 # J/m³ (lossless strings)
    
    # ───────────────────────────────────────────────────────────────────────────
    # 26-LEVEL POLYNOMIAL STRUCTURE
    # ───────────────────────────────────────────────────────────────────────────
    'E_0': 1e-20,                      # J (vacuum base energy)
    '26_level_formula': 'E_n = E_0 × 10^n (E_0 = 10^{-20} J)',
    '26_level_ranges': {
        'n_1_5': {'E_range': '10^{-19}-10^{-15} J', 'application': 'Sub-quantum ([UA] vortices)'},
        'n_6_10': {'E_range': '10^{-14}-10^{-10} J', 'application': 'Nuclear (PDG bindings)'},
        'n_11_15': {'E_range': '10^{-9}-10^{-5} J', 'application': 'Plasma/molecular (solar wind)'},
        'n_16_20': {'E_range': '10^{-4}-1 J', 'application': 'Higgs/stellar (LHC m_H=125 GeV)'},
        'n_21_26': {'E_range': '10-10^6 J', 'application': 'Galactic (Fermi jets, DPM inflation)'},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # VACUUM DENSITIES
    # ───────────────────────────────────────────────────────────────────────────
    'rho_vac_UA': 7.09e-37,            # J/m³ ([UA] vacuum density)
    'rho_vac_SCm': 7.09e-37,           # J/m³ ([SCm] vacuum density)
    'rho_vac_sw': 8e-21,               # J/m³ (solar wind vacuum)
    
    # ───────────────────────────────────────────────────────────────────────────
    # E_REACT (QUASAR/CORE OUTPUT)
    # ───────────────────────────────────────────────────────────────────────────
    'E_react_0': 1e46,                 # J (initial reactor output)
    'E_react_formula': 'E_react = 10^{46} × e^{-0.0005t}',
    
    # ───────────────────────────────────────────────────────────────────────────
    # KEY EQUATIONS (F_U COMPONENTS)
    # ───────────────────────────────────────────────────────────────────────────
    'equations': {
        'F_U_complete': 'F_U = Σ_i [k_i U_gi - β_i U_gi ω_g M_bh / d_g E_react] + Σ_j [μ_j/r_j (1-e^{-γt cos(πt_n)}) ϕ_j] + g_μν + η T_s^{μν} - Σ_i [δ_i U_i E_react]',
        'Ug1': 'U_g1 = k_1 μ_s (M_s/r) e^{-αt} cos(πt_n) (1+β_def)',
        'Ug2': 'U_g2 = k_2 (λ_vac,[UA]+λ_vac,[SCm]) M_s/r² S(r-R_b) (1+δ_sw v_sw) H_SCm E_react',
        'Ug3': 'U_g3 = k_3 Σ_j B_j cos(ω_s t) P_core E_react',
        'Ug4': 'U_g4 = k_4 λ_vac,[SCm] M_bh/d_g e^{-αt} cos(πt_n) (1+f_feedback)',
        'Ubi': 'U_bi = -β_i U_gi ω_g M_bh/d_g (1+ε_sw ρ_vac,sw) [UA] cos(πt_n)',
        'Um': 'Um = Σ_j [μ_j/r_j (1-e^{-γt cos(πt_n)}) ϕ_j] P_SCm E_react',
        'A_mu_nu': 'A_μν = g_μν + η T_s^{μν}',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # DPM (DI-PSEUDO-MONOPOLE) PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    'dpm': {
        'description': 'DPM birth from pre-Big Bang [SCm]-[UA] reaction in 26-shell EM field',
        'epochs': {
            't_1': 'Fissile formation',
            't_2': 'Nuclear synthesis',
            't_3': 'Stellar formation',
            't_4': 'Galaxy assembly',
            't_5': 'Globular clusters',
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # VARIABLES TABLE (11 KEY PARAMETERS)
    # ───────────────────────────────────────────────────────────────────────────
    'variables_table': {
        'eta': {'value': 1e-22, 'units': 'dimensionless', 'description': 'Aether coupling'},
        'g_mu_nu': {'value': '[1,-1,-1,-1]', 'units': 'dimensionless', 'description': 'Minkowski metric'},
        'beta_i': {'value': 0.6, 'units': 'dimensionless', 'description': 'Buoyancy coupling (60%)'},
        'epsilon_sw': {'value': 0.001, 'units': 'dimensionless', 'description': 'Solar wind modulation'},
        'k_i': {'value': '[1.5,1.2,1.8,1.0]', 'units': 'dimensionless', 'description': 'Ug couplings'},
        'r_j': {'value': 1.496e13, 'units': 'm', 'description': 'String distance (100 AU)'},
        'd_g': {'value': 2.55e20, 'units': 'm', 'description': 'Sun-Sgr A* distance'},
        'f_feedback': {'value': 0.1, 'units': 'per dex', 'description': 'BH feedback factor'},
        'omega_g': {'value': 7.3e-16, 'units': 'rad/s', 'description': 'Galactic spin rate'},
        'M_bh': {'value': 8.15e36, 'units': 'kg', 'description': 'Sgr A* SMBH mass'},
        'kappa': {'value': 0.0005, 'units': '1/day', 'description': 'E_react decay rate'},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # VERIFICATION ERRORS (2025 DATA)
    # ───────────────────────────────────────────────────────────────────────────
    'verification_errors': {
        'd_g': {'error_percent': 4.5, 'source': 'Gaia DR4'},
        'M_bh': {'error_percent': 4.7, 'source': 'GRAVITY/Keck'},
        'omega_g': {'error_percent': 23, 'source': 'Rotation curves'},
        'rho_sw': {'error_percent': 0, 'source': 'Parker Solar Probe'},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # MILLENNIUM PROBLEMS CONNECTION
    # ───────────────────────────────────────────────────────────────────────────
    'millennium_connections': {
        'navier_stokes': 'Flux tube jets in [UA] plasma',
        'yang_mills': 'Mass gap via [SCm] Meissner effect',
        'riemann': 'π cycles in cos(πt_n) oscillations',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS IN CondensedPhysics.py
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'UQFF2025VerificationSummary',
        'UnifiedFieldEquation',
        'CompleteUnifiedFieldModel',
        'UQFFCalibrationSummaryModel',
        'DPMModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 18: VARIABLE EXPLANATIONS REFINEMENT
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_VARIABLE_EXPLANATIONS_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # HEAVISIDE FRACTION (f_Heaviside)
    # ───────────────────────────────────────────────────────────────────────────
    'f_Heaviside': 0.01,                    # Heaviside fraction
    'f_Heaviside_amplification': 1e11,      # 10¹³ × 0.01 = 10¹¹ amplification
    'Um_amplified': 2.28e65,                # 2.28×10⁵⁴ × 10¹¹ ≈ 2.28×10⁶⁵ J/m³
    'Um_raw': 2.28e54,                      # Before Heaviside amplification
    
    # ───────────────────────────────────────────────────────────────────────────
    # HELIOSPHERE THICKNESS FACTOR (H_SCm)
    # ───────────────────────────────────────────────────────────────────────────
    'H_SCm': 0.99,                          # Heliosphere thickness factor ~1
    'H_SCm_min': 0.9,                       # Minimum value
    'H_SCm_max': 1.1,                       # Maximum value
    'heliosphere_radius_AU': 100,           # Heliosphere at ~100 AU
    'heliosphere_boundary_thickness_AU': 0.01,  # ~0.01 AU thick boundary
    'scales': 'Ug2 (charge-reactivity)',    # H_SCm scales Ug2 component
    
    # ───────────────────────────────────────────────────────────────────────────
    # UG INDEX (i = 1-4)
    # ───────────────────────────────────────────────────────────────────────────
    'Ug_indices': {
        1: 'Ug1 - Magnetic dipole term',
        2: 'Ug2 - Charge-reactivity term',
        3: 'Ug3 - String rotation term',
        4: 'Ug4 - Vacuum concentration term',
    },
    'Ug_i_count': 4,                        # Four Ug components
    
    # ───────────────────────────────────────────────────────────────────────────
    # STRING INDEX (j = billions/trillions)
    # ───────────────────────────────────────────────────────────────────────────
    'j_min': 1,                             # First string index
    'j_max_conceptual': 1e12,               # Trillions of strings (conceptual)
    'n_layers': 26,                         # 26 quantum layers for computation
    'Σ_j_interpretation': 'Sum over all string contributions',
    
    # ───────────────────────────────────────────────────────────────────────────
    # INERTIA COUPLING (λ_i)
    # ───────────────────────────────────────────────────────────────────────────
    'lambda_i': 1.0,                        # Uniform inertia coupling
    'lambda_i_interpretation': 'λ_i = 1.0 → uniform coupling across all i',
    
    # ───────────────────────────────────────────────────────────────────────────
    # TIME-VARYING DIPOLE MOMENT (μ_j)
    # ───────────────────────────────────────────────────────────────────────────
    'mu_j_formula': '(10³ + 0.4 sin(ω_c t)) × 3.38×10²⁰ T·m³',
    'mu_j_base': 1e3,                       # Base coefficient (unitless)
    'mu_j_amplitude': 0.4,                  # Oscillation amplitude
    'mu_j_scale': 3.38e20,                  # Scale factor T·m³
    'mu_j_at_t0': 3.38e23,                  # μ_j at t=0: 10³ × 3.38×10²⁰
    'omega_c': 1.587e-8,                    # Solar cycle frequency rad/s
    'solar_cycle_years': 12.55,             # ~12.55 year solar cycle
    
    # ───────────────────────────────────────────────────────────────────────────
    # INERTIAL TERM (U_i)
    # ───────────────────────────────────────────────────────────────────────────
    'U_i_formula': 'U_i ≈ 1.38×10⁻⁴⁷ × λ_i',
    'U_i_value': 1.38e-47,                  # Inertial term value
    'U_i_times_lambda': 1.38e-47,           # U_i × λ_i (with λ_i=1)
    'U_i_components': ['rho_vac_SCm', 'rho_vac_UA', 'omega_s', 'cos_tn', 'f_TRZ'],
    
    # ───────────────────────────────────────────────────────────────────────────
    # NEGATIVE TIME (t_n)
    # ───────────────────────────────────────────────────────────────────────────
    't_n_formula': 't_n = t - t_0',
    't_n_negative_allowed': True,           # t_n < 0 permitted
    'cos_pi_tn_purpose': 'Enables time reversals in field equations',
    't_n_oscillation_period_days': 2,       # cos(πt_n) period ~2 days in t_n
    't_minus_formula': 't⁻ = -t_n × e^(π - t_n)',
    
    # ───────────────────────────────────────────────────────────────────────────
    # GALACTIC SPIN (Ω_g)
    # ───────────────────────────────────────────────────────────────────────────
    'Omega_g_observed': 7.3e-16,            # Observed galactic spin (rad/s)
    'Omega_g_calculated': 8.9e-16,          # Calculated from v=220 km/s at r=8 kpc
    'Omega_g_discrepancy_percent': 22,      # ~22% discrepancy (DM drag?)
    'v_galactic': 220e3,                    # 220 km/s in m/s
    'r_galactic_kpc': 8,                    # 8 kpc
    'r_galactic_m': 2.47e20,                # 8 kpc in meters
    'Omega_g_inner_pattern': 4e-16,         # Inner pattern speed (rad/s)
    
    # ───────────────────────────────────────────────────────────────────────────
    # TIME FACTOR
    # ───────────────────────────────────────────────────────────────────────────
    'time_factor_formula': '1 - e^(-γt × cos(πt_n))',
    'gamma': 5e-5,                          # Decay rate κ = 0.0005/day
    'time_factor_behavior': {
        'cos_tn_positive': 'Exponent negative → DECAY',
        'cos_tn_negative': 'Exponent positive → GROWTH',
        't_0': 'Factor = 0',
        't_infinity': 'Factor → 1 (for t_n near 0)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # SGR A* BLACK HOLE REFERENCE
    # ───────────────────────────────────────────────────────────────────────────
    'M_bh_SgrA': 4.3e6,                     # Sgr A* mass in M_☉
    'M_bh_SgrA_kg': 8.55e36,                # Sgr A* mass in kg
    'd_g': 2.47e20,                         # 8 kpc in meters
    'BH_spin_effect': 'Squashes spacetime via frame dragging',
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS IN CondensedPhysics.py
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'UniversalMagnetismModel',          # Um with Heaviside factor
        'UniversalInertiaModel',            # U_i with λ_i, μ_j, t_n
        'DPMModel',                         # Contains H_SCm, Ω_g
        'UQFF2025VerificationSummary',      # Core F_U verification
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 19: UQFF PARAMETER REFINEMENTS
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_PARAMETER_REFINEMENTS_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # CORE PENETRATION FACTOR (P_core)
    # ───────────────────────────────────────────────────────────────────────────
    'P_core_Sun': 1.0,                      # Sun: full penetration (plasma core)
    'P_core_planet': 1e-3,                  # Planets: reduced (solid/liquid cores)
    'P_core_scales': 'Ug3',                 # Scales Ug3 string rotation term
    
    # ───────────────────────────────────────────────────────────────────────────
    # QUASI-LONGITUDINAL WAVE FACTOR (f_quasi)
    # ───────────────────────────────────────────────────────────────────────────
    'f_quasi': 0.01,                        # +1% to Um wave contribution
    'f_quasi_effect': '(1 + f_quasi) = 1.01',
    
    # ───────────────────────────────────────────────────────────────────────────
    # BUBBLE RADIUS (R_b) - HELIOSPHERE BOUNDARY
    # ───────────────────────────────────────────────────────────────────────────
    'R_b_m': 1.496e13,                      # 100 AU in meters
    'R_b_AU': 100,                          # Heliosphere boundary
    'R_b_observed_range_AU': [100, 122],    # Observed range
    'R_b_scales': 'Ug2',                    # Step function in Ug2
    'R_b_step_function': 'S(r - R_b)',      # S=1 outside, S=0 inside
    
    # ───────────────────────────────────────────────────────────────────────────
    # EFFICIENCY DECAY (E_react with κ)
    # ───────────────────────────────────────────────────────────────────────────
    'E_react_base': 1e46,                   # Base efficiency (unitless)
    'kappa': 0.0005,                        # κ: decay rate (day⁻¹)
    'tau_kappa_days': 2000,                 # 1/κ = 2000 days
    'tau_kappa_years': 5.5,                 # ~5.5 year timescale
    'E_react_formula': 'E_react = 10⁴⁶ × e^(-0.0005t)',
    
    # ───────────────────────────────────────────────────────────────────────────
    # RECIPROCATION DECAY (γ for Um)
    # ───────────────────────────────────────────────────────────────────────────
    'gamma': 5e-5,                          # γ: Um decay rate (day⁻¹)
    'tau_gamma_days': 20000,                # 1/γ = 20000 days
    'tau_gamma_years': 55,                  # ~55 year timescale
    'gamma_in_formula': '1 - e^(-γt × cos(πt_n))',
    
    # ───────────────────────────────────────────────────────────────────────────
    # SCM PENETRATION FACTOR (P_SCm)
    # ───────────────────────────────────────────────────────────────────────────
    'P_SCm_Sun': 1.0,                       # Sun: full [SCm] penetration
    'P_SCm_planet': 1e-3,                   # Planets: reduced
    'P_SCm_scales': 'Um',                   # Scales Universal Magnetism
    
    # ───────────────────────────────────────────────────────────────────────────
    # SOLAR CYCLE FREQUENCY (ω_c)
    # ───────────────────────────────────────────────────────────────────────────
    'omega_c': 1.587e-8,                    # rad/s
    'T_solar_cycle_s': 3.96e8,              # Period in seconds
    'T_solar_cycle_years': 12.55,           # UQFF value
    'T_solar_cycle_observed': 11,           # Observed average ~11 yr
    'omega_c_in': 'μ_j(t) oscillation',     # Used in time-varying dipole
    
    # ───────────────────────────────────────────────────────────────────────────
    # SOLAR WIND MODULATION (δ_sw)
    # ───────────────────────────────────────────────────────────────────────────
    'delta_sw': 0.01,                       # Solar wind modulation factor
    'delta_sw_scales': 'Ug2',               # Mods Ug2 by v_sw
    'delta_sw_effect': '1 + δ_sw × v_sw = 5001',
    
    # ───────────────────────────────────────────────────────────────────────────
    # SOLAR WIND VELOCITY (v_sw)
    # ───────────────────────────────────────────────────────────────────────────
    'v_sw': 5e5,                            # 500 km/s in m/s
    'v_sw_km_s': 500,                       # km/s
    'v_sw_observed_range': [400, 500],      # Observed avg at 1 AU (km/s)
    
    # ───────────────────────────────────────────────────────────────────────────
    # SUN REFERENCE VALUES (t=0, t_n=0)
    # ───────────────────────────────────────────────────────────────────────────
    'Ug1_Sun': 1.39e26,                     # Internal dipole (J/m³)
    'Ug2_Sun': 1.18e53,                     # Outer field bubble (J/m³)
    'Ug3_Sun': 1.8e49,                      # Magnetic strings disk (J/m³)
    'Ug4_Sun': 2.5e-20,                     # Star-BH interactions (J/m³)
    'Ubi_Sun': -1.94e27,                    # Buoyancy opposition (J/m³)
    'Um_Sun': 2.28e65,                      # Universal Magnetism (J/m³)
    'A_mu_nu_Sun': [1, -1, -1, -1],         # Metric tensor + 1e-15
    'UI_Sun': 1.38e-47,                     # Universal Inertia (J/m³)
    
    # ───────────────────────────────────────────────────────────────────────────
    # F_U COMPLETE EQUATION
    # ───────────────────────────────────────────────────────────────────────────
    'F_U_formula': 'F_U = Σ_i[k_i Ug_i - β_i Ug_i Ω_g M_bh/d_g E_react] + Um + g_μν + ηT_s^μν - Σ_i[λ_i U_i E_react]',
    'Um_formula': 'Um = Σ_j[μ_j/r_j(1-e^{-γt cos(πt_n)})φ̂_j] × P_SCm × E_react × (1+10¹³f_H)(1+f_quasi)',
    
    # ───────────────────────────────────────────────────────────────────────────
    # 26-LEVEL POLYNOMIAL (unchanged)
    # ───────────────────────────────────────────────────────────────────────────
    'n_levels': 26,
    'energy_levels': {
        '1-5': 'Sub-quantum (10⁻¹⁹ to 10⁻¹⁵ J)',
        '6-10': 'Nuclear (10⁻¹⁴ to 10⁻¹⁰ J)',
        '11-15': 'Plasma/solar wind (10⁻⁹ to 10⁻⁵ J, heliosphere n=13)',
        '16-20': 'Higgs/stellar (10⁻⁴ to 1 J)',
        '21-26': 'Galactic (10 to 10⁶ J)',
    },
    'n_heliosphere': 13,                    # n=13 for heliosphere/solar wind
    
    # ───────────────────────────────────────────────────────────────────────────
    # VERIFICATION STATUS
    # ───────────────────────────────────────────────────────────────────────────
    'verification': {
        'M_bh_SgrA': {'UQFF': 8.15e36, 'observed': 8.55e36, 'error_percent': 5},
        'Omega_g': {'UQFF': 7.3e-16, 'calculated': 8.9e-16, 'error_percent': 22},
        'R_b': {'UQFF': 100, 'observed_range': [100, 122], 'status': 'ALIGNED'},
        'solar_cycle': {'UQFF': 12.55, 'observed': 11, 'status': 'CLOSE'},
        'v_sw': {'UQFF': 500, 'observed_range': [400, 500], 'status': 'ALIGNED'},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'UnifiedFieldEquation',
        'UniversalMagnetismModel',
        'UniversalInertiaModel',
        'DPMModel',
        'UQFF2025VerificationSummary',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 20: UQFF SOLAR/STELLAR VARIABLES
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_SOLAR_STELLAR_VARIABLES_PARAMS: Dict[str, Any] = {
    # ───────────────────────────────────────────────────────────────────────────
    # STELLAR MASS (M_s)
    # ───────────────────────────────────────────────────────────────────────────
    'M_s_Sun': 1.989e30,                    # Solar mass (kg)
    'M_s_units': 'kg',
    'M_s_in_equations': ['Ug1 (∇M_s/r)', 'Ug2 (M_s/r²)'],
    
    # ───────────────────────────────────────────────────────────────────────────
    # ROTATION RATE (ω_s)
    # ───────────────────────────────────────────────────────────────────────────
    'omega_s_Sun': 2.5e-6,                  # rad/s (~29 day period)
    'omega_s_equator_calc': 2.83e-6,        # Calculated from 25.67 day equator
    'omega_s_period_days': 29,              # UQFF value
    'omega_s_equator_days': 25.67,          # Observed equator rotation
    'omega_s_differential': [25, 33],       # Differential rotation range (days)
    
    # ───────────────────────────────────────────────────────────────────────────
    # HEAVISIDE STEP FUNCTION (S)
    # ───────────────────────────────────────────────────────────────────────────
    'S_definition': 'S(r - R_b) = 1 if r > R_b, 0 otherwise',
    'S_purpose': 'Boundary activation for Ug2',
    'R_b_AU': 100,                          # Heliosphere boundary
    
    # ───────────────────────────────────────────────────────────────────────────
    # STRESS-ENERGY TENSOR (T_s^μν)
    # ───────────────────────────────────────────────────────────────────────────
    'T_s_mu_nu': 1.123e7,                   # Speculative stress-energy (J/m³)
    'T_s_mu_nu_SCm': 1.11e7,                # [SCm] contribution
    'T_s_mu_nu_UA': 1.27e3,                 # [UA] contribution
    'eta_coupling': 1e-22,                  # η coupling to metric
    'A_mu_nu_formula': 'A_μν = g_μν + η × T_s^μν',
    'eta_T_perturbation': 1.123e-15,        # η × T_s^μν perturbation
    
    # ───────────────────────────────────────────────────────────────────────────
    # SURFACE MAGNETIC FIELD (B_s)
    # ───────────────────────────────────────────────────────────────────────────
    'B_s_quiet_Sun': 1e-4,                  # Quiet Sun (T) = 1 G
    'B_s_sunspot_max': 0.4,                 # Sunspots max (T) = 4000 G
    'B_s_range': [1e-4, 0.4],               # Full range (T)
    'B_s_units': 'T',
    
    # ───────────────────────────────────────────────────────────────────────────
    # SURFACE TEMPERATURE (T_s)
    # ───────────────────────────────────────────────────────────────────────────
    'T_s_Sun': 5778,                        # Effective photospheric (K)
    'T_s_observed_range': [5772, 5800],     # Observed range (K)
    'T_s_core': 1.5e7,                      # Core temperature (K)
    
    # ───────────────────────────────────────────────────────────────────────────
    # TIME-REVERSAL ZONE FACTOR (f_TRZ)
    # ───────────────────────────────────────────────────────────────────────────
    'f_TRZ': 0.1,                           # +10% to Ui for negentropy
    'f_TRZ_effect': '(1 + f_TRZ) = 1.1',
    'f_TRZ_interpretation': 'Negentropy contribution in TRZs',
    'f_TRZ_status': 'SPECULATIVE',
    
    # ───────────────────────────────────────────────────────────────────────────
    # UG1 DEFECT OSCILLATION (δ_def)
    # ───────────────────────────────────────────────────────────────────────────
    'delta_def_amplitude': 0.01,            # ±1% to Ug1
    'delta_def_frequency': 0.001,           # Angular frequency (day⁻¹)
    'delta_def_formula': 'δ_def = 0.01 × sin(0.001 × t)',
    'delta_def_period_years': 17.2,         # 2π/0.001 ÷ 365.25
    'delta_def_status': 'SPECULATIVE',
    
    # ───────────────────────────────────────────────────────────────────────────
    # DISK UNIT VECTOR (φ̂_j)
    # ───────────────────────────────────────────────────────────────────────────
    'phi_hat_j': 1.0,                       # Unit vector magnitude
    'phi_hat_interpretation': 'Direction of j-th string contribution in Ug3 disk',
    'phi_hat_units': 'Unitless (normalized)',
    
    # ───────────────────────────────────────────────────────────────────────────
    # SUN REFERENCE VALUES (t=0, t_n=0, r>R_b)
    # ───────────────────────────────────────────────────────────────────────────
    'sun_reference': {
        'conditions': 't=0, t_n=0, r>R_b, S=1',
        'Ug1': '1.39e26 × (M_s/r) × (1 + δ_def=0)',
        'Ug2': '1.18e53 × S=1 × (M_s/r²) × (1 + δ_sw × v_sw)',
        'Ug3': '1.8e49 × P_core × cos(ω_s × t × π)',
        'Ug4': '2.5e-20 × (M_bh/d_g)',
        'Ubi': '-1.94e27',
        'Um': '2.28e65 with φ̂_j ~1',
        'A_mu_nu': '[1,-1,-1,-1] + η × T_s^μν ~1e-15',
        'Ui': '1.38e-47 × (1 + f_TRZ=0.1)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'UnifiedFieldEquation',
        'UniversalMagnetismModel',
        'UniversalInertiaModel',
        'DPMModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 21: VACUUM DENSITIES & 26-LEVEL POLYNOMIAL FRAMEWORK
# ═══════════════════════════════════════════════════════════════════════════════
# Source: UQFF Proof Set Document 21
# Topic: Vacuum Energy Densities (ρ_vac) and 26-Level Polynomial E_n Hierarchy
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_VAC_DENSITIES_26LEVEL_PARAMS: Dict[str, Any] = {
    'document_id': 21,
    'document_name': 'UQFF Vacuum Densities & 26-Level Polynomial Framework',
    'created': '2026-02-10',
    
    # ───────────────────────────────────────────────────────────────────────────
    # VAC DENSITY HIERARCHY (ρ_vac,X)
    # ───────────────────────────────────────────────────────────────────────────
    # Four vacuum energy density components from [SCm]-[UA] DPM inflation framework
    
    'rho_vac_A': {
        'value': 1e-23,
        'units': 'J/m³',
        'description': 'Universal Cosmic Aether vacuum energy density',
        'status': 'SPECULATIVE',
        'notes': 'No empirical match 2025; denominator in E_react equation',
    },
    
    'rho_vac_Ui': {
        'value': 2.84e-36,
        'units': 'J/m³',
        'description': 'Universal Inertia vacuum density (Sun)',
        'status': 'SPECULATIVE',
        'notes': 'Contributes to Ui; unverified empirically',
    },
    
    'rho_vac_SCm': {
        'value': 7.09e-37,
        'units': 'J/m³',
        'description': '[SCm] Superconductive Material vacuum density (Sun)',
        'status': 'SPECULATIVE',
        'notes': 'Part of [SCm]-[UA] DPM; scales with system',
    },
    
    'rho_vac_UA': {
        'value': 7.09e-36,
        'units': 'J/m³',
        'description': '[UA] Universal Aether vacuum density (Sun)',
        'status': 'SPECULATIVE',
        'notes': 'ρ_vac,[UA]/ρ_vac,[SCm] = 10 at all scales',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # v_SCm PROPAGATION VELOCITY
    # ───────────────────────────────────────────────────────────────────────────
    
    'v_SCm': {
        'value': 1e8,
        'units': 'm/s',
        'description': 'SCm propagation velocity (~c/3)',
        'status': 'SPECULATIVE',
        'notes': 'c/3 unverified; used in E_react numerator',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # E_react REACTOR EFFICIENCY FORMULA
    # ───────────────────────────────────────────────────────────────────────────
    # E_react = (ρ_vac,[SCm] × v_SCm²) / ρ_vac,A × e^(-κt) ≈ 10^46 × e^(-κt)
    
    'E_react_formula': '(ρ_vac,[SCm] × v_SCm²) / ρ_vac,A × e^(-κt)',
    'E_react_initial': 1e46,
    'E_react_decay_rate': 0.0005,  # κ = 0.0005/day
    
    # ───────────────────────────────────────────────────────────────────────────
    # 26-LEVEL POLYNOMIAL ENERGY STRUCTURE
    # ───────────────────────────────────────────────────────────────────────────
    # E_n = E_0 × 10^n (26 quantum levels spanning 25 orders of magnitude)
    
    'level_ranges': {
        'sub_quantum': {
            'levels': [1, 2, 3, 4, 5],
            'E_range_J': [1e-19, 1e-15],
            'description': 'Sub-quantum phenomena',
        },
        'nuclear': {
            'levels': [6, 7, 8, 9, 10],
            'E_range_J': [1e-14, 1e-10],
            'description': 'Nuclear binding/reactions',
        },
        'plasma': {
            'levels': [11, 12, 13, 14, 15],
            'E_range_J': [1e-9, 1e-5],
            'description': 'Plasma/molecular states',
        },
        'higgs_stellar': {
            'levels': [16, 17, 18, 19, 20],
            'E_range_J': [1e-4, 1],
            'description': 'Higgs exotic/stellar scales',
        },
        'galactic': {
            'levels': [21, 22, 23, 24, 25, 26],
            'E_range_J': [10, 1e6],
            'description': 'Galactic/universal scales',
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # UNIFIED FIELD EQUATION F_U COMPONENTS (Sun at t=0, t_n=0)
    # ───────────────────────────────────────────────────────────────────────────
    
    'F_U_components_sun': {
        'Ug1': {'value': 1.39e26, 'formula': '(M_s/r) × (1 + δ_def=0)'},
        'Ug2': {'value': 1.18e53, 'formula': 'S=1 × (M_s/r²) × (1 + δ_sw × v_sw)'},
        'Ug3': {'value': 1.8e49, 'formula': 'P_core × cos(ω_s × t × π)', 'phi_j': 1},
        'Ug4': {'value': 2.5e-20, 'formula': '(M_bh/d_g)'},
        'Ubi': {'value': -1.94e27, 'formula': '-β × Ug4'},
        'Um': {'value': 2.28e65, 'formula': 'φ̂_j ~1 (dominant term)'},
        'A_mu_nu': {'value': 1e-15, 'formula': '[1,-1,-1,-1] + η × T_s^μν'},
        'Ui': {'value': 1.38e-47, 'formula': '(1 + f_TRZ=0.1)'},
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # EMPIRICAL VERIFICATION (2025)
    # ───────────────────────────────────────────────────────────────────────────
    
    'verification_summary': {
        'empirical_confirmed': ['M_s', 'T_s', 'B_s', 'R_b', 'v_sw', 'ω_c'],
        'minor_discrepancies': {
            'omega_s': 'Model 2.5e-6 vs measured 2.83e-6 (~12% low)',
            'solar_cycle': 'Model 12.55 yr vs observed ~11 yr',
        },
        'speculative_unverified': [
            'ρ_vac,A', 'ρ_vac,Ui', 'ρ_vac,[SCm]', 'ρ_vac,[UA]',
            'v_SCm', 'f_TRZ', 'δ_def',
        ],
        'cosmological_context': 'ρ_vac << cosmological Λ (~1e-9 J/m³)',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # FRAMEWORK CONTEXT
    # ───────────────────────────────────────────────────────────────────────────
    
    'framework_notes': {
        'inflation': '[SCm]-[UA] DPM births inflation',
        'unification': 'Unifies Ug_i, Um, Ub_i, Ui',
        '26_level': 'Unchanged polynomial structure from core theory',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'VacuumEnergyDensityModel',
        'UnifiedFieldEquation',
        'UQFF_Triadic',
        'PolynomialEnergyLevelModel',
        'ReactorEfficiencyModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 22: UQFF EQUATIONS ACROSS ASTROPHYSICAL SYSTEMS
# ═══════════════════════════════════════════════════════════════════════════════
# Source: UQFF Equations Across Astrophysical Systems_22Sept2025.docx
# Topic: Magnetar, Fokker-Planck, CRP, GW170817, r-process synthesis
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_ASTROPHYSICAL_SYSTEMS_PARAMS: Dict[str, Any] = {
    'document_id': 22,
    'document_name': 'UQFF Equations Across Astrophysical Systems',
    'created': '2025-09-22',    
    
    # ───────────────────────────────────────────────────────────────────────────
    # 1. MAGNETAR SGR 1745-2900 (8-TERM EQUATION)
    # ───────────────────────────────────────────────────────────────────────────
    
    'magnetar_SGR1745': {
        'equation': 'g_Magnetar(r,t) = (G·M)/r² × (1+H(z)·t)(1-B/B_crit) + G·M_BH/r_BH² + Ug_sum + Λc²/3 + ℏ⟨H⟩/√(ΔxΔp)·2π/t_H + q(v×B) + ρ',
        'M_kg': 1.4 * 1.989e30,  # 1.4 M_sun
        'r_m': 10e3,             # 10 km radius
        'B_T': 1e15,             # 10^15 T surface field
        'B_crit': 4.4e13,        # QED critical field
        'terms': [
            'Newtonian + Hubble expansion + magnetic suppression',
            'Sgr A* SMBH contribution',
            'UQFF 4-layer sum (Ug1+Ug2+Ug3+Ug4)',
            'Cosmological constant Λc²/3',
            'Quantum uncertainty term',
            'Lorentz force q(v × B)',
            'Density contribution ρ',
        ],
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 2. η EFFICIENCY EQUATION
    # ───────────────────────────────────────────────────────────────────────────
    
    'eta_efficiency': {
        'equation': 'η = k_η × exp(-[SSq]n/26) × exp(-(π-t)) × Um / ρ_vac,[UA]',
        'k_eta': 1e-113,
        'SSq': 0.57,
        'n_default': 13,
        'Um': 0.99,
        'rho_vac_UA': 7.09e-36,
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 3. CRP FOKKER-PLANCK STEADY-STATE
    # ───────────────────────────────────────────────────────────────────────────
    
    'fokker_planck': {
        'equation': '∂n/∂t = ∂/∂p[(dp/dt)n] + ∂²/∂p²[Dn] + Q - n/t_esc',
        'steady_state': 'n(p) ~ p^{-2.2} × exp(-p/p_max)',
        'p_max_eV': 1e16,
        'spectral_index': 2.2,
        'chi_sq_mock': 0.05,
        'SED_peak_eV': 1e15,
        'pp_dominant': '<0.1 PeV',
        'LLAGN_flux': 'matches IceCube background',
        'neutrino_outflow': 0.70,
        'neutrino_inflow': 0.30,
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 4. CRP EXTENSION TO F_U
    # ───────────────────────────────────────────────────────────────────────────
    
    'crp_extension': {
        'equation': 'F_U += Σ D_E × ∂²n/∂p² × exp(-γt)',
        'D_E_formula': 'D_E ∝ E^{0.5}',
        'D_E_exponent': 0.5,
        'gamma_day': 0.00005,
        'gamma_units': 'day⁻¹',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 5. KILONOVA R-PROCESS (GW170817)
    # ───────────────────────────────────────────────────────────────────────────
    
    'kilonova_GW170817': {
        'event': 'GW170817 NS-NS merger',
        'M_ej_fraction': 0.40,     # 40% of ejecta at high velocity
        'v_ej_c': 0.1,             # 0.1c ejecta velocity
        'r_process_solar': 0.95,   # 95% matches solar r-process
        'Ye_midplane': 0.1,        # Electron fraction midplane
        'Ye_outflow': 0.2,         # Electron fraction outflow
        'A_predicted': 254,        # Predicted mass number from exp term
        'UFE_contribution': 0.05,  # ~5% toward Ultra-Fe elements
        'Ye_chi_sq_note': 'χ² to solar abundances predicts A=254',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 6. VERIFIED PARAMETERS
    # ───────────────────────────────────────────────────────────────────────────
    
    'verified_parameters': {
        'beta_i': 0.61,
        'chi_sq_fit': 0.05,
        'neutrino_unification': 0.995,  # ~99.5%
        'flux_match': '~2% predict flux from ρ_vac ratios ~10^{-38}',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS
    # ───────────────────────────────────────────────────────────────────────────
    'related_models': [
        'UQFFAstrophysicalSystemsModel',
        'IceCubeNeutrinoModel',
        'FokkerPlanckCRPModel',
        'KilonovaRProcessModel',
        'MagnetarGravityModel',
    ],
}


# ═══════════════════════════════════════════════════════════════════════════════
# DOCUMENT 23: UQFF EQUATIONS EXTRACTION (22Sept2025.docx)
# Mathematical extraction from astrophysical systems review (Documents 1-9)
# ═══════════════════════════════════════════════════════════════════════════════

UQFF_EQUATIONS_EXTRACTION_PARAMS: Dict[str, Any] = {
    'document_id': 23,
    'document_name': 'UQFF Equations Extraction from Astrophysical Systems',
    'source_file': 'UQFF Equations Across Astrophysical Systems_22Sept2025.docx',
    'created': '2025-09-22',
    
    # ───────────────────────────────────────────────────────────────────────────
    # 1. FOKKER-PLANCK TIME-DEPENDENT EQUATION
    # ───────────────────────────────────────────────────────────────────────────
    
    'fokker_planck_pde': {
        'equation': '∂n/∂t = ∂/∂p[(dp/dt)n] + ∂²/∂p²[Dn] + Q - n/t_esc',
        'steady_state': 'n(p) ~ p^{-2.2} × exp(-p/p_max)',
        'type': 'time-dependent partial differential equation',
        'terms': {
            'advection': '∂/∂p[(dp/dt)n]',
            'diffusion': '∂²/∂p²[Dn]',
            'source': 'Q (injection rate)',
            'escape': '-n/t_esc (cosmic ray escape)',
        },
        'parameters': {
            'spectral_index': 2.2,
            'p_max_eV': 1e16,
            't_esc_days': 1e6,
            'D_E_exponent': 0.5,
        },
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # 2. EXTRACTED NUMERICAL VALUES FROM DOCUMENT
    # ───────────────────────────────────────────────────────────────────────────
    
    'extracted_values': {
        'p_max': '~10^{16} eV',
        'n_p_spectral': 'p^{-2.2}',
        'chi_sq_mock': 0.05,
        'SED_peak': '~10^{15} eV',
        'beta_i': 0.61,
        'gamma_day': 0.00005,
        'D_E_exponent': 0.5,
        'M_ej_fraction': 0.40,
        'v_ej_c': 0.1,
        'r_process_solar': 0.95,
        'UFE_fraction': 0.05,
        'neutrino_unification': 0.995,
        'A_predicted': 254,
        'pp_dominant': '<0.1 PeV SED',
        'flux_match': '~IceCube background for LLAGNs',
        'outflow_neutrino': 0.70,
        'inflow_neutrino': 0.30,
        'rho_vac_ratio': '~10^{-38}',
        'flux_predict': '~2%',
    },
    
    # ───────────────────────────────────────────────────────────────────────────
    # RELATED MODELS
    # ───────────────────────────────────────────────────────────────────────────
    
    'related_models': [
        'UQFFAstrophysicalSystemsModel',
        'compute_fokker_planck_spectrum',
        'CRPTermModel',
    ],
}


def get_tsp_qscope_superconductive_params() -> Dict[str, Any]:
    """
    Get parameters formatted for TSP Q-Scope Superconductive Framework.
    
    Returns dictionary compatible with:
        - GinzburgLandauModel
        - BogoliubovDeGennesModel
        - DPMAttractionModel
        - QWaveResonanceModel
    """
    p = TSP_QSCOPE_SUPERCONDUCTIVE_PARAMS
    return {
        # Q-scope measurements
        'A1': p['A1_initial'],
        'A2': p['A2_channel'],
        'f_primary': p['f_primary_final'],
        'dT': p['dT_final_ms'] / 1000,      # Convert to seconds
        'f_dT': p['f_dT_final'],
        
        # Superconducting parameters
        'Phi_0': p['Phi_0'],
        'coherence_length': p['coherence_length_m'],
        'Delta_gap': p['Delta_gap_final'],
        
        # DPM parameters
        'k_DPM': p['k_DPM'],
        'f_dp': p['f_dp'],
        'U_dp': p['U_dp_computed'],
        
        # Brain wave mapping
        'brainwave_bands': p['brainwave_bands'],
        
        # Resolution status
        'verification': p['verification'],
    }


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ENTRY (for testing)
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("CondensedPhysics_InputData.py")
    print("=" * 70)
    print("\nAvailable Event Parameters:")
    print(f"  - GW170817: {len(GW170817_PARAMS)} parameters")
    print(f"  - JCAP_COSMOLOGY: {len(JCAP_COSMOLOGY_PARAMS)} parameters")
    print(f"  - RACS_J0320_35: {len(RACS_J0320_35_PARAMS)} parameters")
    print(f"  - YANG_MILLS: {len(YANG_MILLS_PARAMS)} parameters")
    print(f"  - RIEMANN_HYPOTHESIS: {len(RIEMANN_HYPOTHESIS_PARAMS)} parameters")
    print(f"  - P_VS_NP: {len(P_VS_NP_PARAMS)} parameters")
    print(f"  - DUST_YIELD: {len(DUST_YIELD_PARAMS)} parameters")
    print(f"  - SGRA_STAR: {len(SGRA_STAR_PARAMS)} parameters")
    print(f"  - SHEAR_CHI_SQUARED: {len(SHEAR_CHI_SQUARED_PARAMS)} parameters")
    print(f"  - TIME_ASYMMETRY: {len(TIME_ASYMMETRY_PARAMS)} parameters")
    print(f"  - NUCLEAR_BINDING_SHELL: {len(NUCLEAR_BINDING_SHELL_PARAMS)} parameters")
    print(f"  - SOLAR_WIND_PARKER_PROBE: {len(SOLAR_WIND_PARKER_PROBE_PARAMS)} parameters")
    print(f"  - ALPHA_BEC_LENR: {len(ALPHA_BEC_LENR_PARAMS)} parameters")
    print(f"  - COSMIC_EGG_HYPERGRAPH: {len(COSMIC_EGG_HYPERGRAPH_PARAMS)} parameters")
    print(f"  - UQFF_COMPRESSED_SUMMARY: {len(UQFF_COMPRESSED_SUMMARY_PARAMS)} parameters")
    print(f"  - UQFF_MASTER_FRAMEWORK: {len(UQFF_MASTER_FRAMEWORK_PARAMS)} parameters")
    print(f"  - DATASET_VERIFICATION_2025: {len(DATASET_VERIFICATION_2025_PARAMS)} parameters")
    print(f"  - RPROCESS_DEFAULTS: {len(RPROCESS_DEFAULTS)} parameters")
    print("\nExample - GW170817 ejecta:")
    print(f"  M_ej_total = {GW170817_PARAMS['M_ej_total']} M_☉")
    print(f"  v_ej = {GW170817_PARAMS['v_ej']} c")
    print(f"  Ye_effective = {GW170817_PARAMS['Ye_effective']}")
    print(f"  r_process_solar = {GW170817_PARAMS['r_process_solar'] * 100}%")
    print("\nExample - JCAP Cosmology:")
    print(f"  ρ_Λ = {JCAP_COSMOLOGY_PARAMS['rho_Lambda_J']:.2e} J/m³")
    print(f"  ρ_DM_local = {JCAP_COSMOLOGY_PARAMS['rho_DM_local']:.2e} J/m³")
    print(f"  λ_vac_cosmic = {JCAP_COSMOLOGY_PARAMS['lambda_vac_cosmic']:.2e} J/m³")
    print("\nExample - RACS J0320-35 Quasar Jets:")
    print(f"  v_jet = {RACS_J0320_35_PARAMS['v_jet_c']} c")
    print(f"  Eddington ratio = {RACS_J0320_35_PARAMS['eddington_ratio']}×")
    print(f"  Asymmetry ratio = {RACS_J0320_35_PARAMS['jet_asymmetry_ratio']}")
    print("\nExample - Yang-Mills Mass Gap:")
    print(f"  Mass gap E_1 = {YANG_MILLS_PARAMS['E1_MeV']} MeV")
    print(f"  Gauge group = {YANG_MILLS_PARAMS['gauge_group']}")
    print(f"  Higgs VEV = {YANG_MILLS_PARAMS['v_H_GeV']} GeV")
