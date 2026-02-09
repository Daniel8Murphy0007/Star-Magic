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
    print(f"  - UQFF_MASTER_FRAMEWORK: {len(UQFF_MASTER_FRAMEWORK_PARAMS)} parameters")
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
