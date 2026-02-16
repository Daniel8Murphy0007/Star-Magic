# =============================================================================
# CONSTANTS - Extracted from source*.js (163 files)
# =============================================================================

EXTRACTED_CONSTANTS = {
    # --- Fundamental ---
    'G': 6.67430e-11,
    'PI': 3.141592653589793,
    'c': 299792458.0,
    'c_light': 3e8,  # Speed of light
    'hbar': 1.054571817e-34,
    'proton_mass': 1.673e-27,  # Proton mass
    'q_charge': 1.602e-19,  # Proton charge

    # --- Magnetic ---
    'B': 1e-6,  # Static magnetic field (T)
    'B0': 2e10,  # Initial B field
    'B0_G': 1e4,  # Initial magnetic field (Gauss)
    'B_crit': 1e11,  # Critical B field
    'B_ref': 0.4,  # T (reference max)
    'B_s_max': 0.4,  # T (sunspot)
    'B_s_min': 1e-4,  # T (quiet Sun)
    'B_t': 1e-9,  # Current magnetic field

    # --- Vacuum ---
    'C_vac': 1.0,
    'Delta_E_vac': 1e-8,  # Vacuum energy differential
    'Delta_Evac': 6.381e-36,
    'E_vac_ISM': 1e-10,  # Vacuum energy in ISM
    'E_vac_neb': 1e-9,  # Vacuum energy in nebula
    'Evac_ISM': 7.09e-37,
    'Evac_neb': 7.09e-36,
    'delta_rho_over_rho': 1e-5,  # Density perturbation
    'rhoUa': 7.09e-36,
    'rho_A': 1e-23,
    'rho_ICM': 1e-26,  # Intracluster medium density (kg/m^3)
    'rho_dust': 1e-23,  # Dust density (kg/m^3)
    'rho_fluid': 1e17,  # Fluid density
    'rho_gas': 1e-24,  # Gas density (kg/m^3)
    'rho_sw': 1e-21,
    'rho_v': 6e-27,
    'rho_vac_A': 1e-30,
    'rho_vac_SCm': 7.09e-37,  # SCm vacuum density
    'rho_vac_UA': 7.09e-36,  # UA vacuum density
    'rho_wind': 1e-21,  # Wind density (kg/m^3)

    # --- Time ---
    'M_dot_0': 0.01,  # Initial mass accretion rate factor
    'M_dot_factor': 1e5,
    'T_ICM': 3e7,  # ICM temperature (K)
    'T_merger': 3e9 * 365.25,
    'T_val': 1e7,
    'dt_ns': 0.1,
    't_Hubble': 13.8e9 * 3.15576e7,  # Hubble time in seconds
    't_Hubble_gyr': 13.8,  # Hubble time in Gyr
    't_example': 5000 * 3.15576e7,
    't_outburst': 1e7 * 365.25,
    't_scale': 1e16,  # Constants
    'tau_B': 4000 * 3.15576e7,  # B decay timescale (not used)
    'tau_Omega': 10000 * 3.15576e7,  # Omega decay timescale
    'tau_SF': 5e6 * 3.156e7,  # Star formation timescale (5 Myr in seconds)
    'tau_SN': 1 * 3.156e7,  # SN decay timescale (s)
    'tau_acc': 9e9 * 3.15576e7,  # Accretion timescale (9 billion years)
    'tau_decay': 3.5 * 365.25,
    'tau_erosion': 1e6 * 3.156e7,  # Erosion timescale
    'tau_exp': 1e6 * 3.156e7,  # Expansion timescale (s)

    # --- Mass ---
    'M': 1.4 * 1.989e30,  # 1.4 M_sun in kg
    'M0': 400000.0 * 1.989e30,  # Initial mass (kg) - 400,000 Msun
    'M_BH': 4e6 * 1.989e30,  # Sgr A* mass
    'M_DM_default': 1e40,  # Default DM mass
    'M_DM_factor': 0.1,  # Dark matter mass fraction
    'M_M31': 1.5e12 * 1.989e30,  # Andromeda mass (M31 host)
    'M_M87': 6.5e12 * 1.989e30,  # M87 mass (central giant)
    'M_NGC1275': 1e12 * 1.989e30,  # NGC 1275 (Perseus A) mass
    'M_SN0': 1.4,
    'M_bh': 8.15e36,
    'M_bulge': 5e11 * 1.989e30,  # Bulge mass
    'M_bullet': 1e14 * 1.989e30,  # Bullet subcluster mass
    'M_cluster': 1.2e15 * 1.989e30,  # Total cluster mass (1.2 quadrillion M☉)
    'M_companion': 8e10 * 1.989e30,  # NGC 5195 mass (80 billion M☉)
    'M_disk': 3e11 * 1.989e30,  # Disk mass
    'M_halo': 1e12 * 1.989e30,
    'M_initial': 4.3e6 * 1.989e30,  # 4.3 million solar masses in kg
    'M_initial_sun': 240.0,
    'M_mag': 1e40,
    'M_main': 2e15 * 1.989e30,  # Main cluster mass
    'M_shell': 2e10 * 1.989e30,  # Shell mass from merger
    'M_star': 20 * 1.989e30,
    'M_sun': 1.989e30,
    'M_sun_val': 1.989e30,
    'Mbh': 8.15e36,

    # --- Distance ---
    'd_M31': 2e23,  # Distance to M31 (200 kpc)
    'd_g': 8.178e3,
    'd_sep': 1e22,  # Separation distance (10 kpc)
    'density': 100,  # Particles/cm³
    'dg': 2.55e20,
    'distance': 1.5e5 * 3.086e16,  # Light years to meters
    'dpm_curv': 1e-22,
    'dpm_phase': 2.36e-3,
    'r': 1e4,  # 10 km radius
    'r_BH': 2.83e16,  # Distance to Sgr A*
    'r_HII': 5e21,  # HII region scale (5 kpc)
    'r_core': 1e23,  # Core radius (100 kpc)
    'r_j': 1e10,
    'r_s': 20000,
    'r_shell': 2e22,  # Shell radius (20 kpc)
    'r_star': 1e10,
    'radius': 200 * 6.96e8,  # Solar radii
    'relative_velocity': 3e5,  # m/s

    # --- Cosmological ---
    'H0': 2.184e-18,  # Different Hubble constant
    'Hz': 2.269e-18,  # H(z) in s^-1
    'Lambda': 1.1e-52,  # Cosmological constant
    'Omega_g': 7.3e-16,
    'fTHz': 1e12,

    # --- Quantum ---
    'Delta_x_Delta_p': 1e-68,  # Uncertainty principle
    'delta_x': 1e-10,  # Position uncertainty
    'epsilon': 0.01,  # Constants
    'epsilon_sw': 0.001,
    'integral_psi': 1.0,  # Wavefunction integral

    # --- Other ---
    'A': 1,  # Mass number
    'A_osc': 1e10,  # Oscillatory amplitude
    'C_concentration': 1.0,
    'DPM_gravity': 1.0,
    'DPM_momentum': 0.93,
    'DPM_stability': 0.01,
    'Delta': 1.3e6 * 1.602e-19,
    'E_0': 0.1,  # Initial erosion factor
    'E_cm': 2.18e-6,
    'E_cm_astro': 1.24e24,  # Buoyancy parameters
    'E_react': 1e46,  # J (reaction energy)
    'F0': 1.83e71,
    'F_CNB': 9.07e-42,  # Cosmic Neutrino Background force
    'F_super': 1e12,  # Superconductor frequency factor
    'Fsuper': 6.287e-19,
    'HSCm': 1.0,
    'H_SCM': 1.0,
    'L0_W': 5e28,  # Initial luminosity (W)
    'L_Halpha': 1e40,  # H-alpha luminosity (W)
    'L_UV': 1e43,  # UV luminosity (W)
    'L_X': 1e41,  # X-ray luminosity (W)
    'L_gamma': 1e39,  # Gamma-ray luminosity (W)
    'L_jet': 1e24,  # Jet length (100 kpc)
    'L_radio': 1e39,  # Radio luminosity (W)
    'L_sun_val': 3.826e26,
    'N': 32,
    'N_galaxies': 1300,  # Number of galaxies
    'P0': 4e-8,  # Initial pressure (Pa)
    'P_AGN': 1e43,  # AGN outburst power (W)
    'P_core': 1.0,  # Core pressure factor
    'P_init': 3.76,  # Pulse period in seconds
    'P_jet': 1e40,  # Jet power (W)
    'QA': 1e-10,
    'Q_A': 1.602e-19,
    'Q_s': 1.602e-19,
    'SFR': 3 * 1.989e30,
    'S_wind': 1.0,
    'UA_SCM': 10,
    'UA_SC_m': 1e-20,  # Superconductive correction
    'UUA': 1.0,
    'U_UA': 1.0,
    'Ug1_proxy': 1.0,
    'V': 1e5,
    'V_infall': 500e3,  # Infall velocity (m/s)
    'V_infl_UA': 1e-6,
    'V_jet': 0.99 * 299792458,  # Relativistic jet velocity
    'V_merger': 4500e3,  # Merger velocity (4500 km/s)
    'V_rot': 200e3,  # Rotation velocity (m/s)
    'Z': 1,  # Atomic number
    'a_universal': 1e12,  # Superconductive parameters
    'age': 5e6,  # years
    'alpha': 0.001,
    'alpha_t': 1e-9,
    'beta_i': 0.6,
    'c_res': 3e8,
    'delta_def': 0.01,
    'delta_sw': 0.01,
    'ejecta_mass': 15 * 1.989e30,  # Solar masses
    'eta': 1e-22,
    'eta_A': 0.01,
    'evaporation_term': 1e-12,
    'evolution_timescale': 8e14,
    'explosion_energy': 1e44,  # Joules
    'fAether': 1.576e-35,
    'fDPM': 1e12,
    'fTRZ': 0.1,
    'f_Aether': 1e13,  # Aether frequency
    'f_DM': 0.85,  # Dark matter fraction
    'f_TRZ': 0.1,  # Time-reversal factor (unique to SGR 0501+4516)
    'f_dp': 0.1,  # Deep pairing factor
    'f_feedback': 1.0,
    'f_quantum': 1e14,  # Quantum frequency
    'f_res': 1.0,  # Resonance factor
    'f_sc': 1,
    'force_jet': 10.0,
    'fosc': 4.57e14,
    'fov': 45.0,
    'fquantum': 1.445e-17,
    'freact': 1e10,
    'g_earth': 9.80665,
    'gamma': 0.00005,
    'gas_mass': 1e4 * 1.989e30,  # Solar masses
    'gas_mass_sun': 10000.0,
    'gas_v': 1e5,  # Gas velocity for EM (m/s)
    'gridPoints': 10,
    'h_disk': 1e21,  # Disk scale height (1 kpc)
    'ionizing_flux': 1e8,  # Photons/m²/s
    'ionizing_stars': 4,  # Theta1 Ori C and companions
    'k1': 1.5,
    'k2': 1.2,
    'k3': 1.8,
    'k4': 2.0,
    'k4_res': 1.0,
    'k_3': 1.8,  # Coupling constant
    'k_4': 1e-40,  # Aether constant
    'k_AGN': 1e-13,  # AGN coupling constant
    'k_B': 1.380649e-23,
    'k_DE': 1e-16,
    'k_LENR': 1e-10,
    'k_LG': 1e-15,  # Local Group coupling
    'k_act': 1e-14,
    'k_asym': 1e-14,  # Asymmetry coupling
    'k_cluster': 1e-16,  # Cluster coupling constant
    'k_dust': 1e-15,  # Dust lane coupling
    'k_jet': 1e-12,  # Jet coupling constant
    'k_merger': 1e-14,  # Merger coupling constant
    'k_neutron': 1e10,
    'k_osc': 1.0,
    'k_rel': 1e-10,
    'k_tidal': 1e-13,  # Tidal coupling constant
    'kappa': 0.0005,
    'kappa_t': 5e-4,
    'kpc_val': 3.086e19,
    'lambda_i': 1.0,
    'learningRate': 0.001,
    'learning_rate': 0.01,
    'length': 2e17,  # meters
    'ly_to_m': 9.461e15,
    'ly_val': 9.461e15,
    'm': 2.0,
    'm_e': 9.10938356e-31,
    'm_factor': 1.0,
    'm_n': 1.674927498e-27,
    'm_p': 1.672621898e-27,  # Initialize variable storage
    'magnetic_field': 1e-4,  # T
    'mass': 9 * 1.989e30,  # Solar masses
    'metallicity': 0.2,  # Solar units
    'mouseSensitivity': 0.1,
    'movementSpeed': 2.5,
    'mu0': 4,
    'n': 1,
    'n_HII': 3000,  # Number of HII regions
    'n_e': 1e6,
    'n_t': 1.0,
    'ns_mass': 1.4 * 1.989e30,
    'num_magnetic_strings': 10,
    'num_strings': 1e9,
    'omega0': 1e-16,  # Angular frequency (rad/s)
    'omega_dot': 1e-3,
    'omega_i': 1e-8,
    'omega_osc': 2,
    'omega_s': 2.5e-6,  # rad/s (solar cycle)
    'p_core': 1.0,
    'pc_val': 3.086e16,
    'period': 41.4 * 24,
    'phi_hat': 1.0,
    'pi_val': 3.141592653589793,
    'pillar_count': 3,  # Famous pillars
    'pitch': 89.0,
    'precession_angle_deg': 30.0,  # Precession angle (degrees)
    'primary_mass': 5e10 * 1.989e30,  # Solar masses
    'progenitor_mass': 20 * 1.989e30,  # Solar masses
    'pulsation_amplitude': 0.1,  # Magnitude
    'q': 1.602176634e-19,  # Initialize variable storage
    'scale_EM': 1e-12,  # EM scaling factor
    'secondary_mass': 1e10 * 1.989e30,  # Solar masses
    'separation': 5e20,  # meters
    'sigma_n': 1e-4,
    'sigma_v': 700e3,  # Velocity dispersion (m/s)
    'size': 24 * 3.086e16,  # Light years across
    'spin_factor': 0.3,  # Spin factor
    'stars_count': 10000,
    't': 1000.0,
    'temperature': 50,  # K
    'tidal_radius': 1e21,  # meters
    'tn': 1e-10,
    'trap_mass': 1e3 * 1.989e30,  # Solar masses
    'trigger_term': 1e-10,
    'uv_luminosity': 1e6,  # Solar luminosities
    'v': 1.7e6,
    'v_surf': 1e6,  # Surface velocity
    'v_sw': 4e5,
    'v_wind': 2e6,  # Wind velocity (m/s)
    'visc': 0.0001,
    'x2': -1.35e172,  # UQFF coupling constant
    'yaw': -90.0,
    'z_gal': 0.016,  # Galaxy redshift

}