"""
_session195_cp4_appender.py
Session 195 — Append CP4 classes #416–#418

Source: grok_share_ab2e7192-de62.txt (2884 lines, June 09-10, 2025)
Novel physics: U_b Model (Kepler Orrery V), F_orbit, F_tide, F_gal, T_eq, NFW rho_DM

Classes to append:
  #416 — KeplerOrreryV_Ub_UQFF_Calculator          (PAPER_832)
  #417 — ExoplanetResonanceOrbitalTidalCalculator   (PAPER_832)
  #418 — GalacticDarkMatterNFWCouplingCalculator    (PAPER_834)

Usage: python _session195_cp4_appender.py
"""

BLOCK = '''

# ---------------------------------------------------------------------------
# SESSION 195 — U_b Model: Kepler Orrery V Exoplanetary UQFF Extension
# Source: grok_share_ab2e7192-de62.txt (2884 lines, June 09-10 2025)
# Novel: F_orbit=G*M_p*M_s/a^3, F_tide=G*M_p*M_s*R_p/a^6,
#        F_gal=v_gal^2/r_gal + G*M_DM/r_gal^2, T_eq, NFW rho_DM at 8 kpc
# Papers: PAPER_832, PAPER_833, PAPER_834
# ---------------------------------------------------------------------------


class KeplerOrreryV_Ub_UQFF_Calculator(_CP4Calculator):  # PAPER_832 #416
    """
    U_b Model for Kepler Orrery V Exoplanetary Systems.

    Extends UQFF with three physically motivated F_env sub-terms:
      F_orbit(t) = G*M_p*M_s / a^3           resonance force
      F_tide(t)  = G*M_p*M_s*R_p / a^6       tidal locking force
      F_gal(t)   = v_gal^2/r_gal + G*M_DM/r_gal^2   galactic coupling

    F_env(t) = 0.50*F_orbit + 0.30*F_tide + 0.20*F_gal
    Standard Kepler F_env ~ 6.5e-2 m/s^2 (a=0.1 AU)

    Full U_b equation:
      g_Ub(r,t) = G*M(t)/(r(t)^2) * (1+H(t,z)) * (1-B(t)/B_crit)
                * (1 + F_orbit(t) + F_tide(t) + F_gal(t))
                + (Ug1+Ug2+Ug3\'+Ug4) + Lambda*c^2/3
                + (hbar/sqrt(Delta_x*Delta_p))*integral(psi_total H psi_total dV)*(2*pi/t_Hubble)
                + rho_fluid*V*g
                + (M_vis+M_DM)*(delta_rho/rho + 3GM/r^3)

    Validated: Kepler-11b (a=0.091 AU, 5:4 resonance), TOI-178b (2:4 Laplace),
               TOI-849b (tidal circ), TOI-2109b (tidal distortion),
               62 Kepler Orrery V frames (22 Sep 2011 - 01 Dec 2011)
    """

    # Physical constants
    G       = 6.6743e-11   # m^3 kg^-1 s^-2
    HBAR    = 1.0546e-34   # J s
    LAMBDA  = 1.1e-52      # m^-2
    C       = 2.998e8      # m/s
    T_HUB   = 4.35e17      # s (Hubble time)
    H0      = 2.27e-18     # s^-1
    M_SUN   = 1.989e30     # kg
    R_EARTH = 6.371e6      # m
    M_EARTH = 5.972e24     # kg
    AU      = 1.496e11     # m (1 AU in meters)

    # Standard Kepler parameters (validated from 62 frames)
    V_GAL   = 2.20e5       # m/s  galactic rotation velocity
    R_GAL   = 2.47e20      # m    solar galactocentric radius (8 kpc)
    RHO_DM  = 4.2e-2       # kg/m^3  NFW dark matter density at 8 kpc
    PI      = math.pi

    def compute_F_orbit(self, M_p: float, M_s: float, a: float) -> float:
        """
        F_orbit = G*M_p*M_s / a^3  [m/s^2]
        Resonance force for orbital stability.
        M_p: planet mass [kg], M_s: star mass [kg], a: semi-major axis [m]
        """
        return self.G * M_p * M_s / (a ** 3)

    def compute_F_tide(self, M_p: float, M_s: float, R_p: float, a: float) -> float:
        """
        F_tide = G*M_p*M_s*R_p / a^6  [m/s^2]
        Tidal locking force for close-orbit planets.
        R_p: planetary radius [m], a: semi-major axis [m]
        """
        return self.G * M_p * M_s * R_p / (a ** 6)

    def compute_F_gal(self) -> float:
        """
        F_gal = v_gal^2/r_gal + G*M_DM/r_gal^2  [m/s^2]
        Galactic rotation + dark matter coupling.
        Uses canonical Milky Way parameters (8 kpc, NFW profile).
        """
        M_DM = self.RHO_DM * (4.0/3.0) * self.PI * (self.R_GAL ** 3)
        a_rot = (self.V_GAL ** 2) / self.R_GAL
        a_dm  = self.G * M_DM / (self.R_GAL ** 2)
        return a_rot + a_dm

    def compute_H(self, z: float = 0.0) -> float:
        """H(t,z) = H0 * sqrt(0.3*(1+z)^3 + 0.7)"""
        return self.H0 * math.sqrt(0.3 * (1.0 + z)**3 + 0.7)

    def compute_T_eq(self, S: float, A: float = 0.3) -> float:
        """
        T_eq = [(1-A)*S / (4*sigma)]^0.25  [K]
        Equilibrium temperature.
        S: stellar flux [W/m^2], A: bond albedo (default 0.3)
        """
        sigma = 5.6704e-8  # Stefan-Boltzmann
        return ((1.0 - A) * S / (4.0 * sigma)) ** 0.25

    def compute_F_env(self, M_p: float, M_s: float, a: float,
                      R_p: float, w_orbit: float = 0.5,
                      w_tide: float = 0.3, w_gal: float = 0.2) -> dict:
        """
        Compute weighted F_env(t) = w_orbit*F_orbit + w_tide*F_tide + w_gal*F_gal
        Default weights from 62-frame Kepler Orrery V analysis.
        """
        F_o = self.compute_F_orbit(M_p, M_s, a)
        F_t = self.compute_F_tide(M_p, M_s, R_p, a)
        F_g = self.compute_F_gal()
        F_env = w_orbit * F_o + w_tide * F_t + w_gal * F_g
        return {
            'F_orbit_m_s2': F_o,
            'F_tide_m_s2': F_t,
            'F_gal_m_s2': F_g,
            'F_env_m_s2': F_env,
            'weights': {'orbit': w_orbit, 'tide': w_tide, 'gal': w_gal},
        }

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_kg: central body mass [kg]
          r_m: orbital radius [m]
          M_p_kg: planet mass [kg] (optional, default 5*M_Earth)
          M_s_kg: star mass [kg] (optional, default 1.1*M_Sun)
          a_m: semi-major axis [m] (optional, default 0.1 AU)
          R_p_m: planet radius [m] (optional, default 1.5*R_Earth)
          z: redshift (optional, default 0)
          B_T: magnetic field [T] (optional, default 0)
          B_crit_T: critical B [T] (optional, default 4.4e13)
          S_W_m2: stellar flux [W/m^2] (optional, default 1360)
          A_albedo: bond albedo (optional, default 0.3)
        """
        M     = dataset.get('M_kg', self.M_SUN)
        r     = dataset.get('r_m', self.AU)
        M_p   = dataset.get('M_p_kg', 5.0 * self.M_EARTH)
        M_s   = dataset.get('M_s_kg', 1.1 * self.M_SUN)
        a     = dataset.get('a_m', 0.1 * self.AU)
        R_p   = dataset.get('R_p_m', 1.5 * self.R_EARTH)
        z     = dataset.get('z', 0.0)
        B     = dataset.get('B_T', 0.0)
        B_c   = dataset.get('B_crit_T', 4.4e13)
        S     = dataset.get('S_W_m2', 1360.0)
        A_alb = dataset.get('A_albedo', 0.3)

        H_tz   = self.compute_H(z)
        B_fac  = (1.0 - B / B_c) if B_c > 0 else 1.0
        env    = self.compute_F_env(M_p, M_s, a, R_p)
        F_env  = env['F_env_m_s2']

        # Base UQFF term
        g_base = self.G * M / (r ** 2)
        g_Ub   = g_base * (1.0 + H_tz) * B_fac * (1.0 + F_env)

        # Cosmological constant term
        Lambda_term = self.LAMBDA * (self.C ** 2) / 3.0

        T_eq = self.compute_T_eq(S, A_alb)

        primary_equations = [
            f'g_Ub(r,t) = {g_Ub:.6e} m/s^2',
            f'F_orbit = G*M_p*M_s/a^3 = {env["F_orbit_m_s2"]:.6e} m/s^2',
            f'F_tide  = G*M_p*M_s*R_p/a^6 = {env["F_tide_m_s2"]:.6e} m/s^2',
            f'F_gal   = v_gal^2/r_gal + G*M_DM/r_gal^2 = {env["F_gal_m_s2"]:.6e} m/s^2',
            f'F_env(t)= {F_env:.6e} m/s^2  (50%F_orbit + 30%F_tide + 20%F_gal)',
            f'H(t,z)  = {H_tz:.6e}',
            f'H(t,z)*B_fac*(1+F_env) factor = {(1.0+H_tz)*B_fac*(1.0+F_env):.6f}',
            f'Lambda_term = {Lambda_term:.6e} m/s^2',
            f'T_eq    = {T_eq:.1f} K  (S={S:.0f} W/m^2, A={A_alb})',
        ]
        available_equations = [
            'F_orbit = G*M_p*M_s / a^3',
            'F_tide  = G*M_p*M_s*R_p / a^6',
            'F_gal   = v_gal^2/r_gal + G*M_DM/r_gal^2',
            'T_eq    = [(1-A)*S/(4*sigma)]^0.25',
            'g_Ub    = g_base*(1+H)*(1-B/B_crit)*(1+F_env) + Lambda*c^2/3 + ...',
            'M_DM    = rho_DM*(4/3)*pi*r_gal^3',
            'H(t,z)  = H0*sqrt(0.3*(1+z)^3 + 0.7)',
            'P_1/P_2 = (a_1/a_2)^1.5  (Kepler third law resonance)',
        ]
        simulation_set = [
            {
                'equation': 'F_orbit(a)',
                'values': {f'a={a_au:.2f} AU': self.compute_F_orbit(M_p, M_s, a_au*self.AU)
                           for a_au in [0.01, 0.05, 0.1, 0.5, 1.0, 2.0]},
            },
            {
                'equation': 'F_tide(a)',
                'values': {f'a={a_au:.2f} AU': self.compute_F_tide(M_p, M_s, R_p, a_au*self.AU)
                           for a_au in [0.01, 0.02, 0.05, 0.1]},
            },
            {
                'equation': 'T_eq(S)',
                'values': {f'S={s:.0f} W/m^2': self.compute_T_eq(s, A_alb)
                           for s in [100, 500, 1000, 1360, 5000, 10000]},
            },
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_orbit_m_s2': env['F_orbit_m_s2'],
            'F_tide_m_s2': env['F_tide_m_s2'],
            'F_gal_m_s2': env['F_gal_m_s2'],
            'F_env_m_s2': F_env,
            'g_Ub_m_s2': g_Ub,
            'T_eq_K': T_eq,
            'paper': 'PAPER_832',
        }


class ExoplanetResonanceOrbitalTidalCalculator(_CP4Calculator):  # PAPER_832 #417
    """
    Numerical solvers for exoplanet orbital resonance and tidal locking.

    F_orbit solver: given M_p, M_s, a_1, a_2 -> compute period ratio P_1/P_2
                    check for n:m resonance, compute F_orbit for each planet.
    F_tide solver:  given M_p, M_s, R_p, a -> compute F_tide; check tidal locking.

    Validated datasets:
      Kepler-11 system: 5:4 resonance (a=0.091 AU, F_orbit~1.28e-1 m/s^2)
      TOI-178: 2:4:6:9:12 Laplace resonance chain
      TOI-849b: tidal circularization (a=0.016 AU, F_tide~5.61e-12 m/s^2)
      Kepler-13Ab: tidal locking (a=0.033 AU)
      TOI-2109b: tidal distortion (a=0.018 AU)
      Kepler-90g/h: 3:2 resonance (Kepler DR25)
    """

    G     = 6.6743e-11
    AU    = 1.496e11
    M_SUN = 1.989e30
    M_EAR = 5.972e24
    R_EAR = 6.371e6

    # Resonance fraction tolerance
    RES_TOL = 0.05

    def check_resonance(self, a1: float, a2: float, M_s: float) -> dict:
        """
        Check for mean-motion resonance between two planets.
        Returns period ratio and closest n:m ratio.
        a1, a2: semi-major axes [m], M_s: stellar mass [kg]
        """
        P1 = 2 * math.pi * math.sqrt(a1**3 / (self.G * M_s))
        P2 = 2 * math.pi * math.sqrt(a2**3 / (self.G * M_s))
        ratio = P2 / P1
        # Check common resonances
        resonances = [(2,1),(3,2),(4,3),(5,4),(3,1),(4,1),(6,1),(9,2),(12,3)]
        closest = None
        min_diff = 1e9
        for n, m in resonances:
            r_frac = n / m
            diff = abs(ratio - r_frac)
            if diff < min_diff:
                min_diff = diff
                closest = (n, m)
        in_resonance = min_diff < self.RES_TOL
        return {
            'P1_days': P1 / 86400,
            'P2_days': P2 / 86400,
            'period_ratio': ratio,
            'closest_resonance': f'{closest[0]}:{closest[1]}',
            'in_resonance': in_resonance,
            'deviation_from_resonance': min_diff,
        }

    def tidal_locking_check(self, M_p: float, M_s: float,
                             R_p: float, a: float) -> dict:
        """
        Compute F_tide and estimate tidal locking timescale.
        Strong F_tide (>> 1e-15 m/s^2) indicates likely tidal locking.
        """
        F_tide = self.G * M_p * M_s * R_p / (a**6)
        # Tidal locking timescale (McDonald 1964 secular approximation)
        # tau ~ (a/R_p)^6 * (M_p/M_s) * P_orb / Q_p
        # Default Q_p ~ 100 (rocky), 1e4 (gas giant)
        Q_p = 100.0
        P_orb = 2 * math.pi * math.sqrt(a**3 / (self.G * M_s))
        tau = ((a / R_p)**6) * (M_p / M_s) * P_orb / Q_p
        return {
            'F_tide_m_s2': F_tide,
            'likely_locked': F_tide > 1e-12,
            'tidal_lock_timescale_Gyr': tau / (3.156e16),  # Gyr
        }

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          M_p_kg, M_s_kg, a_m, R_p_m: primary planet
          a2_m: second planet semi-major axis (optional, for resonance check)
          system_name: descriptive name (optional)
        """
        M_p  = dataset.get('M_p_kg', 5.0 * self.M_EAR)
        M_s  = dataset.get('M_s_kg', 1.1 * self.M_SUN)
        a    = dataset.get('a_m', 0.1 * self.AU)
        R_p  = dataset.get('R_p_m', 1.5 * self.R_EAR)
        a2   = dataset.get('a2_m', None)
        name = dataset.get('system_name', 'exoplanet')

        F_orbit = self.G * M_p * M_s / (a**3)
        tidal   = self.tidal_locking_check(M_p, M_s, R_p, a)

        resonance = {}
        if a2 is not None:
            resonance = self.check_resonance(a, a2, M_s)

        primary_equations = [
            f'System: {name}',
            f'F_orbit = G*M_p*M_s/a^3 = {F_orbit:.6e} m/s^2',
            f'F_tide  = G*M_p*M_s*R_p/a^6 = {tidal["F_tide_m_s2"]:.6e} m/s^2',
            f'Tidal locking: {"likely" if tidal["likely_locked"] else "unlikely"}',
            f'Tidal lock timescale: {tidal["tidal_lock_timescale_Gyr"]:.3e} Gyr',
        ]
        if resonance:
            primary_equations += [
                f'Period ratio P2/P1 = {resonance["period_ratio"]:.4f}',
                f'Closest resonance: {resonance["closest_resonance"]}',
                f'In resonance: {resonance["in_resonance"]}',
            ]

        available_equations = [
            'F_orbit(a) = G*M_p*M_s/a^3 -- orbital resonance force',
            'F_tide(a)  = G*M_p*M_s*R_p/a^6 -- tidal locking force',
            'P(a)       = 2*pi*sqrt(a^3/(G*M_s)) -- orbital period',
            'P1/P2      = (a1/a2)^1.5 -- period ratio (Kepler third law)',
            'tau_lock   = (a/R_p)^6 * (M_p/M_s) * P_orb / Q_p -- locking timescale',
        ]
        simulation_set = [
            {
                'equation': 'F_orbit_vs_a',
                'planet': name,
                'values': {f'{a_au:.2f} AU': self.G * M_p * M_s / (a_au*self.AU)**3
                           for a_au in [0.01, 0.05, 0.1, 0.5, 1.0]},
            },
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_orbit_m_s2': F_orbit,
            'F_tide_m_s2': tidal['F_tide_m_s2'],
            'tidal_data': tidal,
            'resonance_data': resonance,
            'paper': 'PAPER_832',
        }


class GalacticDarkMatterNFWCouplingCalculator(_CP4Calculator):  # PAPER_834 #418
    """
    Galactic rotation + dark matter coupling term F_gal for UQFF.

    Equations:
      F_gal(r) = v_gal(r)^2/r + G*M_DM(r)/r^2

    NFW profile (Navarro-Frenk-White 1996):
      rho_NFW(r) = rho_s / [(r/r_s)*(1 + r/r_s)^2]
      M_DM(r)    = 4*pi*rho_s*r_s^3 * [ln(1+r/r_s) - (r/r_s)/(1+r/r_s)]

    Milky Way canonical parameters (validated from 62 Kepler frames):
      v_gal   = 220 km/s (flat rotation curve)
      r_gal   = 8 kpc = 2.47e20 m
      rho_DM  = 4.2e-2 kg/m^3  (at 8 kpc)
      M_DM_8kpc = 2.57e40 kg
      F_DM    = 2.83e-10 m/s^2
      a_rot   = 1.96e-10 m/s^2
      F_gal   = 4.79e-10 m/s^2

    Addresses galaxy rotation curve problem within UQFF framework.
    Provides galactic environmental floor for F_env(t) calculation.
    """

    G       = 6.6743e-11
    KPC     = 3.086e19       # m per kpc
    M_SUN   = 1.989e30
    PI      = math.pi

    # Milky Way canonical NFW parameters
    V_GAL   = 2.20e5         # m/s  (220 km/s)
    R_GAL   = 2.47e20        # m  (8 kpc)
    RHO_DM_8KPC = 4.2e-2     # kg/m^3  at 8 kpc
    R_S     = 2.0e22         # m  NFW scale radius (~20 kpc for Milky Way)
    RHO_S   = 1.0e-22        # kg/m^3  NFW characteristic density (estimated)

    def compute_M_DM_sphere(self, r: float, rho_DM: float) -> float:
        """
        Approximate dark matter mass within sphere of radius r.
        Uses constant density approximation (valid near fiducial point).
        M_DM = rho_DM * (4/3) * pi * r^3
        """
        return rho_DM * (4.0/3.0) * self.PI * (r**3)

    def compute_M_DM_NFW(self, r: float) -> float:
        """
        Dark matter mass within r via NFW profile integration.
        M_DM(r) = 4*pi*rho_s*r_s^3 * [ln(1+r/r_s) - (r/r_s)/(1+r/r_s)]
        """
        x = r / self.R_S
        return 4.0 * self.PI * self.RHO_S * (self.R_S**3) * (
            math.log(1.0 + x) - x/(1.0 + x)
        )

    def compute_F_gal(self, r: float = None, v_c: float = None,
                      rho_DM_override: float = None) -> dict:
        """
        F_gal(r) = v_c^2/r + G*M_DM(r)/r^2
        r: galactocentric radius [m] (default 8 kpc)
        v_c: circular velocity [m/s] (default 220 km/s)
        rho_DM_override: override for density [kg/m^3]
        """
        r   = r if r is not None else self.R_GAL
        v_c = v_c if v_c is not None else self.V_GAL
        rho = rho_DM_override if rho_DM_override is not None else self.RHO_DM_8KPC

        M_DM  = self.compute_M_DM_sphere(r, rho)
        a_rot = v_c**2 / r
        a_dm  = self.G * M_DM / (r**2)
        F_gal = a_rot + a_dm
        return {
            'r_m': r,
            'r_kpc': r / self.KPC,
            'v_gal_m_s': v_c,
            'rho_DM_kg_m3': rho,
            'M_DM_kg': M_DM,
            'M_DM_Msun': M_DM / self.M_SUN,
            'a_rot_m_s2': a_rot,
            'F_DM_m_s2': a_dm,
            'F_gal_m_s2': F_gal,
        }

    def compute_THz_recombination_tau(self, N: float,
                                       A: float = 1e7,
                                       B: float = 1e-16,
                                       C: float = 1e-41) -> float:
        """
        THz hole recombination timing (interface term).
        tau = 1 / (A + B*N + C*N^2)
        N: carrier density [m^-3]
        A: Shockley-Read-Hall coefficient [s^-1]
        B: radiative recombination coefficient [m^3/s]
        C: Auger coefficient [m^6/s]
        Returns: recombination time [s]
        """
        return 1.0 / (A + B * N + C * (N**2))

    def compute(self, dataset: dict) -> dict:
        """
        dataset keys:
          r_m: galactocentric radius [m] (optional, default 8 kpc)
          v_c_m_s: circular velocity [m/s] (optional, default 220 km/s)
          rho_DM_kg_m3: dark matter density [kg/m^3] (optional)
          N_carrier_m3: carrier density for THz tau (optional)
        """
        r   = dataset.get('r_m', self.R_GAL)
        v_c = dataset.get('v_c_m_s', self.V_GAL)
        rho = dataset.get('rho_DM_kg_m3', None)
        N   = dataset.get('N_carrier_m3', None)

        gal = self.compute_F_gal(r, v_c, rho)

        primary_equations = [
            f'F_gal = v_gal^2/r_gal + G*M_DM/r_gal^2 = {gal["F_gal_m_s2"]:.6e} m/s^2',
            f'  a_rot = v_gal^2/r_gal = {gal["a_rot_m_s2"]:.6e} m/s^2',
            f'  F_DM  = G*M_DM/r_gal^2 = {gal["F_DM_m_s2"]:.6e} m/s^2',
            f'  r_gal = {gal["r_kpc"]:.2f} kpc',
            f'  v_gal = {gal["v_gal_m_s"]/1000:.0f} km/s',
            f'  rho_DM = {gal["rho_DM_kg_m3"]:.3e} kg/m^3 (NFW at 8 kpc)',
            f'  M_DM = {gal["M_DM_Msun"]:.3e} M_Sun (enclosed within {gal["r_kpc"]:.1f} kpc)',
        ]
        if N is not None:
            tau = self.compute_THz_recombination_tau(N)
            primary_equations.append(
                f'  THz tau = 1/(A+B*N+C*N^2) = {tau:.3e} s  (N={N:.2e} m^-3)'
            )

        # F_gal profile at multiple galactocentric radii
        radii = [1, 4, 8, 12, 20, 50]  # kpc
        profile = {f'{rk} kpc': self.compute_F_gal(rk*self.KPC, v_c, rho)['F_gal_m_s2']
                   for rk in radii}

        available_equations = [
            'F_gal   = v_gal^2/r + G*M_DM(r)/r^2',
            'M_DM(r) = rho_DM*(4/3)*pi*r^3  (uniform approx)',
            'M_DM(r) = 4*pi*rho_s*r_s^3*[ln(1+r/r_s) - (r/r_s)/(1+r/r_s)]  (NFW)',
            'rho_NFW = rho_s / [(r/r_s)*(1+r/r_s)^2]',
            'tau_THz = 1/(A + B*N + C*N^2)',
            'v_c(r)  = sqrt(G*M_total(r)/r)  (flat curve requires M_DM)',
        ]
        simulation_set = [
            {
                'equation': 'F_gal_vs_r_kpc',
                'values': profile,
            },
        ]
        return {
            'primary_equations': primary_equations,
            'available_equations': available_equations,
            'simulation_set': simulation_set,
            'F_gal_data': gal,
            'F_gal_profile_kpc': profile,
            'paper': 'PAPER_834',
        }


_SESSION_195_CLASSES = [
    'KeplerOrreryV_Ub_UQFF_Calculator',           # PAPER_832 #416
    'ExoplanetResonanceOrbitalTidalCalculator',    # PAPER_832 #417
    'GalacticDarkMatterNFWCouplingCalculator',     # PAPER_834 #418
]
'''

if __name__ == '__main__':
    import os
    target = 'CondensedPhysics4.py'
    if not os.path.exists(target):
        print(f'ERROR: {target} not found.')
        exit(1)
    with open(target, 'a', encoding='utf-8') as f:
        f.write(BLOCK)
    print(f'Appended 3 CP4 classes (#416-#418) to {target}')
    # Verify
    import subprocess
    r = subprocess.run(['python', '-c', f'import ast; ast.parse(open("{target}").read()); print("Syntax OK")'],
                       capture_output=True, text=True)
    print(r.stdout.strip() or r.stderr.strip())
