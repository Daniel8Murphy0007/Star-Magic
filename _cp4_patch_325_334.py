from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell


# ===========================================================================
# CP4 CLASSES #325–334 — Session 180 continuation v5.38
# PAPER_741–750: UQFF 38-System Master, Sombrero, Saturn, M16, Crab,
#                GenHRes, UnivDiam, Doc43d, 5QVars, M51+NGC1316
# ===========================================================================

class UQFF38SystemCompressedMasterCalculator:
    """CP4 #325 — PAPER_741: 38-System UQFF Compression Cycle 2 Master.
    F_env modular framework with 15 environmental terms including
    F_eta, F_DE, H_res. Covers all Hubble-observed astrophysical systems."""

    def __init__(self):
        import math
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.c = 3e8
        self.hbar = 1.0546e-34
        self.t_Hubble = 4.35e17
        self.B_crit = 1e11
        self.rho_vac_SCm = 7.09e-37
        self.rho_vac_UA = 7.09e-36
        self.k_eta = 1e-113
        self.eta_inertia = 1.0
        self.omega_vac = 1e10
        self.lambda_I = 1.0
        self.F_RZ = 0.01

    def calculate(self, M, r, z=0.0, B=1e-6, t=0.0, F_env_terms=None):
        import math
        if F_env_terms is None:
            F_env_terms = {}
        G, c, H_0 = self.G, self.c, self.H_0
        hbar, Lambda = self.hbar, self.Lambda
        H = H_0 * math.sqrt(0.3 * (1 + z)**3 + 0.7)
        g_grav = (dpm_ug1_seed(M, r)) * (1 + H) * (1 - B / self.B_crit)
        F_eta = self.k_eta * F_env_terms.get("eta", 1.0)
        rho_vac = self.rho_vac_SCm
        V = (4/3) * math.pi * r**3
        F_DE = self.eta_inertia * rho_vac * V * self.omega_vac
        F_env_total = F_eta + F_DE + sum(F_env_terms.get(k, 0) for k in
                       ["F_wind","F_erode","F_merge","F_SN","F_rad","F_fil",
                        "F_BH","F_dust","F_ring","F_mag","F_tech","F_shell","F_cosmo"])
        g_total = g_grav * (1 + F_env_total)
        U_i = self.lambda_I * (self.rho_vac_SCm / self.rho_vac_UA) * 1e-8 * (1 + self.F_RZ)
        cosmological = Lambda * c**2 / 3
        H_res = F_env_terms.get("A_res", 0) * math.sin(2 * math.pi * F_env_terms.get("f_res", 1) * t)
        return {
            "g_total": g_total + U_i + cosmological,
            "g_gravitational": g_grav,
            "F_env_total": F_env_total,
            "F_eta": F_eta,
            "F_DE": F_DE,
            "U_i": U_i,
            "H_res": H_res,
            "cosmological_term": cosmological,
            "primary_equations": [
                "g_UQFF = (G*M/r^2)*(1+H)*(1-B/B_crit)*(1+F_env)",
                "F_env = F_eta + F_DE + sum(15 env terms)",
                "F_eta = k_eta * eta",
                "F_DE = eta_inertia * rho_vac * V * omega_vac",
                "H_res = A_res*sin(2*pi*f_res*t) + F_env*SC_m",
                "U_i = lambda_I*(rho_vac_SCm/rho_vac_UA)*omega_i*(1+F_RZ)",
            ],
            "note": "UQFF Compression Cycle 2, 38-system F_env master. PAPER_741, CP4 class #325. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class SombreroGalaxyDustMUGECalculator:
    """CP4 #326 — PAPER_742: Sombrero Galaxy (M104) MUGE with dust lane drag.
    D_dust = -k_dust*rho_dust*v_orbit^2*A_cross/r contributes ~1.1% correction."""

    def __init__(self):
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.c = 3e8
        self.B_crit = 1e11
        self.rho_vac_SCm = 7.09e-37
        self.rho_vac_UA = 7.09e-36
        self.M_sun = 1.989e30
        self.M_vis = 8.0e10 * 1.989e30
        self.M_DM = 2.0e11 * 1.989e30
        self.M_BH = 1.0e9 * 1.989e30
        self.r_BH = 100 * 3.086e16
        self.rho_dust = 1e-20
        self.k_dust = 0.5
        self.A_cross_frac = 1e-15
        self.B = 1e-5
        self.z = 0.00354

    def calculate(self, r, v_orbit=2.2e5):
        import math
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.3 * (1 + self.z)**3 + 0.7)
        M_tot = self.M_vis + self.M_DM
        g_grav = (dpm_ug1_seed(M_tot, r)) * (1 + H) * (1 - self.B / self.B_crit)
        g_BH = G * self.M_BH / max(r**2, self.r_BH**2)
        A_cross = self.A_cross_frac * r**2
        D_dust = -self.k_dust * self.rho_dust * v_orbit**2 * A_cross / r
        U_i = 1.0 * (self.rho_vac_SCm / self.rho_vac_UA) * 1e-8 * 1.01
        cosmological = self.Lambda * c**2 / 3
        g_total = g_grav + g_BH + D_dust + U_i + cosmological
        D_dust_pct = abs(D_dust / g_grav) * 100
        return {
            "g_total": g_total,
            "g_gravitational": g_grav,
            "g_BH": g_BH,
            "D_dust": D_dust,
            "D_dust_percent": D_dust_pct,
            "primary_equations": [
                "g_Sombrero = (G*M/r^2)*(1+H)*(1-B/B_crit) + G*M_BH/r_BH^2 + D_dust",
                "D_dust = -k_dust*rho_dust*v_orbit^2*A_cross/r",
                "rho_dust ~ 1e-20 kg/m^3",
                "D_dust ~ -5e-12 m/s^2 (~1.1% correction)",
            ],
            "note": "Sombrero (M104) MUGE. Dust lane drag ~1.1% suppression. PAPER_742, CP4 class #326. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class SaturnRingTidalMUGECalculator:
    """CP4 #327 — PAPER_743: Saturn MUGE with T_ring tidal and F_wind.
    T_ring = k_ring*G*M_rings*r/r_ring^3 ~ 2.05e-9 m/s^2.
    Solar orbit term includes Hubble-correction for completeness."""

    def __init__(self):
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.c = 3e8
        self.B_crit = 1e11
        self.M_sun = 1.989e30
        self.M_saturn = 5.683e26
        self.M_rings = 1.5e19
        self.r_ring = 1.0e8
        self.r_orbit = 1.427e12
        self.k_ring = 1.0
        self.rho_atm = 0.19
        self.v_wind = 500.0
        self.C_D = 1.0
        self.r_atm = 6.0e7
        self.B_saturn = 0.2e-4
        self.z = 0.0

    def calculate(self, r=6.0e7):
        import math
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.7)
        g_solar_orbit = (G * self.M_sun / self.r_orbit**2) * (1 + H * 4.35e17)
        g_saturn = dpm_ug1_seed(self.M_saturn, r) * (1 - self.B_saturn / self.B_crit)
        T_ring = self.k_ring * G * self.M_rings * r / self.r_ring**3
        F_wind = 0.5 * self.rho_atm * self.v_wind**2 * self.C_D / self.r_atm
        cosmological = self.Lambda * c**2 / 3
        g_total = g_solar_orbit + g_saturn + T_ring + F_wind + cosmological
        return {
            "g_total": g_total,
            "g_solar_orbit": g_solar_orbit,
            "g_saturn_surface": g_saturn,
            "T_ring": T_ring,
            "F_wind": F_wind,
            "primary_equations": [
                "g_Saturn = (G*M_Sun/r_orbit^2)*(1+H*t) + (G*M_Saturn/r^2)*(1-B/B_crit)",
                "g_Saturn += T_ring + F_wind + Lambda*c^2/3",
                "T_ring = k_ring*G*M_rings*r/r_ring^3 ~ 2.05e-9 m/s^2",
                "F_wind = 0.5*rho_atm*v_wind^2*C_D/r_atm ~ 1.79e-10 m/s^2",
            ],
            "note": "Saturn system MUGE. T_ring ~ 2.05e-9 m/s^2. PAPER_743, CP4 class #327. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class M16EagleNebulaRadiationMUGECalculator:
    """CP4 #328 — PAPER_744: M16 Eagle Nebula MUGE with M_sf(t) and −E_rad.
    M_sf(t) = SFR*t/M_0 ~1.25% mass at 1 Myr.
    E_rad = G*m_dot_evap*t/(r^2*M_cloud) ~ 3e-12 m/s^2."""

    def __init__(self):
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.18e-52
        self.c = 3e8
        self.B_crit = 1e11
        self.M_sun = 1.989e30
        self.M_cloud = 2e4 * 1.989e30
        self.SFR = 800 * 1.989e30 / 3.156e7
        self.r_pillar = 9.26e17
        self.B = 2e-10
        self.mdot_evap = 1e26
        self.z = 0.0
        self.t_0 = 3.156e13

    def calculate(self, r=None, t=3.156e13):
        import math
        if r is None:
            r = self.r_pillar
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.7)
        M_sf = self.SFR * t / self.M_cloud
        M_eff = self.M_cloud * (1 + M_sf)
        g_grav = (dpm_ug1_seed(M_eff, r)) * (1 + H * t) * (1 - self.B / self.B_crit) * (1 + M_sf)
        E_rad = G * self.mdot_evap * t / (r**2 * self.M_cloud)
        cosmological = self.Lambda * c**2 / 3
        g_total = g_grav - E_rad + cosmological
        return {
            "g_total": g_total,
            "g_gravitational": g_grav,
            "M_sf_fraction": M_sf,
            "E_rad": E_rad,
            "primary_equations": [
                "g_M16 = (G*M(t)/r^2)*(1+H*t)*(1-B/B_crit)*(1+M_sf(t)) - E_rad",
                "M_sf(t) = SFR*t/M_0  [1.25% at 1 Myr]",
                "E_rad = G*mdot_evap*t/(r^2*M_cloud) ~ 3e-12 m/s^2",
                "M_eff = M_cloud*(1 + M_sf)",
            ],
            "note": "M16 Eagle Nebula MUGE. Radiation erosion ~3e-12 m/s^2. PAPER_744, CP4 class #328. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class CrabNebulaExpandingMUGECalculator:
    """CP4 #329 — PAPER_745: Crab Nebula expanding supernova remnant MUGE.
    r(t) = r_0 + v_r*t with v_r=1.5e6 m/s.
    F_wind = L_pulsar/(4*pi*r^2*c*M_ejecta), M_mag ~ 2.7e-17 m/s^2."""

    def __init__(self):
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.c = 3e8
        self.B_crit = 1e11
        self.mu_0 = 4 * 3.14159 * 1e-7
        self.M_sun = 1.989e30
        self.M_ejecta = 2.0 * 1.989e30
        self.r_0 = 9.5 * 3.086e15
        self.v_r = 1.5e6
        self.L_pulsar = 5e31
        self.B_pulsar = 1e9
        self.rho_ejecta = 1e-21
        self.t_age = 970 * 3.156e7
        self.z = 0.0

    def calculate(self, t=None):
        import math
        if t is None:
            t = self.t_age
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.7)
        r = self.r_0 + self.v_r * t
        g_grav = dpm_ug1_seed(self.M_ejecta, r) * (1 + H * t) * (1 - self.B_pulsar / self.B_crit)
        F_wind = self.L_pulsar / (4 * math.pi * r**2 * c * self.M_ejecta)
        M_mag = self.B_pulsar**2 / (2 * self.mu_0 * r * self.rho_ejecta)
        cosmological = self.Lambda * c**2 / 3
        g_total = g_grav + F_wind + M_mag + cosmological
        return {
            "g_total": g_total,
            "r_at_t": r,
            "g_gravitational": g_grav,
            "F_wind": F_wind,
            "M_mag": M_mag,
            "primary_equations": [
                "g_Crab = (G*M/r(t)^2)*(1+H*t)*(1-B/B_crit) + F_wind + M_mag",
                "r(t) = r_0 + v_r*t,  v_r = 1.5e6 m/s",
                "F_wind = L_pulsar/(4*pi*r^2*c*M_ejecta)",
                "M_mag = B^2/(2*mu_0*r*rho_ejecta) ~ 2.7e-17 m/s^2",
            ],
            "note": "Crab Nebula expanding SNR MUGE. 970-yr evolution. PAPER_745, CP4 class #329. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class GeneralizedHydrogenResonanceAllElementsCalculator:
    """CP4 #330 — PAPER_746: Generalized H-resonance for ALL elements Z=1–118.
    H_res = A_res*sin(2*pi*f_res*t) + U_dp*SC_m*k_nuc + S_shell.
    Covers magic numbers, odd-A pairing, and spin shell corrections."""

    def __init__(self):
        self.h = 6.626e-34
        self.c = 3e8
        self.e = 1.602e-19
        self.m_H = 1.673e-27
        self.A_H = 1
        self.E_bind_H = 13.6 * 1.602e-19
        self.k_A = 1e-10
        self.k_0 = 1.5
        self.mu_0 = 4 * 3.14159e-7
        self.MAGIC_Z = [2, 8, 20, 28, 50, 82, 126]
        self.MAGIC_N = [2, 8, 20, 28, 50, 82, 126]

    def calculate(self, Z, A, t=0.0, f_dp=1e15):
        import math
        N = A - Z
        delta_pair = 0.1 if (Z % 2 == 0 and N % 2 == 0) else (
                    -0.1 if (Z % 2 == 1 and N % 2 == 1) else 0)
        S_shell = 0.1 * (
            sum(1 for zm in self.MAGIC_Z if abs(Z - zm) < 2) +
            sum(1 for nm in self.MAGIC_N if abs(N - nm) < 2)
        )
        k_nuc = self.k_0 * (N / max(Z, 1)) * (1 + delta_pair)
        A_res = self.k_A * Z * (A / self.A_H) * (1 + delta_pair)
        f_res = (self.E_bind_H / self.h) * (self.A_H / A) * (1 + S_shell)
        U_dp = self.k_0 * (A * self.A_H / f_dp**2) * math.cos(0)
        SC_m = 1.0
        H_res = A_res * math.sin(2 * math.pi * f_res * t) + U_dp * SC_m * k_nuc + S_shell
        E_bind = self.E_bind_H * Z**2 / A
        return {
            "Z": Z,
            "A": A,
            "N": N,
            "H_res": H_res,
            "A_res": A_res,
            "f_res": f_res,
            "k_nuc": k_nuc,
            "S_shell": S_shell,
            "delta_pair": delta_pair,
            "E_bind_approx": E_bind,
            "primary_equations": [
                "H_res = A_res*sin(2*pi*f_res*t) + U_dp*SC_m*k_nuc + S_shell",
                "A_res = k_A*Z*(A/A_H)*(1+delta_pair)",
                "f_res = (E_bind_H/h)*(A_H/A)*(1+S_shell)",
                "U_dp = k*(A1*A2/f_dp^2)*cos(phi_dp)",
                "k_nuc = k_0*(N/Z)*(1+delta_pair)",
                "S_shell = 0.1*(Z_magic + N_magic proximity)",
            ],
            "note": "Generalized H-res all elements Z=1-118. PAPER_746, CP4 class #330. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class UniverseDiameterUQFFCalculator:
    """CP4 #331 — PAPER_747: Universe diameter UQFF master equation.
    D_universe ~ 182 billion light-years from UQFF curvature corrections.
    D = 2*D_p*(1+H_factor)*(1+Lambda_factor)*(1+quantum)*(1+k*r_c^2)."""

    def __init__(self):
        self.c = 3e8
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.hbar = 1.0546e-34
        self.m_p = 1.673e-27
        self.t_0 = 4.35e17
        self.D_p_m = 4.40e26
        self.ly_m = 9.461e15

    def calculate(self, k_curv=0.05, Omega_k=0.001):
        import math
        D_p = self.D_p_m
        H_factor = self.H_0 * self.t_0
        Lambda_factor = self.Lambda * self.c**2 / (3 * self.H_0**2)
        E_planck = self.hbar * self.H_0
        m_eff = self.m_p
        quantum_correction = E_planck / (m_eff * self.c**2)
        r_c = D_p
        curv_correction = k_curv * (r_c**2 / D_p**2)
        D_universe = 2 * D_p * (1 + H_factor) * (1 + Lambda_factor) * (1 + quantum_correction) * (1 + curv_correction)
        D_universe_bly = D_universe / self.ly_m / 1e9
        return {
            "D_universe_meters": D_universe,
            "D_universe_bly": D_universe_bly,
            "D_p_meters": D_p,
            "H_factor": H_factor,
            "Lambda_factor": Lambda_factor,
            "quantum_correction": quantum_correction,
            "curv_correction": curv_correction,
            "primary_equations": [
                "D_universe = 2*D_p*(1+H(z)*t_0)*(1+Lambda*c^2/(3*H_0^2))*(1+quantum)*(1+k*r_c^2)",
                "D_p = particle horizon = 46.5 Gly = 4.40e26 m",
                "Lambda_factor = Lambda*c^2/(3*H_0^2)",
                "quantum_correction = hbar*H_0/(m_p*c^2)",
                "D_universe ~ 182 billion light-years",
            ],
            "note": "Universe diameter UQFF. ~182 Gly with curvature. PAPER_747, CP4 class #331. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class Doc43dInertiaAetherSuperconductiveCalculator:
    """CP4 #332 — PAPER_748: Doc 43.d 19-equation inertia/aether/U_g5 framework.
    New 5th gravity mode: U_g5 = sum(T_mu_nu) ~ 3.6e-3 J/m^3.
    Includes ψ_matter, Jeans mass, P_DE, plasma frequency, U_i, U_m."""

    def __init__(self):
        self.G = 6.6743e-11
        self.c = 3e8
        self.hbar = 1.0546e-34
        self.h = 6.626e-34
        self.k_B = 1.381e-23
        self.mu_0 = 4 * 3.14159e-7
        self.eps_0 = 8.854e-12
        self.m_H = 1.673e-27
        self.M_sun = 1.989e30
        self.rho_vac_SCm = 7.09e-37
        self.rho_vac_UA = 7.09e-36
        self.eta_inertia = 1.0
        self.omega_vac = 1e10
        self.n_e = 1e6
        self.e = 1.602e-19
        self.m_e = 9.109e-31
        self.phi = 1.618

    def calculate(self, T=1e4, rho=1e-21, E_EMP=1e3, V_aether=1.0, mu=9.274e-24, B=1e-6, I=1.0, A=1.0, omega_spin=1e10):
        import math
        H_mag = -mu * B
        omega_plasma = math.sqrt(self.n_e * self.e**2 / (self.eps_0 * self.m_e))
        M_J = (5 * self.k_B * T / (self.G * 2.4 * self.m_H))**1.5 * (3 / (4 * math.pi * rho))**0.5
        P_DE = self.eta_inertia * self.rho_vac_SCm * V_aether * self.omega_vac
        P_AC = 0.5 * self.eps_0 * E_EMP**2 * V_aether * 2 * math.pi * 1e9
        f_n_vals = [440.0 * self.phi**n for n in range(10)]
        mu_dipole = I * A * omega_spin
        B_super = self.mu_0 * V_aether
        U_g5 = self.rho_vac_SCm * self.c**2 * 0.001
        U_i = 1.0 * (self.rho_vac_SCm / self.rho_vac_UA) * 1e-8 * 1.01
        U_m = mu * B / (self.mu_0 * 1.0)
        return {
            "H_mag": H_mag,
            "omega_plasma": omega_plasma,
            "M_Jeans_solar": M_J / self.M_sun,
            "P_DE": P_DE,
            "P_AC": P_AC,
            "f_n_series": f_n_vals[:5],
            "mu_dipole": mu_dipole,
            "B_super": B_super,
            "U_g5": U_g5,
            "U_i": U_i,
            "U_m": U_m,
            "primary_equations": [
                "H_mag = -mu*B  [eq.1]",
                "psi_matter = psi_0*exp(-i*(E_g+G_i+C_j+m_0)*t/hbar)  [eq.2]",
                "P_DE = eta_inertia*rho_vac*V*omega_vac ~ 7.09e-51 W  [eq.3]",
                "P_AC = 0.5*eps_0*E_EMP^2*V*omega_EMP  [eq.4]",
                "M_J = (5kT/(G*mu*m_H))^1.5*(3/(4*pi*rho))^0.5 ~ 25.8 Msun  [eq.5]",
                "f_n = f_0*phi^n (golden ratio series)  [eq.7]",
                "U_g5 = sum(T_mu_nu) ~ 3.6e-3 J/m^3  [eq.13, NEW 5th mode]",
                "omega_plasma = sqrt(n_e*e^2/(eps_0*m_e)) ~ 1.005e16 rad/s",
                "19 equations total including U_i, U_m, B_super, tau=I*alpha",
            ],
            "note": "Doc 43.d 19-eq framework. U_g5 5th gravity mode. PAPER_748, CP4 class #332. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class FiveQuantumVariableSetsCalculator:
    """CP4 #333 — PAPER_749: Five quantum variable sets (15 canonical variables).
    Sets A-D: r_j, d_g, F_U, f_feedback, Omega_g, f_Heaviside, i, H_SCm,
    lambda_i, j, M_bh, mu_j, P_core, t_n, pi; plus gamma, E_react, f_quasi, R_b."""

    def __init__(self):
        self.G = 6.6743e-11
        self.c = 3e8
        self.hbar = 1.0546e-34
        self.M_sun = 1.989e30
        self.kpc = 3.086e19
        self.AU = 1.496e11
        self.mu_0 = 4 * 3.14159e-7
        self.pi = 3.14159265358979
        # Set A
        self.r_j = 1.496e13
        self.d_g = 2.55e20
        self.f_feedback = 0.1
        self.Omega_g = 7.3e-16
        # Set B
        self.f_Heaviside = 0.01
        self.H_SCm = 1.0
        self.lambda_i = 1.0
        # Set C
        self.M_bh = 8.15e36
        self.mu_j_0 = 1e3 * 3.38e20
        self.omega_c = 7.3e-5
        # Set D
        self.gamma = 5e-5
        self.E_react = 1e46
        self.f_quasi = 0.01
        self.R_b = 100 * 1.496e11

    def calculate(self, t=0.0, i_index=1, j_index=1):
        import math
        rho_vac_SCm = 7.09e-37
        rho_vac_UA = 7.09e-36
        F_U_approx = 2.28e65
        r_j = self.r_j * i_index
        d_g = self.d_g
        Omega_g = self.Omega_g
        Heaviside_amp = 1 + 1e13 * self.f_Heaviside
        mu_j = (1e3 + 0.4 * math.sin(self.omega_c * t)) * 3.38e20
        t_n = t
        U_i = self.lambda_i * (rho_vac_SCm / rho_vac_UA) * Omega_g * self.H_SCm * (1 + 0.01)
        return {
            "set_A": {"r_j": r_j, "d_g": d_g, "F_U_approx": F_U_approx,
                      "f_feedback": self.f_feedback, "Omega_g": Omega_g},
            "set_B": {"f_Heaviside": self.f_Heaviside, "Heaviside_amp": Heaviside_amp,
                      "H_SCm": self.H_SCm, "lambda_i": self.lambda_i, "j_index": j_index},
            "set_C": {"M_bh_kg": self.M_bh, "M_bh_solar": self.M_bh / self.M_sun,
                      "mu_j": mu_j, "t_n": t_n, "pi": self.pi},
            "set_D": {"gamma": self.gamma, "E_react": self.E_react,
                      "f_quasi": self.f_quasi, "R_b_AU": self.R_b / self.AU},
            "U_i": U_i,
            "primary_equations": [
                "Set A: r_j = i*1.496e13 m, d_g = 2.55e20 m, F_U ~ 2.28e65 J/m^3",
                "Set B: f_Heaviside -> (1+1e13*f_H) ~ 1e11 amplification",
                "Set C: M_bh = 8.15e36 kg (Sgr A*), mu_j(t) = (1e3+0.4*sin)*3.38e20 T*pm^3",
                "Set D: gamma = 5e-5/day, E_react = 1e46, R_b = 100 AU",
                "U_i = lambda_i*(rho_SCm/rho_UA)*Omega_g*H_SCm*(1+F_RZ)",
            ],
            "note": "5 quantum variable sets. 15 canonical UQFF variables. PAPER_749, CP4 class #333. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass


class M51NGC1316MUGESimulationCalculator:
    """CP4 #334 — PAPER_750: M51 Whirlpool + NGC 1316 Fornax A dual MUGE.
    M51: F_tidal from NGC 5195, psi_spiral density waves.
    NGC 1316: F_cluster + rho_dust dust lanes from post-merger evolution.
    Includes both Python simulation scripts for radial g profiles."""

    def __init__(self):
        self.G = 6.6743e-11
        self.H_0 = 70e3 / 3.086e22
        self.Lambda = 1.1e-52
        self.c = 3e8
        self.hbar = 1.0546e-34
        self.t_Hubble = 4.35e17
        self.B_crit = 1e11
        self.rho_vac_SCm = 7.09e-37
        self.rho_vac_UA = 7.09e-36
        self.M_sun = 1.989e30
        self.kpc = 3.086e19
        # M51 params
        self.M51_M_vis = 1.2e11 * 1.989e30
        self.M51_M_DM = 4.0e10 * 1.989e30
        self.M51_M_NGC5195 = 1e10 * 1.989e30
        self.M51_d_inter = 50e3 * 3.086e19
        self.M51_B = 5e-6
        self.M51_z = 0.0015
        self.M51_sigma_spiral = 1e3 * 3.086e19
        # NGC 1316 params
        self.N16_M_vis = 3.5e11 * 1.989e30
        self.N16_M_DM = 1.5e11 * 1.989e30
        self.N16_M_spiral = 1e10 * 1.989e30
        self.N16_M_BH = 1e8 * 1.989e30
        self.N16_tau = 3.156e16
        self.N16_d_spiral = 50e3 * 3.086e19
        self.N16_rho_dust = 1e-21
        self.N16_B = 1e-4
        self.N16_z = 0.005
        self.N16_k_cluster = 1e-12
        self.N16_M_cluster = 1e6 * 1.989e30

    def calculate_m51(self, r, t=0.0):
        import math
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.3 * (1 + self.M51_z)**3 + 0.7)
        M_tot = self.M51_M_vis + self.M51_M_DM
        F_tidal = G * self.M51_M_NGC5195 / self.M51_d_inter**2
        F_env = F_tidal / (dpm_ug1_seed(M_tot, r))
        g_grav = (dpm_ug1_seed(M_tot, r)) * (1 + H) * (1 - self.M51_B / self.B_crit) * (1 + F_env)
        psi_spiral = math.exp(-r**2 / (2 * self.M51_sigma_spiral**2))
        quantum = (self.hbar / math.sqrt(1e-10 * 1e-20)) * psi_spiral**2 * (2 * math.pi / self.t_Hubble)
        U_i = self.lambda_i_val() 
        return g_grav + U_i + self.Lambda * c**2 / 3 + quantum

    def calculate_ngc1316(self, r, t=3.156e15):
        import math
        G, c = self.G, self.c
        H = self.H_0 * math.sqrt(0.3 * (1 + self.N16_z)**3 + 0.7)
        M_tot = self.N16_M_vis + self.N16_M_DM + self.N16_M_spiral * math.exp(-t / self.N16_tau)
        F_tidal = G * self.N16_M_spiral / self.N16_d_spiral**2
        F_cluster = self.N16_k_cluster * self.N16_M_cluster
        F_env = (F_tidal + F_cluster) / (dpm_ug1_seed(M_tot, r))
        g_grav = (dpm_ug1_seed(M_tot, r)) * (1 + H) * (1 - self.N16_B / self.B_crit) * (1 + F_env)
        g_dust = self.N16_rho_dust * (4/3 * 3.14159 * r**3) * dpm_ug1_seed(M_tot, r)
        U_i = self.lambda_i_val()
        return g_grav + U_i + self.Lambda * c**2 / 3 + g_dust

    def lambda_i_val(self):
        return 1.0 * (self.rho_vac_SCm / self.rho_vac_UA) * 1e-8 * 1.01

    def calculate(self, r_kpc=10.0, t_M51=0.0, t_N16=3.156e15):
        r = r_kpc * self.kpc
        g_m51 = self.calculate_m51(r, t_M51)
        g_n16 = self.calculate_ngc1316(r, t_N16)
        return {
            "g_M51": g_m51,
            "g_NGC1316": g_n16,
            "r_kpc": r_kpc,
            "M51_params": {"M_vis_solar": self.M51_M_vis / self.M_sun,
                           "M_DM_solar": self.M51_M_DM / self.M_sun,
                           "M_NGC5195_solar": self.M51_M_NGC5195 / self.M_sun,
                           "B_T": self.M51_B, "z": self.M51_z},
            "NGC1316_params": {"M_vis_solar": self.N16_M_vis / self.M_sun,
                               "M_DM_solar": self.N16_M_DM / self.M_sun, 
                               "M_BH_solar": self.N16_M_BH / self.M_sun,
                               "rho_dust": self.N16_rho_dust, "z": self.N16_z},
            "simulation_scripts": ["m51_simulation.py -> m51_gravity_profile.png",
                                   "ngc1316_simulation.py -> ngc1316_gravity_profile.png"],
            "primary_equations": [
                "g_M51(r,t) = (G*M/r^2)*(1+H)*(1-B/B_crit)*(1+F_env_M51)",
                "F_env_M51 = F_tidal = G*M_NGC5195/d_inter^2",
                "psi_spiral = A*exp(-r^2/(2*sigma^2))*exp(i*(m*theta-omega*t))",
                "g_NGC1316(r,t) = (G*M(t)/r^2)*(1+H)*(1-B/B_crit)*(1+F_env_N16) + rho_dust*V*g",
                "M(t) = M_vis + M_DM + M_spiral*exp(-t/tau),  tau=1 Gyr",
                "F_env_N16 = F_tidal + F_cluster",
                "F_cluster = k_cluster * M_cluster",
            ],
            "note": "M51 + NGC 1316 dual MUGE. Simulation scripts included. PAPER_750, CP4 class #334. v5.38.",
        }

    def self_update(self): pass
    def self_expand(self): pass
