from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell
"""
Append CP4 classes #335–#344 to CondensedPhysics4.py
Session 181 | PAPER_751–760 | v5.39
"""

BLOCK = '''

# ========================================================================
# CP4 #335 — THzQScopeEarthCoreSig1to50Calculator
# PAPER_751 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class THzQScopeEarthCoreSig1to50Calculator:
    ENTRY = 335
    PAPER = "PAPER_751"
    CPP_MODULE = "THzQScopeEarthCoreSig1to50"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("f_sampling_Hz", "float", 1.246, "Displayed sampling frequency (Hz)"),
        ("f_THz", "float", 1.25e12, "Actual THz resonance frequency (Hz)"),
        ("dA_A", "float", 6.205, "Full-scale current amplitude (A)"),
        ("V_pp_max_V", "float", 1.00, "Peak-to-peak voltage max (V)"),
        ("V_eff_V", "float", 0.35, "Effective (RMS) voltage (V)"),
        ("Z_imp_ohm", "float", 50.0, "Instrument impedance (Ohm)"),
        ("N_channels", "int", 50, "Number of resonance channels"),
    ]
    PRIMARY_OUTPUT = "P_peak_W"
    PRIMARY_INPUT = "V_eff_V"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        V_eff = p.get("V_eff_V", 0.35)
        Z_imp = p.get("Z_imp_ohm", 50.0)
        f_THz = p.get("f_THz", 1.25e12)
        dA    = p.get("dA_A", 6.205)
        N     = p.get("N_channels", 50)
        omega = 2 * math.pi * f_THz
        P_peak = V_eff**2 / Z_imp
        I_eff  = V_eff / Z_imp
        # UQFF resonance ratio [{U_m:SM_m}/Ug1^SCm]
        SCm = 0.99
        rho_UA, rho_SCm = self.UQFF_CONSTANTS["RHO_UA"], self.UQFF_CONSTANTS["RHO_SCM"]
        resonance_ratio = (rho_UA / rho_SCm) ** SCm
        P_coupled = P_peak * resonance_ratio * N
        return {
            "omega_THz_rad_s": omega,
            "P_peak_W": P_peak,
            "I_eff_A": I_eff,
            "dA_full_scale_A": dA,
            "resonance_ratio": resonance_ratio,
            "P_coupled_all_channels_W": P_coupled,
            "N_channels": N,
            "primary_equations": [
                "omega = 2*pi*f_THz = 7.854e12 rad/s",
                "P_peak = V_eff^2 / Z_imp = 0.35^2/50 = 2.45e-3 W",
                "P_THz = [{U_m:SM_m}/Ug1^SCm] * V_eff^2/Z_imp * sum_k delta(f-f_k)",
                "resonance_ratio = (rho_UA/rho_SCm)^SCm",
            ],
            "note": "THz hole between two pseudo-monopoles. PAPER_751, CP4 #335. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        import math
        results = []
        channels = sweep if sweep else [i * 2.5e9 + 1.2e12 for i in range(50)]
        for f_k in channels:
            omega_k = 2 * math.pi * f_k
            P_k = (0.35**2 / 50.0) * (1 + self.UQFF_CONSTANTS["F_TRZ"])
            results.append({"f_k_Hz": f_k, "omega_k": omega_k, "P_k_W": P_k})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #336 — V838MonLightEchoUQFFCalculator
# PAPER_752 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class V838MonLightEchoUQFFCalculator:
    ENTRY = 336
    PAPER = "PAPER_752"
    CPP_MODULE = "V838MonLightEchoUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("L_outburst_W", "float", 2.3e38, "Outburst luminosity (W)"),
        ("rho_0_kg_m3", "float", 1e-22, "Reference dust density (kg/m3)"),
        ("beta", "float", 1.0, "Ug1 attenuation factor"),
        ("sigma_scatter_m2", "float", 1e-27, "Scatter cross-section (m2)"),
        ("t_years", "float", 3.0, "Observation epoch (years)"),
    ]
    PRIMARY_OUTPUT = "I_echo_W_m2"
    PRIMARY_INPUT = "t_years"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        c = 3e8
        G = 6.674e-11
        M_star = 1.989e30
        L  = p.get("L_outburst_W", 2.3e38)
        r0 = p.get("rho_0_kg_m3", 1e-22)
        beta = p.get("beta", 1.0)
        sigma = p.get("sigma_scatter_m2", 1e-27)
        t_yr = p.get("t_years", 3.0)
        f_TRZ = self.UQFF_CONSTANTS["F_TRZ"]
        rho_UA  = self.UQFF_CONSTANTS["RHO_UA"]
        rho_SCm = self.UQFF_CONSTANTS["RHO_SCM"]
        t_s = t_yr * 3.156e7
        r_echo = c * t_s
        Ug1 = dpm_ug1_seed(M_star, r_echo)
        rho_dust = r0 * math.exp(-beta * Ug1)
        vac_correction = 1 + rho_UA / rho_SCm
        I_echo = (L / (4 * math.pi * r_echo**2)) * sigma * rho_dust * (1 + f_TRZ) * vac_correction
        return {
            "r_echo_m": r_echo,
            "Ug1_m_s2": Ug1,
            "rho_dust_kg_m3": rho_dust,
            "I_echo_W_m2": I_echo,
            "vac_correction_factor": vac_correction,
            "primary_equations": [
                "r_echo(t) = c*t",
                "rho_dust(r,t) = rho_0 * exp(-beta*Ug1(r,t))",
                "I_echo = [L/(4*pi*(c*t)^2)] * sigma * rho_dust * (1+fTRZ) * (1+rho_UA/rho_SCm)",
            ],
            "note": "V838 Mon UQFF light echo intensity. PAPER_752, CP4 #336. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) for i in range(1, 11)]
        for t_yr in t_vals:
            res = self.compute({"t_years": t_yr})
            results.append({"t_yr": t_yr, "I_echo": res["I_echo_W_m2"], "r_echo_m": res["r_echo_m"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #337 — MagnetarEvolutionUQFFCalculator
# PAPER_753 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class MagnetarEvolutionUQFFCalculator:
    ENTRY = 337
    PAPER = "PAPER_753"
    CPP_MODULE = "MagnetarEvolutionUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_kg", "float", 2.785e30, "Neutron star mass (kg)"),
        ("r_m", "float", 2.0e4, "Neutron star radius (m)"),
        ("B0_T", "float", 1e10, "Initial B-field (T)"),
        ("tau_B_s", "float", 1.262e11, "B-field decay timescale (s)"),
        ("B_crit_T", "float", 1e11, "Critical B-field (T)"),
        ("Omega0_rad_s", "float", 1.2566, "Initial spin rate (rad/s)"),
        ("tau_spin_s", "float", 3.156e11, "Spin-down timescale (s)"),
        ("t_s", "float", 1.578e11, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_Magnetar_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G, c = 6.674e-11, 3e8
        H0 = 2.184e-18
        M  = p.get("M_kg", 2.785e30)
        r  = p.get("r_m", 2.0e4)
        B0 = p.get("B0_T", 1e10)
        tau_B = p.get("tau_B_s", 1.262e11)
        B_crit = p.get("B_crit_T", 1e11)
        Omega0 = p.get("Omega0_rad_s", 1.2566)
        tau_spin = p.get("tau_spin_s", 3.156e11)
        t  = p.get("t_s", 1.578e11)
        B_t = B0 * math.exp(-t / tau_B)
        Omega_t = Omega0 * math.exp(-t / tau_spin)
        g_grav = (dpm_ug1_seed(M, r)) * (1 + H0 * t) * (1 - B_t / B_crit)
        # Ug1 + Ug4 floor
        Ug_floor = 1.007e12
        # GW quadrupole spin-down contribution
        GW_term = (32 * G**4 * M**3 * r**2 * Omega_t**4) / (5 * c**5 * r**4)
        g_total = g_grav + Ug_floor + GW_term
        return {
            "B_t_T": B_t,
            "Omega_t_rad_s": Omega_t,
            "g_grav_m_s2": g_grav,
            "GW_term": GW_term,
            "g_Magnetar_m_s2": g_total,
            "primary_equations": [
                "g(t) = (G*M/r^2)*(1+H0*t)*(1-B(t)/B_crit) + Ug_floor + GW_term",
                "B(t) = B0*exp(-t/tau_B)",
                "Omega(t) = Omega0*exp(-t/tau_spin)",
                "GW_term = 32*G^4*M^3*r^2*Omega^4/(5*c^5*r^4)",
            ],
            "note": "g_Magnetar(t=5kyr) ~ 4.474e12 m/s^2. PAPER_753, CP4 #337. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 3.156e10 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_Magnetar": res["g_Magnetar_m_s2"], "B_t": res["B_t_T"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #338 — SgrAStarEvolutionUQFFCalculator
# PAPER_754 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class SgrAStarEvolutionUQFFCalculator:
    ENTRY = 338
    PAPER = "PAPER_754"
    CPP_MODULE = "SgrAStarEvolutionUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_kg", "float", 8.552e36, "SMBH mass (kg)"),
        ("r_s_m", "float", 1.27e10, "Schwarzschild radius (m)"),
        ("M_dot0_solar_yr", "float", 0.01, "Initial accretion rate (M_sun/yr)"),
        ("tau_acc_s", "float", 2.84e17, "Accretion timescale (s)"),
        ("tau_spin_s", "float", 2.84e17, "Spin-down timescale (s)"),
        ("theta_prec_deg", "float", 30.0, "Precession angle (deg)"),
        ("t_s", "float", 1.42e17, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_SgrA_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G, c = 6.674e-11, 3e8
        H0, Lambda = 2.184e-18, 1.089e-52
        M_sun = 1.989e30
        M    = p.get("M_kg", 8.552e36)
        r_s  = p.get("r_s_m", 1.27e10)
        Mdot0 = p.get("M_dot0_solar_yr", 0.01) * M_sun / 3.156e7
        tau_acc = p.get("tau_acc_s", 2.84e17)
        tau_sp  = p.get("tau_spin_s", 2.84e17)
        theta   = math.radians(p.get("theta_prec_deg", 30.0))
        t = p.get("t_s", 1.42e17)
        Mdot_t = Mdot0 * math.exp(-t / tau_acc)
        # Ω(t) = (0.3c/r_s)*exp(-t/tau_sp)
        Omega_t = (0.3 * c / r_s) * math.exp(-t / tau_sp)
        g_grav  = (dpm_ug1_seed(M, r_s)) * (1 + H0 * t) * math.sin(theta)
        g_Ug    = (Lambda * c**2 / 3) + 1.0e3   # Ug1+Ug2 small floor
        g_total = g_grav + g_Ug
        return {
            "Mdot_t_kg_s": Mdot_t,
            "Omega_t_rad_s": Omega_t,
            "g_grav_m_s2": g_grav,
            "g_SgrA_m_s2": g_total,
            "primary_equations": [
                "g(r_s,t) = (G*M/r_s^2)*(1+H0*t)*sin(theta_prec) + Ug_floor",
                "M_dot(t) = M_dot0*exp(-t/tau_acc)",
                "Omega(t) = (0.3c/r_s)*exp(-t/tau_spin)",
                "r_s = 2GM/c^2 = 1.27e10 m",
            ],
            "note": "g_SgrA*(t=4.5Gyr) ~ 1.250e7 m/s^2. PAPER_754, CP4 #338. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 3.156e16 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_SgrA": res["g_SgrA_m_s2"], "Mdot": res["Mdot_t_kg_s"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #339 — TapestryBlazingStarbirthNGC2014Calculator
# PAPER_755 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class TapestryBlazingStarbirthNGC2014Calculator:
    ENTRY = 339
    PAPER = "PAPER_755"
    CPP_MODULE = "TapestryBlazingStarbirthNGC2014"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_initial_kg", "float", 4.774e32, "Initial cluster mass (kg)"),
        ("r_m", "float", 9.461e16, "Cluster radius (m)"),
        ("v_wind_mps", "float", 2e6, "Wind velocity (m/s)"),
        ("rho_ISM", "float", 1e-21, "ISM density (kg/m3)"),
        ("M_dot0_solar_yr", "float", 41.67, "Mean SF accretion rate (M_sun/yr)"),
        ("tau_SF_s", "float", 1.578e14, "SF timescale (s)"),
        ("B_T", "float", 1e-6, "Magnetic field (T)"),
        ("t_s", "float", 7.89e13, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_Starbirth_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.674e-11
        M_sun = 1.989e30
        M0  = p.get("M_initial_kg", 4.774e32)
        r   = p.get("r_m", 9.461e16)
        v_w = p.get("v_wind_mps", 2e6)
        rho = p.get("rho_ISM", 1e-21)
        Md0 = p.get("M_dot0_solar_yr", 41.67) * M_sun / 3.156e7
        tau = p.get("tau_SF_s", 1.578e14)
        B   = p.get("B_T", 1e-6)
        t   = p.get("t_s", 7.89e13)
        # Star formation mass loading
        M_SF = Md0 * t * math.exp(-t / tau)
        M_t  = M0 + M_SF
        H0 = 2.184e-18
        B_crit = 4.4e9
        g_grav  = (dpm_ug1_seed(M_t, r)) * (1 + H0 * t) * (1 - B / B_crit)
        g_ram   = rho * v_w**2 / r
        # EM Aether term: q*(v*B)*11e-12
        g_EM    = 1.0 * (v_w * B) * 11 * 1e-12
        g_total = g_grav + g_ram + g_EM
        return {
            "M_t_kg": M_t,
            "g_grav_m_s2": g_grav,
            "g_ram_m_s2": g_ram,
            "g_EM_m_s2": g_EM,
            "g_Starbirth_m_s2": g_total,
            "primary_equations": [
                "g(r,t) = (G*M(t)/r^2)*(1+H0*t)*(1-B/Bcrit) + rho*v^2/r + q*(v*B)*11e-12",
                "M(t) = M_initial*(1 + M_SF(t))",
                "g_EM = q*(v_wind x B) * A_aeth * A_scale  [dominant term]",
            ],
            "note": "g_Starbirth(t=2.5Myr) ~ 1.053e-4 m/s^2. PAPER_755, CP4 #339. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 1.578e13 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_total": res["g_Starbirth_m_s2"], "g_EM": res["g_EM_m_s2"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #340 — Westerlund2SuperClusterUQFFCalculator
# PAPER_756 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class Westerlund2SuperClusterUQFFCalculator:
    ENTRY = 340
    PAPER = "PAPER_756"
    CPP_MODULE = "Westerlund2SuperClusterUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_initial_kg", "float", 5.967e34, "Initial cluster mass (kg)"),
        ("r_m", "float", 9.461e16, "Cluster radius (m)"),
        ("v_wind_mps", "float", 2e6, "Wind velocity (m/s)"),
        ("rho_ISM", "float", 1e-20, "ISM density (kg/m3)"),
        ("M_dot0_solar_yr", "float", 3.333, "Mean SF accretion rate (M_sun/yr)"),
        ("tau_SF_s", "float", 6.312e13, "SF timescale (s)"),
        ("B_T", "float", 1e-5, "Magnetic field (T)"),
        ("t_s", "float", 3.156e13, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_Westerlund2_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.674e-11
        M_sun = 1.989e30
        M0  = p.get("M_initial_kg", 5.967e34)
        r   = p.get("r_m", 9.461e16)
        v_w = p.get("v_wind_mps", 2e6)
        rho = p.get("rho_ISM", 1e-20)
        Md0 = p.get("M_dot0_solar_yr", 3.333) * M_sun / 3.156e7
        tau = p.get("tau_SF_s", 6.312e13)
        B   = p.get("B_T", 1e-5)
        t   = p.get("t_s", 3.156e13)
        M_SF = Md0 * t * math.exp(-t / tau)
        M_t  = M0 + M_SF
        H0   = 2.184e-18
        B_crit = 4.4e9
        g_grav  = (dpm_ug1_seed(M_t, r)) * (1 + H0 * t) * (1 - B / B_crit)
        g_ram   = rho * v_w**2 / r
        g_EM    = 1.0 * (v_w * B) * 11 * 1e-12
        g_total = g_grav + g_ram + g_EM
        return {
            "M_t_kg": M_t,
            "g_grav_m_s2": g_grav,
            "g_ram_m_s2": g_ram,
            "g_EM_m_s2": g_EM,
            "g_Westerlund2_m_s2": g_total,
            "primary_equations": [
                "g(r,t) = (G*M(t)/r^2)*(1+H0*t)*(1-B/Bcrit) + rho*v^2/r + q*(v*B)*11e-12",
                "rho_ISM=1e-20 kg/m3 (10x NGC2014), B=1e-5 T -> g_EM x10 larger",
            ],
            "note": "g_W2(t=1Myr) ~ 1.053e-3 m/s^2. PAPER_756, CP4 #340. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 6.312e12 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_W2": res["g_Westerlund2_m_s2"], "g_EM": res["g_EM_m_s2"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #341 — PillarsOfCreationM16ErosionCalculator
# PAPER_757 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class PillarsOfCreationM16ErosionCalculator:
    ENTRY = 341
    PAPER = "PAPER_757"
    CPP_MODULE = "PillarsOfCreationM16Erosion"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_kg", "float", 2.009e34, "Total pillar mass (kg)"),
        ("r_m", "float", 4.731e16, "Pillar half-length (m)"),
        ("rho_ISM", "float", 1e-21, "ISM density (kg/m3)"),
        ("B_T", "float", 1e-6, "Magnetic field (T)"),
        ("v_wind_mps", "float", 2e6, "Wind velocity (m/s)"),
        ("E0", "float", 0.1, "Erosion amplitude"),
        ("tau_erode_s", "float", 3.156e13, "Erosion timescale (s)"),
        ("t_s", "float", 1.578e13, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_Pillars_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.674e-11
        M  = p.get("M_kg", 2.009e34)
        r  = p.get("r_m", 4.731e16)
        rho = p.get("rho_ISM", 1e-21)
        B   = p.get("B_T", 1e-6)
        v_w = p.get("v_wind_mps", 2e6)
        E0  = p.get("E0", 0.1)
        tau = p.get("tau_erode_s", 3.156e13)
        t   = p.get("t_s", 1.578e13)
        H0 = 2.184e-18
        B_crit = 4.4e9
        E_t = E0 * math.exp(-t / tau)
        surv = 1.0 - E_t
        g_grav = (dpm_ug1_seed(M, r)) * (1 + H0 * t) * (1 - B / B_crit) * surv
        g_ram  = rho * v_w**2 / r
        g_EM   = 1.0 * (v_w * B) * 11 * 1e-12 * surv
        g_total = g_grav + g_ram + g_EM
        return {
            "E_t": E_t,
            "survival_factor": surv,
            "g_grav_m_s2": g_grav,
            "g_EM_m_s2": g_EM,
            "g_Pillars_m_s2": g_total,
            "primary_equations": [
                "g(r,t) = (G*M/r^2)*(1+H0*t)*(1-B/Bcrit)*(1-E(t)) + ram + g_EM*(1-E(t))",
                "E(t) = E0*exp(-t/tau_erode)",
                "(1-E(t=0.5Myr)) ~ 0.93935",
            ],
            "note": "g_Pillars(t=0.5Myr) ~ 1.053e-4 m/s^2. PAPER_757, CP4 #341. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 3.156e12 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_Pillars": res["g_Pillars_m_s2"], "E_t": res["E_t"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #342 — RingsOfRelativityEinsteinRingCalculator
# PAPER_758 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class RingsOfRelativityEinsteinRingCalculator:
    ENTRY = 342
    PAPER = "PAPER_758"
    CPP_MODULE = "RingsOfRelativityEinsteinRing"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_kg", "float", 1.989e44, "Cluster mass (kg)"),
        ("r_E_m", "float", 3.086e20, "Einstein radius (m)"),
        ("z_lens", "float", 0.5, "Lens redshift"),
        ("D_LS_over_DS", "float", 0.5, "D_LS/D_S ratio"),
        ("B_ICM_T", "float", 1e-5, "ICM magnetic field (T)"),
        ("v_ICM_mps", "float", 1e6, "ICM velocity (m/s)"),
        ("t_Gyr", "float", 5.0, "Evaluation epoch (Gyr)"),
    ]
    PRIMARY_OUTPUT = "g_Rings_m_s2"
    PRIMARY_INPUT = "t_Gyr"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G, c = 6.674e-11, 3e8
        H0_km = 70.0
        Omega_m, Omega_L = 0.3, 0.7
        M   = p.get("M_kg", 1.989e44)
        r_E = p.get("r_E_m", 3.086e20)
        z   = p.get("z_lens", 0.5)
        DLS_DS = p.get("D_LS_over_DS", 0.5)
        B   = p.get("B_ICM_T", 1e-5)
        v   = p.get("v_ICM_mps", 1e6)
        t_s = p.get("t_Gyr", 5.0) * 3.156e16
        B_crit = 4.4e9
        # Hubble at z
        H_z = H0_km * math.sqrt(Omega_m * (1 + z)**3 + Omega_L)  # km/s/Mpc
        H_z_si = H_z * 1e3 / 3.086e22  # s^-1
        # Lensing efficiency
        L_t = (G * M / (c**2 * r_E)) * DLS_DS
        g_grav = (dpm_ug1_seed(M, r_E)) * (1 + H_z_si * t_s) * (1 - B / B_crit) * (1 + L_t)
        g_EM   = 1.0 * (v * B) * 11 * 1e-12
        g_total = g_grav + g_EM
        return {
            "H_z_km_s_Mpc": H_z,
            "L_t_lensing": L_t,
            "g_grav_m_s2": g_grav,
            "g_EM_m_s2": g_EM,
            "g_Rings_m_s2": g_total,
            "primary_equations": [
                "g(r_E,t) = (G*M/r_E^2)*(1+H(z)*t)*(1-B/Bcrit)*(1+L(t)) + g_EM",
                "H(z=0.5) = 70*sqrt(0.3*3.375+0.7) = 91.63 km/s/Mpc",
                "L(t) = G*M/(c^2*r_E) * D_LS/D_S = 2.388e-4",
            ],
            "note": "g_Rings ~ 1.053e-2 m/s^2. PAPER_758, CP4 #342. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        z_vals = sweep if sweep else [i * 0.1 for i in range(1, 11)]
        for z in z_vals:
            res = self.compute({"z_lens": z})
            results.append({"z_lens": z, "g_Rings": res["g_Rings_m_s2"], "H_z": res["H_z_km_s_Mpc"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #343 — HorseheadNebulaBarnard33UQFFCalculator
# PAPER_759 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class HorseheadNebulaBarnard33UQFFCalculator:
    ENTRY = 343
    PAPER = "PAPER_759"
    CPP_MODULE = "HorseheadNebulaBarnard33UQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("r_m", "float", 1.182e16, "Pillar tip radius (m)"),
        ("L_star_W", "float", 3.826e31, "Ionising star luminosity (W)"),
        ("rho_ISM", "float", 1e-21, "ISM density (kg/m3)"),
        ("B_T", "float", 1e-5, "Magnetic field (T)"),
        ("v_mps", "float", 1e5, "Velocity (m/s)"),
        ("E0", "float", 0.1, "Erosion amplitude"),
        ("tau_erode_s", "float", 3.156e13, "Erosion timescale (s)"),
        ("t_s", "float", 1.578e13, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_Horsehead_m_s2"
    PRIMARY_INPUT = "r_m"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        c = 3e8
        G = 6.674e-11
        m_H = 1.67e-27
        M_cloud = 1.989e30 * 500  # ~500 M_sun Barnard 33 core
        r   = p.get("r_m", 1.182e16)
        L_s = p.get("L_star_W", 3.826e31)
        rho = p.get("rho_ISM", 1e-21)
        B   = p.get("B_T", 1e-5)
        v   = p.get("v_mps", 1e5)
        E0  = p.get("E0", 0.1)
        tau = p.get("tau_erode_s", 3.156e13)
        t   = p.get("t_s", 1.578e13)
        H0  = 2.184e-18
        B_crit = 4.4e9
        E_t  = E0 * math.exp(-t / tau)
        surv = 1.0 - E_t
        # Radiation pressure acceleration
        P_rad = (L_s / (4 * math.pi * r**2 * c)) * (rho / m_H)
        # UQFF gravity
        g_grav = (dpm_ug1_seed(M_cloud, r)) * (1 + H0 * t) * (1 - B / B_crit) * surv
        # EM Aether term
        g_EM   = 1.0 * (v * B) * 11 * 1e-12 * surv
        g_total = g_grav + P_rad + g_EM
        return {
            "E_t": E_t,
            "survival_factor": surv,
            "P_rad_m_s2": P_rad,
            "g_grav_m_s2": g_grav,
            "g_EM_m_s2": g_EM,
            "g_Horsehead_m_s2": g_total,
            "primary_equations": [
                "g(r,t) = (G*M/r^2)*(1+H0*t)*(1-B/Bcrit)*(1-E(t)) + P_rad + g_EM*(1-E(t))",
                "P_rad = L_star/(4*pi*r^2*c) * rho/m_H = 4.347e-5 m/s^2",
                "g_EM ~ 1.097e-3 m/s^2 [dominant]",
            ],
            "note": "g_Horsehead ~ 1.097e-3 m/s^2. PAPER_759, CP4 #343. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 3.156e12 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_HH": res["g_Horsehead_m_s2"], "P_rad": res["P_rad_m_s2"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)


# ========================================================================
# CP4 #344 — NGC1275MagneticMonsterPerseusACalculator
# PAPER_760 | Session 181 | thread_06Jun2025.txt
# ========================================================================
class NGC1275MagneticMonsterPerseusACalculator:
    ENTRY = 344
    PAPER = "PAPER_760"
    CPP_MODULE = "NGC1275MagneticMonsterPerseusA"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_total_kg", "float", 1.991e42, "Total cluster mass (kg)"),
        ("r_m", "float", 9.46e20, "Evaluation radius (m)"),
        ("M_SMBH_kg", "float", 1.592e39, "SMBH mass (kg)"),
        ("tau_BH_s", "float", 3.156e15, "AGN feedback timescale (s)"),
        ("B_fil_T", "float", 1e-8, "Filament B-field (T)"),
        ("V_fil_m3", "float", 1.42e50, "Filament volume (m3)"),
        ("M_fil_kg", "float", 1.989e36, "Filament mass (kg)"),
        ("v_merger_mps", "float", 3e6, "Merger velocity (m/s)"),
        ("z", "float", 0.0176, "Cluster redshift"),
        ("t_s", "float", 1.578e15, "Evaluation epoch (s)"),
    ]
    PRIMARY_OUTPUT = "g_NGC1275_m_s2"
    PRIMARY_INPUT = "t_s"

    def compute(self, params=None):
        import math
        p = {k: v for _, k, v, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.674e-11
        mu0 = 4 * math.pi * 1e-7
        H0_km = 70.0
        Omega_m, Omega_L = 0.3, 0.7
        M   = p.get("M_total_kg", 1.991e42)
        r   = p.get("r_m", 9.46e20)
        M_BH = p.get("M_SMBH_kg", 1.592e39)
        tau_BH = p.get("tau_BH_s", 3.156e15)
        B_fil = p.get("B_fil_T", 1e-8)
        V_fil = p.get("V_fil_m3", 1.42e50)
        M_fil = p.get("M_fil_kg", 1.989e36)
        v_m  = p.get("v_merger_mps", 3e6)
        z    = p.get("z", 0.0176)
        t    = p.get("t_s", 1.578e15)
        B_crit = 4.4e9
        F0 = 0.1
        F_BH = F0 * (1 - math.exp(-t / tau_BH))
        H_z_km = H0_km * math.sqrt(Omega_m * (1 + z)**3 + Omega_L)
        H_z_si = H_z_km * 1e3 / 3.086e22
        g_grav  = (dpm_ug1_seed(M, r)) * (1 + H_z_si * t) * (1 - B_fil / B_crit) * (1 - F_BH)
        # Filament magnetic support
        a_fil = (B_fil**2 * V_fil) / (2 * mu0 * M_fil * r)
        # EM Aether merger term
        g_EM  = 1.0 * (v_m * B_fil) * 11 * 1e-12
        g_total = g_grav + a_fil + g_EM
        return {
            "F_BH": F_BH,
            "H_z_km_s_Mpc": H_z_km,
            "g_grav_m_s2": g_grav,
            "a_fil_m_s2": a_fil,
            "g_EM_m_s2": g_EM,
            "g_NGC1275_m_s2": g_total,
            "primary_equations": [
                "g(r,t) = (G*M_total/r^2)*(1+H(z)*t)*(1-B_fil/Bcrit)*(1-F_BH(t)) + a_fil + g_EM",
                "F_BH(t) = F0*(1-exp(-t/tau_BH)) = 0.03935  at t=50Myr",
                "a_fil = B^2*V_fil/(2*mu0*M_fil*r) = 2.840e-9 m/s^2",
                "H(z=0.0176) = 70.56 km/s/Mpc",
            ],
            "note": "g_NGC1275(t=50Myr) ~ 3.160e-5 m/s^2. PAPER_760, CP4 #344. v5.39.",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        t_vals = sweep if sweep else [float(i) * 3.156e14 for i in range(1, 11)]
        for t in t_vals:
            res = self.compute({"t_s": t})
            results.append({"t_s": t, "g_N1275": res["g_NGC1275_m_s2"], "F_BH": res["F_BH"]})
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath) as f:
            data = json.load(f)
        for k, v in data.items():
            setattr(self, k, v)
'''

with open("CondensedPhysics4.py", "a", encoding="utf-8") as f:
    f.write(BLOCK)

print("Appended CP4 #335-344 to CondensedPhysics4.py")
