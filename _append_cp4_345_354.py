"""
_append_cp4_345_354.py
Session 181 — Batch 2: PAPER_761–770
Appends CP4 classes #345–354 to CondensedPhysics4.py
"""

BLOCK = '''
# ========================================================================
# CP4 #345 — HubbleUltraDeepFieldUQFFCalculator
# PAPER_761 | Session 181 | thread_06Jun2025.txt
# Hubble Ultra Deep Field — UQFF Galaxy Evolution
# g_HUDF ≈ 1.053e-3 m/s²
# ========================================================================
class HubbleUltraDeepFieldUQFFCalculator:
    ENTRY = 345
    PAPER = "PAPER_761"
    CPP_MODULE = "HubbleUltraDeepFieldUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e42, "Total HUDF field mass [kg] (10^12 M_sun)"),
        ("r", "float", 1.5e22, "HUDF field radius [m] (1.5 Mpc)"),
        ("z", "float", 3.0, "Average galaxy redshift"),
        ("t", "float", 4.103e17, "Lookback time [s] (13 Gyr)"),
        ("M_evo", "float", 0.13, "Evolutionary mass fraction"),
        ("M_merge", "float", 0.2, "Merger mass loss fraction"),
        ("v_EM", "float", 1e6, "Aether EM velocity [m/s]"),
        ("B_EM", "float", 1e-6, "Intergalactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_HUDF_m_s2"
    PRIMARY_INPUT = "M"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        c = 3e8
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_evo = 1.0 + p["M_evo"]
        factor_merge = 1.0 - p["M_merge"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_HUDF = g_grav * factor_Hz * factor_evo * factor_merge * factor_TRZ + a_EM
        return {
            "g_HUDF_m_s2": g_HUDF,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "Hz_s": Hz,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e41, 5e41, 1e42, 5e42, 1e43]
        for val in sweep:
            p = {sweep_param or "M": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #346 — NGC1792StellarForgeUQFFCalculator
# PAPER_762 | Session 181 | thread_06Jun2025.txt
# NGC 1792 "The Stellar Forge" — UQFF Starburst Evolution
# g_NGC1792 ≈ 1.053e-2 m/s²
# ========================================================================
class NGC1792StellarForgeUQFFCalculator:
    ENTRY = 346
    PAPER = "PAPER_762"
    CPP_MODULE = "NGC1792StellarForgeUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e40, "Galaxy mass [kg] (10^10 M_sun)"),
        ("r", "float", 3.78e20, "Galaxy radius [m] (~40 kly)"),
        ("z", "float", 0.0095, "Redshift"),
        ("t", "float", 3.156e15, "SFR integration time [s] (100 Myr)"),
        ("SFR", "float", 10.0, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.1, "Stellar mass fraction"),
        ("F_sn", "float", 0.031605, "Supernova feedback fraction"),
        ("v_EM", "float", 1e6, "Aether EM velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC1792_m_s2"
    PRIMARY_INPUT = "SFR"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_sn = 1.0 - p["F_sn"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_NGC1792 = g_grav * factor_Hz * factor_sf * factor_sn * factor_TRZ + a_EM
        return {
            "g_NGC1792_m_s2": g_NGC1792,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1.0, 5.0, 10.0, 20.0, 50.0]
        for val in sweep:
            p = {sweep_param or "SFR": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #347 — SombreroGalaxyM104UQFFCalculator
# PAPER_763 | Session 181 | thread_06Jun2025.txt
# Sombrero Galaxy M104 — UQFF SMBH + Dust Lane
# g_Sombrero ≈ 5.351e-1 m/s²
# ========================================================================
class SombreroGalaxyM104UQFFCalculator:
    ENTRY = 347
    PAPER = "PAPER_763"
    CPP_MODULE = "SombreroGalaxyM104UQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 2.009e41, "Total galaxy mass [kg]"),
        ("r", "float", 2.36e20, "Galaxy radius [m] (25 kly)"),
        ("M_BH", "float", 1.989e39, "SMBH mass [kg] (10^9 M_sun)"),
        ("r_BH", "float", 1e15, "SMBH influence radius [m]"),
        ("z", "float", 0.0063, "Redshift"),
        ("t", "float", 8.086e16, "Age [s] (2.56 Gyr)"),
        ("a_dust", "float", 0.4, "Dust lane absorption acceleration [m/s^2]"),
        ("v_EM", "float", 2e5, "Aether EM velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_Sombrero_m_s2"
    PRIMARY_INPUT = "M_BH"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        g_BH = G * p["M_BH"] / p["r_BH"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_Sombrero = (g_grav * factor_Hz * factor_TRZ
                      + g_BH + p["a_dust"] + a_EM)
        return {
            "g_Sombrero_m_s2": g_Sombrero,
            "g_grav": g_grav,
            "g_BH": g_BH,
            "a_dust": p["a_dust"],
            "a_EM": a_EM,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e38, 5e38, 1e39, 5e39, 1e40]
        for val in sweep:
            p = {sweep_param or "M_BH": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #348 — Saturn26DUQFFCalculator
# PAPER_764 | Session 181 | thread_06Jun2025.txt
# Saturn Ring System 26D UQFF Planetary Evolution
# g_Saturn ≈ 10.44 m/s²
# ========================================================================
class Saturn26DUQFFCalculator:
    ENTRY = 348
    PAPER = "PAPER_764"
    CPP_MODULE = "Saturn26DUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_Saturn", "float", 5.683e26, "Saturn mass [kg]"),
        ("r_surface", "float", 6.0268e7, "Saturn equatorial radius [m]"),
        ("M_Sun", "float", 1.989e30, "Solar mass [kg]"),
        ("r_orbit", "float", 1.43e12, "Saturn orbital radius [m]"),
        ("M_ring", "float", 1.5e19, "Ring system mass [kg]"),
        ("r_ring", "float", 7e7, "Main ring outer radius [m]"),
        ("v_wind", "float", 400.0, "Saturn wind velocity [m/s]"),
        ("B_Saturn", "float", 2e-5, "Saturn surface B field [T]"),
        ("t", "float", 1.420e17, "System age [s] (4.5 Gyr)"),
    ]
    PRIMARY_OUTPUT = "g_Saturn_m_s2"
    PRIMARY_INPUT = "M_Saturn"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        z = 0.0
        Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        g_Saturn = G * p["M_Saturn"] / p["r_surface"]**2
        g_Sun_tide = G * p["M_Sun"] / p["r_orbit"]**2
        T_ring = G * p["M_ring"] / p["r_ring"]**2 * 1e-7
        a_wind = p["v_wind"]**2 / p["r_surface"] * 1e-6
        m_e = 9.11e-31
        M_mag = (q * p["v_wind"] * p["B_Saturn"] / m_e) * 1e-12
        g_total = g_Saturn * factor_Hz + g_Sun_tide * 1e-4 + T_ring + a_wind + M_mag
        return {
            "g_Saturn_m_s2": g_total,
            "g_Saturn_surface": g_Saturn,
            "g_Sun_tide": g_Sun_tide,
            "T_ring": T_ring,
            "a_wind": a_wind,
            "factor_Hz": factor_Hz,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e26, 5.683e26, 7e26, 1e27, 2e27]
        for val in sweep:
            p = {sweep_param or "M_Saturn": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #349 — M16EagleNebulaStarsUQFFCalculator
# PAPER_765 | Session 181 | thread_06Jun2025.txt
# M16 Eagle Nebula Pillars of Creation 26D UQFF Star Formation
# g_M16 ≈ 1.053e-3 m/s²
# ========================================================================
class M16EagleNebulaStarsUQFFCalculator:
    ENTRY = 349
    PAPER = "PAPER_765"
    CPP_MODULE = "M16EagleNebulaStarsUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 2.387e33, "Region mass [kg] (1200 M_sun)"),
        ("r", "float", 3.31e17, "Pillar region radius [m] (~35 ly)"),
        ("z", "float", 0.0015, "Redshift"),
        ("t", "float", 1.578e14, "SFR time [s] (5 Myr)"),
        ("M_sf_factor", "float", 4.472, "1+M_sf: stellar mass growth factor"),
        ("E_rad", "float", 0.2433, "Radiation energy loss fraction"),
        ("v_EM", "float", 1e5, "Aether EM velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Pillar B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_M16_m_s2"
    PRIMARY_INPUT = "M"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = p["M_sf_factor"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_M16 = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_M16_m_s2": g_M16,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1.0, 2.0, 3.0, 4.472, 6.0]
        for val in sweep:
            p = {sweep_param or "M_sf_factor": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #350 — CrabNebulaPulsarWindUQFFCalculator
# PAPER_766 | Session 181 | thread_06Jun2025.txt
# Crab Nebula Pulsar Wind 26D UQFF SNR
# g_Crab ≈ 1.481e6 m/s²
# ========================================================================
class CrabNebulaPulsarWindUQFFCalculator:
    ENTRY = 350
    PAPER = "PAPER_766"
    CPP_MODULE = "CrabNebulaPulsarWindUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 9.149e30, "Nebula mass [kg] (4.6 M_sun)"),
        ("r", "float", 5.2e16, "Nebula radius [m] (5.5 ly)"),
        ("P_pulsar", "float", 5e31, "Pulsar spin-down luminosity [W]"),
        ("v_shock", "float", 1.5e6, "Shock expansion velocity [m/s]"),
        ("rho_fil", "float", 1e-21, "Filament gas density [kg/m^3]"),
        ("B_mag", "float", 1e-8, "Average nebula B field [T]"),
        ("z", "float", 0.0015, "Redshift"),
        ("t", "float", 3.064e10, "Nebula age [s] (971 yr)"),
    ]
    PRIMARY_OUTPUT = "g_Crab_m_s2"
    PRIMARY_INPUT = "P_pulsar"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_e = 9.11e-31
        q = 1.602e-19
        c = 3e8
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        F_wind = (p["P_pulsar"] / (4 * math.pi * p["r"]**2)) * (1 + p["v_shock"] / c)
        a_wind = (F_wind / p["rho_fil"]) * 1e-12
        M_mag = (q * p["v_shock"] * p["B_mag"] / m_e) * 1e-12
        g_Crab = g_grav * factor_Hz * factor_TRZ + a_wind + M_mag
        return {
            "g_Crab_m_s2": g_Crab,
            "g_grav": g_grav,
            "a_wind": a_wind,
            "M_mag": M_mag,
            "factor_Hz": factor_Hz,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e31, 5e31, 1e32, 5e32, 1e33]
        for val in sweep:
            p = {sweep_param or "P_pulsar": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #351 — NGC2264ConeNebulaUQFFCalculator
# PAPER_767 | Session 181 | thread_06Jun2025.txt
# NGC 2264 Cone Nebula / Christmas Tree Cluster UQFF Star Formation
# g_NGC2264 ≈ 1.053e-2 m/s²
# ========================================================================
class NGC2264ConeNebulaUQFFCalculator:
    ENTRY = 351
    PAPER = "PAPER_767"
    CPP_MODULE = "NGC2264ConeNebulaUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e33, "Cluster mass [kg] (1000 M_sun)"),
        ("r", "float", 4.73e16, "Region radius [m] (~5 ly)"),
        ("z", "float", 0.0006, "Redshift"),
        ("t", "float", 9.468e13, "Integration time [s] (3 Myr)"),
        ("M_sf", "float", 1.5, "Star-formation mass ratio"),
        ("E_rad", "float", 0.1554, "Radiation energy loss fraction"),
        ("v_EM", "float", 1e6, "Aether EM velocity [m/s]"),
        ("B_EM", "float", 1e-5, "HII region B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC2264_m_s2"
    PRIMARY_INPUT = "SFR"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_NGC2264 = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_NGC2264_m_s2": g_NGC2264,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.5, 1.0, 1.5, 2.0, 3.0]
        for val in sweep:
            p = {sweep_param or "M_sf": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #352 — UGC10214TadpoleGalaxyTidalCalculator
# PAPER_768 | Session 181 | thread_06Jun2025.txt
# UGC 10214 Tadpole Galaxy — UQFF Tidal Interaction Dynamics
# g_Tadpole ≈ 3.160e-3 m/s²
# ========================================================================
class UGC10214TadpoleGalaxyTidalCalculator:
    ENTRY = 352
    PAPER = "PAPER_768"
    CPP_MODULE = "UGC10214TadpoleGalaxyTidalCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (10^11 M_sun)"),
        ("r", "float", 1.3e21, "Galaxy radius [m] (~133 kly)"),
        ("z", "float", 0.028, "Redshift (420 Mly)"),
        ("t", "float", 1.578e16, "Interaction time [s] (500 Myr)"),
        ("SFR", "float", 5.0, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.025, "SFR mass fraction"),
        ("T0_tidal", "float", 0.3, "Tidal stripping amplitude"),
        ("tau_tidal", "float", 3.156e16, "Tidal stripping timescale [s] (1 Gyr)"),
        ("v_tidal", "float", 3e5, "Tidal tail velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_Tadpole_m_s2"
    PRIMARY_INPUT = "v_tidal"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        M_tidal = p["T0_tidal"] * (1.0 - math.exp(-p["t"] / p["tau_tidal"]))
        factor_tidal = 1.0 - M_tidal
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_tidal"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_Tadpole = g_grav * factor_Hz * factor_sf * factor_tidal * factor_TRZ + a_EM
        return {
            "g_Tadpole_m_s2": g_Tadpole,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "M_tidal": M_tidal,
            "factor_Hz": factor_Hz,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 2e5, 3e5, 5e5, 1e6]
        for val in sweep:
            p = {sweep_param or "v_tidal": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #353 — NGC4676MiceGalaxiesDualMergerCalculator
# PAPER_769 | Session 181 | thread_06Jun2025.txt
# NGC 4676 Mice Galaxies — UQFF Dual Galaxy Merger Dynamics
# g_Mice ≈ 1.053e-1 m/s²
# ========================================================================
class NGC4676MiceGalaxiesDualMergerCalculator:
    ENTRY = 353
    PAPER = "PAPER_769"
    CPP_MODULE = "NGC4676MiceGalaxiesDualMergerCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_total", "float", 3.978e41, "Total system mass [kg] (2x10^11 M_sun each)"),
        ("r", "float", 3e20, "Interaction radius [m] (~31 kly)"),
        ("z", "float", 0.022, "Redshift (~290 Mly)"),
        ("t", "float", 9.468e15, "Post-encounter time [s] (300 Myr)"),
        ("T0_merge", "float", 0.5, "Merger mass redistribution amplitude"),
        ("tau_merge", "float", 1.262e16, "Merger timescale [s] (400 Myr)"),
        ("v_EM", "float", 1e6, "Aether EM velocity [m/s]"),
        ("B_starburst", "float", 1e-4, "Starburst-enhanced B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_Mice_m_s2"
    PRIMARY_INPUT = "B_starburst"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        H0 = 2.268e-18
        Om, OL = 0.3, 0.7
        z = p["z"]
        Hz = H0 * math.sqrt(Om * (1 + z)**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        M_merge = p["T0_merge"] * (1.0 - math.exp(-p["t"] / p["tau_merge"]))
        factor_merge = 1.0 - M_merge
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M_total"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_starburst"] / m_p) * 11 * 1e-12
        g_Mice = g_grav * factor_Hz * factor_merge * factor_TRZ + a_EM
        return {
            "g_Mice_m_s2": g_Mice,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "M_merge": M_merge,
            "factor_Hz": factor_Hz,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e-5, 5e-5, 1e-4, 5e-4, 1e-3]
        for val in sweep:
            p = {sweep_param or "B_starburst": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)


# ========================================================================
# CP4 #354 — RedSpiderNebulaNG6537UQFFCalculator
# PAPER_770 | Session 181 | thread_06Jun2025.txt
# Red Spider Nebula NGC 6537 — UQFF Bipolar Outflow Planetary Nebula
# g_RedSpider ≈ 2.107e-2 m/s²
# ========================================================================
class RedSpiderNebulaNG6537UQFFCalculator:
    ENTRY = 354
    PAPER = "PAPER_770"
    CPP_MODULE = "RedSpiderNebulaNG6537UQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e30, "WD + ejecta mass [kg] (1 M_sun)"),
        ("r", "float", 1e16, "Nebula radius [m] (~1.06 ly)"),
        ("v_wind", "float", 2e6, "Stellar wind velocity [m/s] (2000 km/s)"),
        ("L_wd", "float", 3.826e30, "WD luminosity [W] (10^4 L_sun)"),
        ("rho_gas", "float", 1e-21, "Nebula gas density [kg/m^3]"),
        ("B_EM", "float", 1e-5, "Wind shock B field [T]"),
        ("z", "float", 0.0013, "Redshift (~4000 ly)"),
    ]
    PRIMARY_OUTPUT = "g_RedSpider_m_s2"
    PRIMARY_INPUT = "v_wind"

    def compute(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11
        m_p = 1.673e-27
        q = 1.602e-19
        c = 3e8
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        F_rad = p["L_wd"] / (4 * math.pi * p["r"]**2)
        P_rad = F_rad / c
        P_rad_term = (P_rad / p["rho_gas"]) * 1e-6 * 1e-3
        a_EM = (q * p["v_wind"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_RedSpider = g_grav * factor_TRZ + P_rad_term + a_EM
        return {
            "g_RedSpider_m_s2": g_RedSpider,
            "g_grav": g_grav,
            "P_rad_term": P_rad_term,
            "a_EM": a_EM,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e5, 1e6, 1.5e6, 2e6, 3e6]
        for val in sweep:
            p = {sweep_param or "v_wind": val}
            r = self.compute(p)
            r["sweep_val"] = val
            results.append(r)
        return results

    def add_mod(self, fn):
        self._mods = getattr(self, "_mods", [])
        self._mods.append(fn)

    def update_from_file(self, filepath):
        import json
        with open(filepath, "r") as f:
            data = json.load(f)
        for k, v in data.items():
            for i, (name, t, d, desc) in enumerate(self.PARAMETERS):
                if name == k:
                    self.PARAMETERS[i] = (name, t, v, desc)

'''

with open("CondensedPhysics4.py", "a", encoding="utf-8") as f:
    f.write(BLOCK)

print("Appended CP4 #345-354 to CondensedPhysics4.py")
print("PAPER_761 (HUDF) → CP4 #345 HubbleUltraDeepFieldUQFFCalculator")
print("PAPER_762 (NGC1792) → CP4 #346 NGC1792StellarForgeUQFFCalculator")
print("PAPER_763 (Sombrero) → CP4 #347 SombreroGalaxyM104UQFFCalculator")
print("PAPER_764 (Saturn) → CP4 #348 Saturn26DUQFFCalculator")
print("PAPER_765 (M16) → CP4 #349 M16EagleNebulaStarsUQFFCalculator")
print("PAPER_766 (Crab) → CP4 #350 CrabNebulaPulsarWindUQFFCalculator")
print("PAPER_767 (NGC2264) → CP4 #351 NGC2264ConeNebulaUQFFCalculator")
print("PAPER_768 (Tadpole) → CP4 #352 UGC10214TadpoleGalaxyTidalCalculator")
print("PAPER_769 (Mice) → CP4 #353 NGC4676MiceGalaxiesDualMergerCalculator")
print("PAPER_770 (RedSpider) → CP4 #354 RedSpiderNebulaNG6537UQFFCalculator")
