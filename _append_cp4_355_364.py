"""
_append_cp4_355_364.py
Session 181 — Batch 3: PAPER_771–780
Appends CP4 classes #355–364 to CondensedPhysics4.py
"""

BLOCK = '''
# ========================================================================
# CP4 #355 — NGC3372EtaCarinaeNebulaUQFFCalculator
# PAPER_771 | Session 181 | thread_06Jun2025.txt
# NGC 3372 Eta Carinae — UQFF LBV Stellar Wind
# g_EtaCar ≈ 5.267e-3 m/s²
# ========================================================================
class NGC3372EtaCarinaeNebulaUQFFCalculator:
    ENTRY = 355
    PAPER = "PAPER_771"
    CPP_MODULE = "NGC3372EtaCarinaeNebulaUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e35, "Nebula mass [kg] (~10^5 M_sun)"),
        ("r", "float", 2e17, "Nebula radius [m] (~21 ly)"),
        ("z", "float", 0.0025, "Redshift"),
        ("t", "float", 3.156e13, "SFR integration time [s] (1 Myr)"),
        ("SFR", "float", 2.0, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.02, "Stellar mass fraction"),
        ("E_rad", "float", 0.15, "Radiation energy loss fraction"),
        ("v_EM", "float", 5e5, "LBV stellar wind velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Nebula B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_EtaCar_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_EtaCar = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_EtaCar_m_s2": g_EtaCar,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 5e5, 1e6, 2e6, 5e6]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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
# CP4 #356 — AGCarinaeNebulaUQFFCalculator
# PAPER_772 | Session 181 | thread_06Jun2025.txt
# AG Carinae — UQFF LBV Eruptive Wind
# g_AGCar ≈ 1.053e-2 m/s²
# ========================================================================
class AGCarinaeNebulaUQFFCalculator:
    ENTRY = 356
    PAPER = "PAPER_772"
    CPP_MODULE = "AGCarinaeNebulaUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 3.978e31, "LBV mass [kg] (~20 M_sun)"),
        ("r", "float", 1e16, "Wind shell radius [m] (~1 ly)"),
        ("z", "float", 0.002, "Redshift"),
        ("t", "float", 9.468e10, "Age [s] (3,000 yr)"),
        ("E_rad", "float", 0.20, "Radiation energy loss fraction"),
        ("v_EM", "float", 1e6, "LBV eruption wind velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Stellar wind B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_AGCar_m_s2"
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
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_AGCar = g_grav * factor_Hz * factor_rad * factor_TRZ + a_EM
        return {
            "g_AGCar_m_s2": g_AGCar,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_rad": factor_rad,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e5, 1e6, 2e6, 5e6, 1e7]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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
# CP4 #357 — M42OrionNebulaUQFFCalculator
# PAPER_773 | Session 181 | thread_06Jun2025.txt
# M42 Orion Nebula — UQFF HII Star Nursery
# g_M42 ≈ 1.053e-3 m/s²
# ========================================================================
class M42OrionNebulaUQFFCalculator:
    ENTRY = 357
    PAPER = "PAPER_773"
    CPP_MODULE = "M42OrionNebulaUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 3.978e33, "Nebula mass [kg] (~2000 M_sun)"),
        ("r", "float", 2e16, "Nebula radius [m] (~2 ly)"),
        ("z", "float", 0.0004, "Redshift"),
        ("t", "float", 9.468e12, "Age [s] (~300,000 yr)"),
        ("SFR", "float", 0.3, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.045, "Stellar mass fraction"),
        ("E_rad", "float", 0.12, "Radiation energy loss fraction"),
        ("v_EM", "float", 1e5, "HII region gas velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Molecular cloud B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_M42_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_M42 = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_M42_m_s2": g_M42,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.1, 0.3, 1.0, 3.0, 10.0]
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
# CP4 #358 — TarantulaNebula30DorUQFFCalculator
# PAPER_774 | Session 181 | thread_06Jun2025.txt
# Tarantula Nebula 30 Doradus — UQFF Extreme Starburst HII
# g_Tarantula ≈ 1.053e-1 m/s²
# ========================================================================
class TarantulaNebula30DorUQFFCalculator:
    ENTRY = 358
    PAPER = "PAPER_774"
    CPP_MODULE = "TarantulaNebula30DorUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e35, "Nebula mass [kg] (~5×10^4 M_sun in ionized region)"),
        ("r", "float", 3e17, "Nebula radius [m] (~31 ly)"),
        ("z", "float", 0.0005, "Redshift (LMC distance)"),
        ("t", "float", 9.468e13, "Age [s] (~3 Myr)"),
        ("SFR", "float", 5.0, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.15, "Stellar mass fraction"),
        ("E_rad", "float", 0.20, "Radiation energy loss fraction"),
        ("v_EM", "float", 1e6, "Starburst wind velocity [m/s]"),
        ("B_EM", "float", 1e-4, "Amplified starburst B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_Tarantula_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_Tar = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_Tarantula_m_s2": g_Tar,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 5e5, 1e6, 5e6, 1e7]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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
# CP4 #359 — NGC2841QuietSpiralUQFFCalculator
# PAPER_775 | Session 181 | thread_06Jun2025.txt
# NGC 2841 Quiet Flocculent Spiral — UQFF Low-SFR Galaxy
# g_NGC2841 ≈ 1.053e-3 m/s²
# ========================================================================
class NGC2841QuietSpiralUQFFCalculator:
    ENTRY = 359
    PAPER = "PAPER_775"
    CPP_MODULE = "NGC2841QuietSpiralUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.02, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 5e20, "Disk radius [m] (~53 kly)"),
        ("z", "float", 0.0031, "Redshift"),
        ("t", "float", 9.468e16, "Evolution time [s] (~3 Gyr)"),
        ("SFR", "float", 0.5, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.015, "Stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC2841_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_NGC2841 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_NGC2841_m_s2": g_NGC2841,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.1, 0.5, 1.0, 2.0, 5.0]
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
# CP4 #360 — MysticMountainCarinaUQFFCalculator
# PAPER_776 | Session 181 | thread_06Jun2025.txt
# Mystic Mountain Carina Pillar — UQFF Protostellar Jet+UV Dual-Forcing
# g_MysticMtn ≈ 1.053e-3 m/s²
# ========================================================================
class MysticMountainCarinaUQFFCalculator:
    ENTRY = 360
    PAPER = "PAPER_776"
    CPP_MODULE = "MysticMountainCarinaUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e32, "Pillar mass [kg] (~100 M_sun)"),
        ("r", "float", 1e16, "Pillar radius [m] (~1 ly)"),
        ("z", "float", 0.0025, "Redshift (Carina distance)"),
        ("t", "float", 1.578e12, "Age [s] (~50,000 yr pillar)"),
        ("SFR", "float", 0.1, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.1, "Stellar mass fraction"),
        ("E_rad", "float", 0.15, "External UV + jet erosion fraction"),
        ("v_EM", "float", 1e5, "Protostellar HH jet velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Molecular cloud B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_MysticMtn_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_MysticMtn = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_MysticMtn_m_s2": g_MysticMtn,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e4, 5e4, 1e5, 5e5, 1e6]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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
# CP4 #361 — NGC6217BarredSpiralUQFFCalculator
# PAPER_777 | Session 181 | thread_06Jun2025.txt
# NGC 6217 Barred Spiral — UQFF SBbc Hubble SM4 First Light
# g_NGC6217 ≈ 1.053e-3 m/s²
# ========================================================================
class NGC6217BarredSpiralUQFFCalculator:
    ENTRY = 361
    PAPER = "PAPER_777"
    CPP_MODULE = "NGC6217BarredSpiralUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.04, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 3e20, "Disk radius [m] (~30 kly)"),
        ("z", "float", 0.0045, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("SFR", "float", 1.0, "Bar-driven star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.045, "Bar-driven stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC6217_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_NGC6217 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_NGC6217_m_s2": g_NGC6217,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.1, 0.5, 1.0, 2.0, 5.0]
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
# CP4 #362 — StephansQuintetGalaxyGroupUQFFCalculator
# PAPER_778 | Session 181 | thread_06Jun2025.txt
# Stephan's Quintet HCG 92 — UQFF Compact Group Shock
# g_SQ ≈ 1.053e-1 m/s²
# ========================================================================
class StephansQuintetGalaxyGroupUQFFCalculator:
    ENTRY = 362
    PAPER = "PAPER_778"
    CPP_MODULE = "StephansQuintetGalaxyGroupUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 9.945e41, "Group mass [kg] (~5×10^11 M_sun, 4 galaxies)"),
        ("r", "float", 1e21, "Group radius [m] (~105 kly)"),
        ("z", "float", 0.022, "Redshift"),
        ("t", "float", 9.468e15, "Starburst age [s] (~300 Myr)"),
        ("SFR", "float", 10.0, "Shock-induced SFR [M_sun/yr]"),
        ("M_sf", "float", 0.05, "SFR mass fraction"),
        ("M_merge", "float", 0.15, "Tidal interaction mass fraction"),
        ("v_EM", "float", 1e6, "Intergalactic shock velocity [m/s]"),
        ("B_EM", "float", 1e-4, "Shock-amplified B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_SQ_m_s2"
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
        factor_sf_merge = 1.0 + p["M_sf"] + p["M_merge"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_SQ = g_grav * factor_Hz * factor_sf_merge * factor_TRZ + a_EM
        return {
            "g_SQ_m_s2": g_SQ,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf_merge": factor_sf_merge,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 5e5, 1e6, 5e6, 1e7]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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
# CP4 #363 — NGC7049LenticularUQFFCalculator
# PAPER_779 | Session 181 | thread_06Jun2025.txt
# NGC 7049 Isolated Lenticular — UQFF Ancient Globular-Rich S0
# g_NGC7049 ≈ 1.053e-3 m/s²
# ========================================================================
class NGC7049LenticularUQFFCalculator:
    ENTRY = 363
    PAPER = "PAPER_779"
    CPP_MODULE = "NGC7049LenticularUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.02, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 5e20, "Galaxy radius [m] (~53 kly)"),
        ("z", "float", 0.0067, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("SFR", "float", 0.2, "Very low SFR [M_sun/yr]"),
        ("M_sf", "float", 0.010, "Stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Quiescent B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC7049_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_NGC7049 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_NGC7049_m_s2": g_NGC7049,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.05, 0.2, 0.5, 1.0, 2.0]
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
# CP4 #364 — CarinaNebulaNGC3324UQFFCalculator
# PAPER_780 | Session 181 | thread_06Jun2025.txt
# Carina Nebula NGC 3324 "Cosmic Cliffs" — UQFF JWST First Light
# g_CosCliffs ≈ 2.107e-3 m/s²
# ========================================================================
class CarinaNebulaNGC3324UQFFCalculator:
    ENTRY = 364
    PAPER = "PAPER_780"
    CPP_MODULE = "CarinaNebulaNGC3324UQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.1, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e35, "Nebula mass [kg] (~10^5 M_sun)"),
        ("r", "float", 2e17, "Cliff radius [m] (~21 ly)"),
        ("z", "float", 0.0025, "Redshift"),
        ("t", "float", 1.578e12, "Age [s] (~50,000 yr)"),
        ("SFR", "float", 2.0, "Active HII SFR [M_sun/yr]"),
        ("M_sf", "float", 0.10, "Stellar mass fraction"),
        ("E_rad", "float", 0.12, "UV photodissociation fraction"),
        ("v_EM", "float", 2e5, "JWST-measured photoionized outflow [m/s]"),
        ("B_EM", "float", 1e-5, "Molecular cloud B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_CosCliffs_m_s2"
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
        factor_sf = 1.0 + p["M_sf"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        g_CosCliffs = g_grav * factor_Hz * factor_sf * factor_rad * factor_TRZ + a_EM
        return {
            "g_CosCliffs_m_s2": g_CosCliffs,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
            "note": "v=2e5 m/s from JWST NIRCam Cosmic Cliffs first-light image",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 2e5, 5e5, 1e6, 2e6]
        for val in sweep:
            p = {sweep_param or "v_EM": val}
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

TARGET = r"c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics4.py"

with open(TARGET, "a", encoding="utf-8") as f:
    f.write(BLOCK)

print(f"Appended CP4 #355-364 to {TARGET}")

# Line count verification
with open(TARGET, "r", encoding="utf-8") as f:
    lines = f.readlines()
print(f"Total lines: {len(lines)}")
