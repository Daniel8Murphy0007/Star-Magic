"""
_append_cp4_365_377.py
Session 181 — Batch 4: PAPER_781–793
Appends CP4 classes #365–377 to CondensedPhysics4.py
(Last 5 standard + 8 Three-UQFF classes)
"""

BLOCK = '''
# ========================================================================
# CP4 #365 — M74PhantomGalaxyUQFFCalculator
# PAPER_781 | Session 181 | thread_06Jun2025.txt
# M74 Phantom Galaxy — UQFF Grand Design Spiral Reference
# g_M74 ≈ 1.053e-3 m/s²
# ========================================================================
class M74PhantomGalaxyUQFFCalculator:
    ENTRY = 365
    PAPER = "PAPER_781"
    CPP_MODULE = "M74PhantomGalaxyUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.04, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 5e20, "Disk radius [m] (~53 kly)"),
        ("z", "float", 0.0022, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("SFR", "float", 1.5, "Star-formation rate [M_sun/yr]"),
        ("M_sf", "float", 0.045, "Stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_M74_m_s2"
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
        g_M74 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_M74_m_s2": g_M74,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.5, 1.0, 1.5, 2.0, 5.0]
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
# CP4 #366 — NGC1672BarredSpiralUQFFCalculator
# PAPER_782 | Session 181 | thread_06Jun2025.txt
# NGC 1672 Barred Spiral — UQFF SB Strong Bar JWST 2023
# g_NGC1672 ≈ 2.107e-3 m/s²
# ========================================================================
class NGC1672BarredSpiralUQFFCalculator:
    ENTRY = 366
    PAPER = "PAPER_782"
    CPP_MODULE = "NGC1672BarredSpiralUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 3e20, "Disk radius [m] (~32 kly)"),
        ("z", "float", 0.004, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("SFR", "float", 3.0, "Bar-driven SFR [M_sun/yr]"),
        ("M_sf", "float", 0.06, "Bar-enhanced mass fraction"),
        ("v_EM", "float", 2e5, "Bar-driven outflow velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC1672_m_s2"
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
        g_NGC1672 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_NGC1672_m_s2": g_NGC1672,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
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


# ========================================================================
# CP4 #367 — NGC5866EdgeOnLenticularUQFFCalculator
# PAPER_783 | Session 181 | thread_06Jun2025.txt
# NGC 5866 Edge-On Lenticular — UQFF Dust Lane S0
# g_NGC5866 ≈ 1.053e-3 m/s²
# ========================================================================
class NGC5866EdgeOnLenticularUQFFCalculator:
    ENTRY = 367
    PAPER = "PAPER_783"
    CPP_MODULE = "NGC5866EdgeOnLenticularUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.02, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 3e20, "Disk radius [m] (~32 kly)"),
        ("z", "float", 0.0029, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("SFR", "float", 0.1, "Very low SFR [M_sun/yr]"),
        ("M_sf", "float", 0.008, "Minimal stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Quiescent B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_NGC5866_m_s2"
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
        g_NGC5866 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_NGC5866_m_s2": g_NGC5866,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_sf": factor_sf,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.01, 0.1, 0.5, 1.0, 2.0]
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
# CP4 #368 — M82CigarStarburstUQFFCalculator
# PAPER_784 | Session 181 | thread_06Jun2025.txt
# M82 Cigar Galaxy — UQFF Starburst Superwind
# g_M82 ≈ 1.053e-1 m/s²
# ========================================================================
class M82CigarStarburstUQFFCalculator:
    ENTRY = 368
    PAPER = "PAPER_784"
    CPP_MODULE = "M82CigarStarburstUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e40, "Galaxy mass [kg] (~10^10 M_sun)"),
        ("r", "float", 2e20, "Disk radius [m] (~21 kly)"),
        ("z", "float", 0.0008, "Redshift"),
        ("t", "float", 3.156e15, "Starburst age [s] (~100 Myr)"),
        ("SFR", "float", 10.0, "Starburst SFR [M_sun/yr]"),
        ("M_sf", "float", 0.15, "Starburst mass fraction"),
        ("v_EM", "float", 1e6, "Superwind velocity [m/s]"),
        ("B_EM", "float", 1e-4, "Starburst-amplified B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_M82_m_s2"
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
        g_M82 = g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM
        return {
            "g_M82_m_s2": g_M82,
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
# CP4 #369 — SpirographNebulaIC418UQFFCalculator
# PAPER_785 | Session 181 | thread_06Jun2025.txt
# Spirograph Nebula IC 418 — UQFF Planetary Nebula Fast Wind
# g_IC418 ≈ 1.580e-2 m/s²
# ========================================================================
class SpirographNebulaIC418UQFFCalculator:
    ENTRY = 369
    PAPER = "PAPER_785"
    CPP_MODULE = "SpirographNebulaIC418UQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.193e30, "Nebula envelope mass [kg] (~0.6 M_sun)"),
        ("r", "float", 1e16, "Nebula radius [m] (~0.1 pc)"),
        ("z", "float", 0.0007, "Redshift (2000 ly distance)"),
        ("t", "float", 9.468e10, "Age [s] (~3000 yr)"),
        ("E_rad", "float", 0.20, "EUV photoionization loss fraction"),
        ("v_EM", "float", 1.5e6, "Central star fast wind [m/s]"),
        ("B_EM", "float", 1e-5, "PN B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_IC418_m_s2"
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
        g_IC418 = g_grav * factor_Hz * factor_rad * factor_TRZ + a_EM
        return {
            "g_IC418_m_s2": g_IC418,
            "g_grav": g_grav,
            "a_EM": a_EM,
            "factor_Hz": factor_Hz,
            "factor_rad": factor_rad,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e5, 1e6, 1.5e6, 2e6, 3e6]
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
# CP4 #370 — NGC4826BlackEyeGalaxyThreeUQFF
# PAPER_786 | Session 181 | thread_06Jun2025.txt
# NGC 4826 Black Eye Galaxy — Three-UQFF Warped Counter-Rotating Disk
# g_primary ≈ 1.053e-3 m/s²
# ========================================================================
class NGC4826BlackEyeGalaxyThreeUQFF:
    ENTRY = 370
    PAPER = "PAPER_786"
    CPP_MODULE = "NGC4826BlackEyeGalaxyThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.04, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 2.83e20, "Effective radius [m] (~30 kly)"),
        ("z", "float", 0.0014, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("M_sf", "float", 0.015, "Stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        R_freq = 1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"]
        return g_c * R_freq

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        g_res = self.compute_resonant(params)
        g_buoy = self.compute_buoyancy(params)
        return {
            "compressed": g_comp,
            "resonant": g_res,
            "buoyancy": g_buoy,
            "g_primary_m_s2": g_comp,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 2e5, 5e5, 1e6]
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
# CP4 #371 — NGC1805LMCClusterThreeUQFF
# PAPER_787 | Session 181 | thread_06Jun2025.txt
# NGC 1805 LMC Cluster — Three-UQFF Dense Young Cluster
# g_primary ≈ 1.053e-3 m/s²
# ========================================================================
class NGC1805LMCClusterThreeUQFF:
    ENTRY = 371
    PAPER = "PAPER_787"
    CPP_MODULE = "NGC1805LMCClusterThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.04, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e34, "Cluster mass [kg] (~10^4 M_sun)"),
        ("r", "float", 9.46e16, "Half-light radius [m] (~3 pc)"),
        ("z", "float", 0.0005, "Redshift (LMC)"),
        ("t", "float", 1.578e16, "Cluster age [s] (~500 Myr)"),
        ("M_sf", "float", 0.05, "Past formation mass fraction"),
        ("v_EM", "float", 1e5, "Cluster velocity dispersion [m/s]"),
        ("B_EM", "float", 1e-5, "LMC B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        R_freq = 1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"]
        return g_c * R_freq

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        g_res = self.compute_resonant(params)
        g_buoy = self.compute_buoyancy(params)
        return {
            "compressed": g_comp,
            "resonant": g_res,
            "buoyancy": g_buoy,
            "g_primary_m_s2": g_comp,
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
# CP4 #372 — NGC6307NGC7027PNPairThreeUQFF
# PAPER_788 | Session 181 | thread_06Jun2025.txt
# NGC 6307 + NGC 7027 PN Pair — Three-UQFF Fast Wind Dual System
# g_primary ≈ 1.580e-2 m/s²
# ========================================================================
class NGC6307NGC7027PNPairThreeUQFF:
    ENTRY = 372
    PAPER = "PAPER_788"
    CPP_MODULE = "NGC6307NGC7027PNPairThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.193e30, "Representative PN mass [kg] (~0.6 M_sun)"),
        ("r", "float", 9.46e15, "Representative radius [m]"),
        ("z", "float", 0.0007, "Redshift"),
        ("t", "float", 9.468e10, "Age [s] (~3000 yr)"),
        ("E_rad", "float", 0.20, "EUV photoionization fraction"),
        ("v_EM", "float", 1.5e6, "Central star fast wind [m/s]"),
        ("B_EM", "float", 1e-5, "PN B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_rad * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        R_freq = 1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"]
        return g_c * R_freq

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        g_res = self.compute_resonant(params)
        g_buoy = self.compute_buoyancy(params)
        return {
            "compressed": g_comp,
            "resonant": g_res,
            "buoyancy": g_buoy,
            "g_primary_m_s2": g_comp,
            "NGC6307_g": g_comp,
            "NGC7027_g": g_comp,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e5, 1e6, 1.5e6, 2e6, 3e6]
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
# CP4 #373 — CassiniRingGapsThreeUQFFCalculator
# PAPER_789 | Session 181 | thread_06Jun2025.txt
# Cassini Ring Gaps — Three-UQFF Saturn Ring Resonance
# g_Cassini_Division = 2.635 m/s²
# ========================================================================
class CassiniRingGapsThreeUQFFCalculator:
    ENTRY = 373
    PAPER = "PAPER_789"
    CPP_MODULE = "CassiniRingGapsThreeUQFFCalculator"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.0, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M_Saturn", "float", 5.683e26, "Saturn mass [kg]"),
        ("r_Encke", "float", 1.335e8, "Encke Gap radius [m]"),
        ("r_Cassini", "float", 1.200e8, "Cassini Division mid radius [m]"),
        ("r_Maxwell", "float", 8.748e7, "Maxwell Gap radius [m]"),
        ("B_ring", "float", 1e-7, "Saturn ring plane B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_Cassini_Division_m_s2"
    PRIMARY_INPUT = "M_Saturn"

    def _g_gap(self, M, r, B):
        import math
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        g_grav = G * M / r**2
        v_orb = math.sqrt(G * M / r)
        a_EM = (q * v_orb * B / m_p) * 11 * 1e-12
        return g_grav + a_EM

    def compute_compressed(self, params=None):
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        return self._g_gap(p["M_Saturn"], p["r_Cassini"], p["B_ring"])

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        R_freq = 1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"]
        return g_c * R_freq

    def compute_buoyancy(self, params=None):
        return self.compute_compressed(params)

    def compute(self, params=None):
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_Encke = self._g_gap(p["M_Saturn"], p["r_Encke"], p["B_ring"])
        g_Cassini = self._g_gap(p["M_Saturn"], p["r_Cassini"], p["B_ring"])
        g_Maxwell = self._g_gap(p["M_Saturn"], p["r_Maxwell"], p["B_ring"])
        R_freq = 1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"]
        return {
            "compressed": g_Cassini,
            "resonant": g_Cassini * R_freq,
            "buoyancy": g_Cassini,
            "g_Cassini_Division_m_s2": g_Cassini,
            "g_Encke_m_s2": g_Encke,
            "g_Maxwell_m_s2": g_Maxwell,
            "g_primary_m_s2": g_Cassini,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [8e7, 1e8, 1.2e8, 1.335e8, 1.5e8]
        for val in sweep:
            p = {sweep_param or "r_Cassini": val}
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
# CP4 #374 — ESO391_12LenticularThreeUQFF
# PAPER_790 | Session 181 | thread_06Jun2025.txt
# ESO 391-12 Lenticular with Dust Ring — Three-UQFF Extended S0
# g_primary ≈ 1.053e-3 m/s²
# ========================================================================
class ESO391_12LenticularThreeUQFF:
    ENTRY = 374
    PAPER = "PAPER_790"
    CPP_MODULE = "ESO391_12LenticularThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.02, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 4.73e20, "Effective radius [m] (~50 kly)"),
        ("z", "float", 0.0067, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("M_sf", "float", 0.008, "Minimal stellar mass fraction"),
        ("v_EM", "float", 1e5, "Disk rotation [m/s]"),
        ("B_EM", "float", 1e-5, "Quiescent B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        return g_c * (1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"])

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        return {
            "compressed": g_comp,
            "resonant": self.compute_resonant(params),
            "buoyancy": self.compute_buoyancy(params),
            "g_primary_m_s2": g_comp,
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [1e5, 2e5, 5e5, 1e6]
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
# CP4 #375 — M57RingNebulaThreeUQFF
# PAPER_791 | Session 181 | thread_06Jun2025.txt
# M57 Ring Nebula — Three-UQFF Planetary Nebula Archetype
# g_primary ≈ 1.580e-2 m/s²
# ========================================================================
class M57RingNebulaThreeUQFF:
    ENTRY = 375
    PAPER = "PAPER_791"
    CPP_MODULE = "M57RingNebulaThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.193e30, "Nebula shell mass [kg] (~0.6 M_sun)"),
        ("r", "float", 1.89e15, "Inner ring radius [m] (~0.2 pc)"),
        ("z", "float", 0.0008, "Redshift (2300 ly)"),
        ("t", "float", 1.262e11, "Age [s] (~4000 yr)"),
        ("E_rad", "float", 0.18, "EUV photoionization fraction"),
        ("v_EM", "float", 1.5e6, "Central star fast wind [m/s]"),
        ("B_EM", "float", 1e-5, "PN B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_rad = 1.0 - p["E_rad"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_rad * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        return g_c * (1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"])

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        return {
            "compressed": g_comp,
            "resonant": self.compute_resonant(params),
            "buoyancy": self.compute_buoyancy(params),
            "g_primary_m_s2": g_comp,
            "note": "JWST 2023 revealed 3D barrel structure; inner ring r=1.89e15 m used",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e5, 1e6, 1.5e6, 2e6, 3e6]
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
# CP4 #376 — LargeMagellanicCloudThreeUQFF
# PAPER_792 | Session 181 | thread_06Jun2025.txt
# Large Magellanic Cloud — Three-UQFF Irregular Satellite Galaxy
# g_primary ≈ 1.053e-3 m/s²
# ========================================================================
class LargeMagellanicCloudThreeUQFF:
    ENTRY = 376
    PAPER = "PAPER_792"
    CPP_MODULE = "LargeMagellanicCloudThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.04, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e40, "LMC mass [kg] (~10^10 M_sun)"),
        ("r", "float", 6.62e19, "LMC radius [m] (~7 kly)"),
        ("z", "float", 0.0005, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("M_sf", "float", 0.05, "LMC-wide SFR mass fraction"),
        ("v_EM", "float", 1e5, "LMC rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "LMC B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_TRZ = 1.0 + self.UQFF_CONSTANTS["F_TRZ"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_sf * factor_TRZ + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        return g_c * (1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"])

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        return {
            "compressed": g_comp,
            "resonant": self.compute_resonant(params),
            "buoyancy": self.compute_buoyancy(params),
            "g_primary_m_s2": g_comp,
            "note": "Compare PAPER_774 Tarantula sub-region: g=1.053e-1 (local starburst v=1e6, B=1e-4)",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [5e4, 1e5, 2e5, 5e5, 1e6]
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
# CP4 #377 — ESO510G13WarpedSpiralThreeUQFF
# PAPER_793 | Session 181 | thread_06Jun2025.txt
# ESO 510-G13 Warped Spiral — Three-UQFF Geometry Invariance
# g_primary ≈ 1.053e-3 m/s²
# UQFF Geometry Invariance Theorem: g independent of disk warp amplitude
# ========================================================================
class ESO510G13WarpedSpiralThreeUQFF:
    ENTRY = 377
    PAPER = "PAPER_793"
    CPP_MODULE = "ESO510G13WarpedSpiralThreeUQFF"
    UQFF_CONSTANTS = {
        "RHO_UA": 7.09e-36, "RHO_SCM": 7.09e-37,
        "F_TRZ": 0.05, "KAPPA": 0.0005, "SSQ": 0.57,
        "MU_J": 3.38e23, "GAMMA": 5.0e-5 / 86400.0,
    }
    PARAMETERS = [
        ("M", "float", 1.989e41, "Galaxy mass [kg] (~10^11 M_sun)"),
        ("r", "float", 3.78e20, "Disk radius [m] (~40 kly)"),
        ("z", "float", 0.010, "Redshift"),
        ("t", "float", 1.578e17, "Evolution time [s] (~5 Gyr)"),
        ("M_sf", "float", 0.03, "Stellar mass fraction"),
        ("f_warp", "float", 0.05, "Warp geometry correction"),
        ("v_EM", "float", 1e5, "Disk rotation velocity [m/s]"),
        ("B_EM", "float", 1e-5, "Galactic B field [T]"),
    ]
    PRIMARY_OUTPUT = "g_primary_m_s2"
    PRIMARY_INPUT = "M"

    def compute_compressed(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        G = 6.6743e-11; m_p = 1.673e-27; q = 1.602e-19
        H0 = 2.268e-18; Om, OL = 0.3, 0.7
        Hz = H0 * math.sqrt(Om * (1 + p["z"])**3 + OL)
        factor_Hz = 1.0 + Hz * p["t"]
        factor_sf = 1.0 + p["M_sf"]
        factor_warp = 1.0 + p["f_warp"]
        g_grav = G * p["M"] / p["r"]**2
        a_EM = (q * p["v_EM"] * p["B_EM"] / m_p) * 11 * 1e-12
        return g_grav * factor_Hz * factor_sf * factor_warp + a_EM

    def compute_resonant(self, params=None):
        g_c = self.compute_compressed(params)
        return g_c * (1.0 + self.UQFF_CONSTANTS["KAPPA"] * self.UQFF_CONSTANTS["SSQ"])

    def compute_buoyancy(self, params=None):
        import math
        p = {k: d for k, _, d, _ in self.PARAMETERS}
        if params:
            p.update(params)
        g_c = self.compute_compressed(params)
        V = (4.0/3.0) * math.pi * p["r"]**3
        a_Ubi = self.UQFF_CONSTANTS["RHO_UA"] * V * g_c / 1.673e-27
        return g_c + a_Ubi

    def compute(self, params=None):
        g_comp = self.compute_compressed(params)
        return {
            "compressed": g_comp,
            "resonant": self.compute_resonant(params),
            "buoyancy": self.compute_buoyancy(params),
            "g_primary_m_s2": g_comp,
            "note": "UQFF Geometry Invariance: 90-degree S-warp does not alter EM ground state",
        }

    def simulate(self, sweep=None, sweep_param=None):
        results = []
        sweep = sweep or [0.0, 0.05, 0.1, 0.2, 0.5]
        for val in sweep:
            p = {sweep_param or "f_warp": val}
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

print(f"Appended CP4 #365-377 to {TARGET}")

with open(TARGET, "r", encoding="utf-8") as f:
    lines = f.readlines()
print(f"Total lines: {len(lines)}")
