#!/usr/bin/env python3
"""
Append TenAstroSystemsMUGECalculator (CP4 #316, PAPER_732)
and EighteenAstroSystemsMUGECalculator (CP4 #317, PAPER_733)
to CondensedPhysics4.py.
Session 178 — v5.35
"""
import math

NEW_CLASSES = '''

# =============================================================================
# CP4 #316 — PAPER_732: 10 Astro-Systems MUGE Compressed UQFF
# Entry #101: (10 astro-systems)_cpp_09May2025
# Systems: NGC 2264, UGC 10214, NGC 4676, Red Spider Nebula, NGC 3372,
#          AG Carinae, M42 (Orion), Tarantula Nebula, NGC 2841, Mystic Mountain
# =============================================================================
class TenAstroSystemsMUGECalculator(object):
    # PAPER_732 | CP4 #316
    # MUGE: g(r,t) = G*M/r^2 * (1+H(z)*t) * (1+M_evo) * (1-E_rad) * (1+f_TRZ) + F_em
    # R(t) = R_grav*cos(omega_grav*t) + R_mag*cos(omega_mag*t) * aether * (1+f_TRZ)

    G          = 6.674e-11
    H_0        = 2.268e-18   # 70 km/s/Mpc in s^-1
    c          = 2.998e8
    hbar       = 1.055e-34
    rho_UA     = 7.09e-36    # aether density (kg/m^3)
    rho_SCm    = 7.09e-37
    q_e        = 1.602e-19
    m_p        = 1.673e-27
    f_TRZ      = 0.1
    em_scale   = 1.0e-12
    M_sun      = 1.989e30    # kg
    aether_ratio = 10.0      # rho_UA / rho_SCm

    SYSTEMS = [
        # (name, M_kg, r_m, z, SFR_sun_per_yr, B_T, v_wind_m_s, tau_erode_s, E_0)
        ("NGC 2264 (Cone Nebula)",          1.989e33, 4.73e16, 0.0006, 0.5,  1.0e-5, 1.0e6,  6.312e13, 0.20),
        ("UGC 10214 (Tadpole Galaxy)",      1.989e41, 1.24e21, 0.028,  1.0,  1.0e-5, 2.0e5,  3.156e16, 0.05),
        ("NGC 4676 (Mice Galaxies)",        3.978e41, 3.0e20,  0.022,  10.0, 1.0e-4, 3.0e5,  3.156e15, 0.10),
        ("Red Spider Nebula",               1.193e30, 1.0e16,  0.0013, 0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15),
        ("NGC 3372 (Carina Nebula)",        1.989e35, 2.0e17,  0.0025, 2.0,  1.0e-5, 1.5e6,  6.312e13, 0.20),
        ("AG Carinae Nebula",               3.978e31, 1.0e16,  0.002,  0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15),
        ("M42 (Orion Nebula)",              3.978e33, 2.0e16,  0.0004, 0.3,  1.0e-5, 1.0e6,  6.312e13, 0.10),
        ("Tarantula Nebula (30 Doradus)",   1.989e35, 3.0e17,  0.0005, 5.0,  1.0e-4, 2.0e6,  3.156e13, 0.25),
        ("NGC 2841 (Spiral Galaxy)",        1.989e41, 5.0e20,  0.0031, 0.5,  1.0e-5, 1.5e5,  3.156e16, 0.05),
        ("Mystic Mountain (Carina pillar)", 1.989e32, 1.0e16,  0.0025, 0.1,  1.0e-5, 1.0e6,  3.156e12, 0.20),
    ]

    def H_z(self, z):
        a = 1.0 + z
        return self.H_0 * math.sqrt(0.3 * a**3 + 0.7)

    def M_evo(self, t, SFR, M_kg):
        t_yr = t / 3.156e7
        M_0_sun = M_kg / self.M_sun
        return SFR * t_yr / M_0_sun if M_0_sun > 0 else 0.0

    def E_rad(self, t, E_0, tau):
        return E_0 * (1.0 - math.exp(-t / tau)) if tau > 0 else 0.0

    def g_muge(self, system_tuple, t=3.156e14):
        name, M, r, z, SFR, B, v_w, tau, E_0 = system_tuple
        g_grav = self.G * M / r**2
        hub    = 1.0 + self.H_z(z) * t
        m_fac  = 1.0 + self.M_evo(t, SFR, M)
        e_fac  = 1.0 - self.E_rad(t, E_0, tau)
        trz    = 1.0 + self.f_TRZ
        F_em   = self.q_e * v_w * B / self.m_p * (1.0 + self.aether_ratio) * self.em_scale
        return g_grav * hub * m_fac * e_fac * trz + F_em

    def R_resonance(self, system_tuple, t=3.156e14):
        name, M, r, z, SFR, B, v_w, tau, E_0 = system_tuple
        g_grav  = self.G * M / r**2
        m_fac   = 1.0 + self.M_evo(t, SFR, M)
        R_grav  = g_grav * m_fac
        R_mag   = self.q_e * v_w * B / self.m_p * self.em_scale
        omega_g = 2.0 * math.pi / tau if tau > 0 else 1.0e-13
        omega_m = omega_g * 100.0
        return (R_grav * math.cos(omega_g * t)
                + R_mag * math.cos(omega_m * t) * self.aether_ratio * (1.0 + self.f_TRZ))

    def compute(self, dataset: dict) -> dict:
        t = dataset.get("t", 3.156e14)
        results = []
        for sys in self.SYSTEMS:
            g = self.g_muge(sys, t)
            R = self.R_resonance(sys, t)
            results.append({"name": sys[0], "g_muge": g, "R_resonance": R})
        return {
            "primary_equation": (
                "g(r,t) = G*M/r^2 * (1+H(z)*t) * (1+M_evo) * (1-E_rad) * (1+f_TRZ) + F_em"
            ),
            "systems": results,
            "paper": "PAPER_732",
            "cp4_index": 316,
        }

    def self_update(self): pass
    def self_expand(self): pass


# =============================================================================
# CP4 #317 — PAPER_733: 18 Astro-Systems MUGE Compressed + 26D UQFF
# Entry #105: (18 astro-systems)_cpp_09May2025
# Systems: 10 from PAPER_732 + NGC 6217, Stephan\'s Quintet, NGC 7049,
#          Carina Nebula (NGC 3324), M74, NGC 1672, NGC 5866, M82,
#          Spirograph Nebula (IC 418)
# =============================================================================
class EighteenAstroSystemsMUGECalculator(object):
    # PAPER_733 | CP4 #317
    # Full MUGE + 26D Ug1-Ug4i quantum state expansion:
    # E_DPM,i = (hbar*c/r_i^2)*Q_i*[SCm]_i  where r_i=r/i, Q_i=i, [SCm]_i=1e-5*i^2
    # Ug1_i = E_DPM,i*(1+H(z)*t)*(1-E_rad)*cos(theta_i)*(1+f_TRZ,i)
    # Ug2_i = E_DPM,i*(1-B/B_crit)*(1+M_sf)*(1+rho_UA/rho_SCm)*cos(omega_j*t)
    # Ug3_i = E_DPM,i*(q*v*B/m_p)*(1-T_lock)*(1+f_TRZ,i)
    # Ug4i_i= (hbar*c/r_THz,i)*(1+f_Um,i)*(1+rho_UA/rho_SCm)

    G          = 6.674e-11
    H_0        = 2.268e-18
    c          = 2.998e8
    hbar       = 1.055e-34
    rho_UA     = 7.09e-36
    rho_SCm    = 7.09e-37
    q_e        = 1.602e-19
    m_p        = 1.673e-27
    M_sun      = 1.989e30
    B_crit     = 4.4e13
    f_TRZ      = 0.1
    em_scale   = 1.0e-12
    aether_ratio = 10.0

    SYSTEMS = [
        # (name, M_kg, r_m, z, SFR, B_T, v_wind, tau_erode, E_0)
        ("NGC 2264 (Cone Nebula)",            1.989e33, 4.73e16, 0.0006, 0.5,  1.0e-5, 1.0e6,  6.312e13, 0.20),
        ("UGC 10214 (Tadpole Galaxy)",        1.989e41, 1.24e21, 0.028,  1.0,  1.0e-5, 2.0e5,  3.156e16, 0.05),
        ("NGC 4676 (Mice Galaxies)",          3.978e41, 3.0e20,  0.022,  10.0, 1.0e-4, 3.0e5,  3.156e15, 0.10),
        ("Red Spider Nebula",                 1.193e30, 1.0e16,  0.0013, 0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15),
        ("NGC 3372 (Carina Nebula)",          1.989e35, 2.0e17,  0.0025, 2.0,  1.0e-5, 1.5e6,  6.312e13, 0.20),
        ("AG Carinae Nebula",                 3.978e31, 1.0e16,  0.002,  0.0,  1.0e-5, 2.0e6,  3.156e11, 0.15),
        ("M42 (Orion Nebula)",                3.978e33, 2.0e16,  0.0004, 0.3,  1.0e-5, 1.0e6,  6.312e13, 0.10),
        ("Tarantula Nebula (30 Doradus)",     1.989e35, 3.0e17,  0.0005, 5.0,  1.0e-4, 2.0e6,  3.156e13, 0.25),
        ("NGC 2841 (Spiral Galaxy)",          1.989e41, 5.0e20,  0.0031, 0.5,  1.0e-5, 1.5e5,  3.156e16, 0.05),
        ("Mystic Mountain (Carina pillar)",   1.989e32, 1.0e16,  0.0025, 0.1,  1.0e-5, 1.0e6,  3.156e12, 0.20),
        ("NGC 6217 (barred spiral)",          1.989e41, 3.0e20,  0.0045, 1.0,  1.0e-5, 1.5e5,  3.156e16, 0.05),
        ("Stephan\'s Quintet",                9.945e41, 1.0e21,  0.022,  10.0, 1.0e-4, 3.0e5,  3.156e15, 0.10),
        ("NGC 7049 (lenticular)",             1.989e41, 5.0e20,  0.0067, 0.2,  1.0e-5, 1.0e5,  3.156e16, 0.05),
        ("Carina Nebula (NGC 3324)",          1.989e35, 2.0e17,  0.0025, 2.0,  1.0e-5, 1.5e6,  6.312e13, 0.20),
        ("M74 (NGC 628, grand-design)",       1.989e41, 5.0e20,  0.0022, 1.0,  1.0e-5, 1.5e5,  3.156e16, 0.05),
        ("NGC 1672 (Seyfert barred)",         1.989e41, 3.0e20,  0.004,  2.0,  1.0e-5, 2.0e5,  3.156e15, 0.08),
        ("NGC 5866 (edge-on lenticular)",     1.989e41, 3.0e20,  0.0029, 0.3,  1.0e-5, 1.0e5,  3.156e16, 0.05),
        ("M82 (Cigar Galaxy, starburst)",     1.989e40, 2.0e20,  0.0008, 10.0, 1.0e-4, 5.0e5,  3.156e14, 0.20),
        ("Spirograph Nebula (IC 418)",        1.193e30, 1.0e16,  0.0007, 0.0,  1.0e-5, 1.5e6,  3.156e11, 0.12),
    ]

    def H_z(self, z):
        a = 1.0 + z
        return self.H_0 * math.sqrt(0.3 * a**3 + 0.7)

    def M_evo(self, t, SFR, M_kg):
        t_yr = t / 3.156e7
        M_0_sun = M_kg / self.M_sun
        return SFR * t_yr / M_0_sun if M_0_sun > 0 else 0.0

    def E_rad(self, t, E_0, tau):
        return E_0 * (1.0 - math.exp(-t / tau)) if tau > 0 else 0.0

    def Ug_26D(self, system_tuple, t):
        name, M, r, z, SFR, B, v_w, tau, E_0 = system_tuple
        hub    = 1.0 + self.H_z(z) * t
        e_fac  = 1.0 - self.E_rad(t, E_0, tau)
        m_fac  = 1.0 + self.M_evo(t, SFR, M)
        B_rat  = B / self.B_crit
        total  = 0.0
        for i in range(1, 27):
            r_i   = r / i
            Q_i   = float(i)
            SCm_i = 1.0e-5 * i * i
            theta_i = i * math.pi / 26.0
            f_TRZ_i = self.f_TRZ * (1.0 + 0.01 * i)
            E_DPM = (self.hbar * self.c / r_i**2) * Q_i * SCm_i
            Ug1 = E_DPM * hub * e_fac * math.cos(theta_i) * (1 + f_TRZ_i)
            Ug2 = E_DPM * (1 - B_rat) * m_fac * (1 + self.aether_ratio) * math.cos(i * 6.28e-13 * t)
            Ug3 = E_DPM * (self.q_e * v_w * B / self.m_p) * (1 + f_TRZ_i)
            r_T = r_i * 1e-3
            f_Um = 0.01 * i
            Ug4i = (self.hbar * self.c / r_T**2) * (1 + f_Um) * (1 + self.aether_ratio)
            total += Ug1 + Ug2 + Ug3 + Ug4i
        return total

    def g_muge(self, system_tuple, t=3.156e14):
        name, M, r, z, SFR, B, v_w, tau, E_0 = system_tuple
        g_grav = self.G * M / r**2
        hub    = 1.0 + self.H_z(z) * t
        m_fac  = 1.0 + self.M_evo(t, SFR, M)
        e_fac  = 1.0 - self.E_rad(t, E_0, tau)
        trz    = 1.0 + self.f_TRZ
        F_em   = self.q_e * v_w * B / self.m_p * (1.0 + self.aether_ratio) * self.em_scale
        ug26   = self.Ug_26D(system_tuple, t)
        return g_grav * hub * m_fac * e_fac * trz + F_em + ug26 * 1.0e-40

    def R_resonance(self, system_tuple, t=3.156e14):
        name, M, r, z, SFR, B, v_w, tau, E_0 = system_tuple
        g_grav  = self.G * M / r**2
        m_fac   = 1.0 + self.M_evo(t, SFR, M)
        R_grav  = g_grav * m_fac
        R_mag   = self.q_e * v_w * B / self.m_p * self.em_scale
        omega_g = 2.0 * math.pi / tau if tau > 0 else 1.0e-13
        omega_m = omega_g * 100.0
        return (R_grav * math.cos(omega_g * t)
                + R_mag * math.cos(omega_m * t) * self.aether_ratio * (1.0 + self.f_TRZ))

    def compute(self, dataset: dict) -> dict:
        t = dataset.get("t", 3.156e14)
        results = []
        for sys in self.SYSTEMS:
            g = self.g_muge(sys, t)
            R = self.R_resonance(sys, t)
            results.append({"name": sys[0], "g_muge": g, "R_resonance": R})
        return {
            "primary_equation": (
                "g(r,t) = G*M/r^2 * (1+H(z)*t) * (1+M_evo) * (1-E_rad) * (1+f_TRZ)"
                " + F_em + Sigma_26D[Ug1+Ug2+Ug3+Ug4i]"
            ),
            "systems": results,
            "paper": "PAPER_733",
            "cp4_index": 317,
        }

    def self_update(self): pass
    def self_expand(self): pass
'''

target = "CondensedPhysics4.py"
with open(target, 'a', encoding='utf-8') as f:
    f.write(NEW_CLASSES)

# Verify
import ast
content = open(target, 'r', encoding='utf-8', errors='replace').read()
classes = [n.name for n in ast.walk(ast.parse(content)) if isinstance(n, ast.ClassDef)]
print(f"CP4 class count: {len(classes)}")
print(f"Last 3 classes: {classes[-3:]}")
