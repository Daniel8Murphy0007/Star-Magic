"""Session 177 -- Append CP4 #315 NGC1316MergerEvolutionCalculator to CondensedPhysics4.py"""
import os, re

CP4 = r"C:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic\CondensedPhysics4.py"

cls = "NGC1316MergerEvolutionCalculator"

UQFF_CONSTANTS = """\
    # UQFF universal constants
    G       = 6.6743e-11
    c       = 3.0e8
    hbar    = 1.0546e-34
    mu_0    = 1.2566e-6
    k_B     = 1.3806e-23
    M_sun   = 1.989e30
    kpc     = 3.086e19
    Mpc     = 3.086e22
    rho_UA  = 7.09e-36
    rho_SCm = 7.09e-37
    f_TRZ   = 0.1
    kappa   = 5.0e-4
    SSq     = 0.57
    mu_J    = 3.38e23
    Lambda  = 1.1e-52
    H_0     = 2.269e-18
    t_H     = 4.35e17"""

ENTRY = f"""

class {cls}(object):
    # PAPER_731: NGC 1316 (Fornax A) merger-driven elliptical galaxy evolution
    # MUGE: g_NGC1316 = G*M(t)/r^2 * (1+H(z)) * (1-B/Bcrit) * (1+F_env)
    #       + (U_g1+U_g2+U_g3prime+U_g4) + U_i + Lambda*c^2/3
    #       + psi_dust quantum term + rho_dust*V*g_local + (M_vis+M_DM)*dm_term
    # Source: grok_share_ba508f76c8e.txt entry #64 | Session 177
{UQFF_CONSTANTS}

    # NGC 1316 system parameters
    M_visible    = 3.5e11   # M_sun  -- stellar visible mass
    M_DM         = 1.5e11   # M_sun  -- dark matter halo
    M_spiral     = 1.0e10   # M_sun  -- merger progenitor
    M_cluster    = 1.0e6    # M_sun  -- globular cluster mass
    r_0_kpc      = 46.0e3   # kpc    -- galaxy radius
    d_spiral_kpc = 50.0e3   # kpc    -- progenitor distance
    sigma_kpc    = 2.0e3    # kpc    -- dust lane sigma
    tau_merge    = 3.156e16 # s      -- merger decay timescale (~1 Gyr)
    z_redshift   = 0.005    # NGC 1316 redshift
    B_AGN        = 1.0e-4   # T      -- AGN magnetic field
    B_crit       = 1.0e11   # T      -- critical field
    omega_spin   = 1.0e-3   # rad/s  -- BH spin
    H_aether     = 1.0e-5   # A/m    -- aether field
    rho_dust     = 1.0e-21  # kg/m^3 -- dust density
    Vol_galaxy   = 1.0e51   # m^3    -- galactic volume
    k_cluster    = 1.0e-12  # N/M_sun
    lambda_I     = 1.0
    F_RZ         = 0.01
    omega_i      = 1.0e-8   # rad/s
    t_n          = 0.0
    k_4          = 1.0

    def __init__(self, dataset=None):
        self.dataset = dataset or {{}}
        self.version = "Session177"

    def compute(self, dataset=None):
        d = dataset or self.dataset
        return {{
            "paper":     "PAPER_731",
            "cp4_entry": 315,
            "class":     "{cls}",
            "domain":    "NGC 1316 Fornax A merger evolution MUGE AGN dust lanes dark matter",
            "equations": self._primary_equations(d),
        }}

    def _primary_equations(self, d):
        import math
        t     = d.get("t", 0.0)
        r_kpc = d.get("r_kpc", 46.0e3)
        r     = r_kpc * self.kpc
        M_total = (self.M_visible + self.M_DM + self.M_spiral * math.exp(-t / self.tau_merge)) * self.M_sun
        r_eff   = r + 1.0e3 * t
        H_z     = self.H_0 * math.sqrt(0.3 * (1.0 + self.z_redshift)**3 + 0.7)
        F_tidal   = self.G * self.M_spiral * self.M_sun / (self.d_spiral_kpc * self.kpc)**2
        F_cluster = self.k_cluster * self.M_cluster * self.M_sun
        U_g4 = self.k_4 * 1.0e46 * math.exp(-self.kappa * t)
        U_i  = self.lambda_I * (self.rho_SCm / self.rho_UA) * self.omega_i * (1.0 + self.F_RZ)
        g_core = self.G * M_total / (r**2)
        return {{
            "M_total_kg":  M_total,
            "r_eff_m":     r_eff,
            "H_z_inv_s":   H_z,
            "F_tidal":     F_tidal,
            "F_cluster":   F_cluster,
            "U_g4":        U_g4,
            "U_i":         U_i,
            "g_core_m_s2": g_core,
            "Lambda_cosm": self.Lambda * self.c**2 / 3.0,
            "UQFF_factor": 1.0 + self.rho_SCm / self.rho_UA,
        }}

    def self_update(self): pass
    def self_expand(self): pass

"""

def main():
    if not os.path.exists(CP4):
        print(f"ERROR: {CP4} not found")
        return
    with open(CP4, "r", encoding="latin-1") as f:
        content = f.read()
    if f"class {cls}" in content:
        print(f"SKIP: {cls} already in CP4")
        return
    with open(CP4, "a", encoding="utf-8") as f:
        f.write(ENTRY)
    print(f"Appended {cls} (CP4 #315)")
    with open(CP4, "r", encoding="latin-1") as f:
        c2 = f.read()
    count = len(re.findall(r"^class\s+\w+", c2, re.MULTILINE))
    print(f"CP4 total classes now: {count}")

if __name__ == "__main__":
    main()
