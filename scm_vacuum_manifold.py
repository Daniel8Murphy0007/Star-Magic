# scm_vacuum_manifold.py
# Generated from clean 27FEB2026_A.docx thread + repo alignment
# SCm Vacuum Manifold, Buoyancy, Phonon, Negative-Time, Primordial Split
# Long-form mathematical derivatives included as comments + sympy

import sympy as sp
import numpy as np
from mpmath import li, polylog  # for VDS Li_26

# ==================== VERBATIM CONSTANTS FROM CLEAN THREAD ====================
SSQ = sp.Rational(57, 100)          # [SSq] = 0.57
KAPPA = sp.Rational(5, 10000)       # κ = 5.0 × 10^{-4} day^{-1}
RHO_VAC_SCM = 7.09e-37              # kg/m³
THZ_PHONON = 1.25e12                # 1.25 THz
NEG_TIME_RANGE = sp.symbols('t_n', negative=True)  # t_n < 0

# ==================== LONG-FORM DERIVATIONS ====================

# 1. SCm Vacuum Manifold (primordial "matter" before gravity)
rho_scm = sp.symbols('rho_vac_SCm', positive=True)  # ρ_{vac,SCm} (safe name — no comma)
phi = sp.Function('Φ')(sp.symbols('ω'), sp.symbols('Γ'))  # Gaussian phonon activation

# 2. Negative-Time Modulation
t_n = sp.symbols('t_n', real=True)
cos_pi_tn = sp.cos(sp.pi * t_n)   # flips sign of A_μν and Ubi

# 3. Buoyancy Force F_U_Bi_i (outside-to-inside)
F_0, G, M, r, Omega_g, d_g, wind_mod, U_UA = sp.symbols('F_0 G M r Omega_g d_g wind_mod U_UA', positive=True)
beta_i = sp.symbols(r'\beta_i', positive=True)  # ≈ 0.61
Ug_k = sp.symbols('Ug_k', real=True)  # 4-component

F_U_Bi_i = sp.Integral(
    -F_0 + (G * M / r**2) + rho_scm * U_UA * cos_pi_tn,
    (r, 0, sp.oo)
)  # full long-form integral (thread master)

# 4. 26D Vacuum Density Series (VDS)
n = sp.symbols('n', integer=True, positive=True)
VDS = sp.Sum( (SSQ**n) / n**26 , (n, 1, sp.oo) )   # = Li_26([SSq])

# 5. Phonon Resonance (Holmlid bridge)
omega = sp.symbols('ω', positive=True)
Gamma = sp.symbols('Γ', positive=True)
Phi_gaussian = sp.exp( - (omega - THZ_PHONON)**2 / (2 * Gamma**2) )   # 1.25 THz Gaussian

# 6. Primordial Split: E_net(t, Γ)
E_net = sp.Function('E_net')(t_n, Gamma)   # positive/negative buoyancy branch

# ==================== NUMERICAL HELPERS (for VS Code testing) ====================
def compute_F_U_Bi_i_numerical(M_bh=1.989e30, r=6.96e8, Gamma=1e12):
    """Numerical F_U_Bi_i at a given radius (linearised integrand × displacement).

    F_U_Bi_i integrand = -F_0 + (G*M/r²) + ρ_UA·DPM_stab + Φ·ρ_SCm
    evaluated at the given r, then multiplied by the Gaussian displacement x₂.
    This is the outside-to-inside buoyancy force from the SCm vacuum manifold.
    """
    import math
    G_N = 6.6743e-11
    rho_ua = 7.09e-36
    rho_scm_val = float(RHO_VAC_SCM)
    F_0_val = 1.0e-10
    t_n_default = -100.0   # canonical negative-time

    cos_pi_tn = math.cos(math.pi * t_n_default)
    delta = (2 * math.pi * float(THZ_PHONON)) - (2 * math.pi * float(THZ_PHONON))
    Phi_ph = math.exp(-delta**2 / (2 * max(float(Gamma), 1.0)**2))   # = 1.0 on-resonance

    grav_proj = G_N * float(M_bh) / (float(r)**2) if float(r) > 0 else 0.0
    DPM_stab  = rho_ua * cos_pi_tn
    integrand = -F_0_val + grav_proj * cos_pi_tn + DPM_stab + Phi_ph * rho_scm_val
    x_2       = float(r) * Phi_ph * abs(cos_pi_tn)
    return integrand * x_2

def vds_numerical(terms=1000):
    return float(polylog(26, float(SSQ)))

# ==================== EXPORT FOR LATEX ====================
def export_all_to_latex():
    latex_dict = {
        'rho_scm': sp.latex(rho_scm),
        'F_U_Bi_i': sp.latex(F_U_Bi_i),
        'VDS': sp.latex(VDS),
        'Phi_gaussian': sp.latex(Phi_gaussian),
        'cos_pi_tn': sp.latex(cos_pi_tn),
        'E_net': sp.latex(E_net)
    }
    return latex_dict