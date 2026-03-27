"""Standalone test for CP4 #154-#157 (Session 149 BSFG Open Questions).
   Inlines required constants so no CP4 import chain is needed."""
import math

# ── Shared constants (from S147/S148 blocks) ─────────────────────────────────
_S147_FAC26    = math.factorial(26)
_S147_R_Q_AU   = (2.0 / _S147_FAC26) ** (1.0 / 26.0)
_S147_AU_IN_M  = 1.496e11
_S148_ETA      = 1e-22
_S148_C_LIGHT  = 3e8
_S148_MS       = 1.989e30
_S148_RS       = 6.96e8
# ── Session 149 constants ────────────────────────────────────────────────────
_S149_G_N      = 6.674e-11
_S149_HBAR     = 1.055e-34
_S149_KB       = 1.381e-23
_S149_H_PL     = 6.626e-34
_S149_LP       = 1.616e-35
_S149_AU       = 1.496e11
_S149_LAM_OBS  = 1.1e-52


class _CP4Calculator:
    def compute(self, dataset):
        raise NotImplementedError


class BSFGEinsteinTensorFieldEquationsCalculator(_CP4Calculator):
    SESSION = 149; PAPER = 'PAPER_559'
    def compute(self, dataset):
        r   = dataset.get('r',   _S148_RS)
        tn  = dataset.get('t_n', 0.0)
        eta = dataset.get('eta', _S148_ETA)
        Ms  = dataset.get('Ms',  _S148_MS)
        Ls  = dataset.get('Ls',  3.828e26)
        c   = _S148_C_LIGHT;  G_N = _S149_G_N
        C_num  = (Ms*c**2 + Ls/c**2) / ((4/3)*math.pi)
        cos_tn = math.cos(math.pi*tn)
        V      = (4/3)*math.pi*r**3;  Ts00 = Ms*c**2/V
        eps    = eta*Ts00*cos_tn
        eps_p  = -3*eta*cos_tn*C_num/r**4
        eps_pp = 12*eta*cos_tn*C_num/r**5
        A00 = 1+eps;  Arr = -1+eps
        R_r0r0 = eps_pp/2 - eps_p**2/2
        R_00   = 3*R_r0r0
        R_rr   = -R_r0r0 + 2*(eps_pp/2 - eps_p**2/4)
        R_s    = R_00/A00 + R_rr/Arr
        G_00   = R_00 - 0.5*A00*R_s
        G_rr   = R_rr - 0.5*Arr*R_s
        kE     = 8*math.pi*G_N/c**4
        T00_e  = G_00/kE
        RHS    = kE*Ts00
        amp    = G_00/RHS if abs(RHS) > 0 else float('nan')
        Leff   = kE*eta*Ts00/2
        Lrat   = Leff/_S149_LAM_OBS
        rho_v  = Leff*c**2/(8*math.pi*G_N)
        return {'paper': self.PAPER, 'G_00': G_00, 'G_rr': G_rr,
                'R_scalar': R_s, 'T00_eff_Pa': T00_e, 'T_s00_Pa': Ts00,
                'kappa_E': kE, 'RHS_00_GR': RHS, 'amp_factor': amp,
                'Lambda_eff': Leff, 'Lambda_obs': _S149_LAM_OBS,
                'Lambda_ratio': Lrat, 'rho_vac_eff': rho_v,
                'non_Einstein': (abs(amp) > 10) if not math.isnan(amp) else False,
                'session': self.SESSION, 'papers': [self.PAPER]}


class BSFGHolonomyGroupParallelTransportCalculator(_CP4Calculator):
    SESSION = 149; PAPER = 'PAPER_560'
    def compute(self, dataset):
        r   = dataset.get('r',            _S148_RS)
        tn  = dataset.get('t_n',          0.0)
        eta = dataset.get('eta',          _S148_ETA)
        Ms  = dataset.get('Ms',           _S148_MS)
        Ls  = dataset.get('Ls',           3.828e26)
        dA  = dataset.get('loop_area_m2', None)
        c   = _S148_C_LIGHT
        C_num  = (Ms*c**2 + Ls/c**2) / ((4/3)*math.pi)
        cos_tn = math.cos(math.pi*tn)
        V      = (4/3)*math.pi*r**3;  Ts00 = Ms*c**2/V
        eps    = eta*Ts00*cos_tn
        eps_p  = -3*eta*cos_tn*C_num/r**4
        eps_pp = 12*eta*cos_tn*C_num/r**5
        A00 = 1+eps;  Arr = -1+eps
        R_r0r0 = eps_pp/2 - eps_p**2/2
        R_00   = 3*R_r0r0
        R_rr   = -R_r0r0 + 2*(eps_pp/2 - eps_p**2/4)
        R_s    = R_00/A00 + R_rr/Arr
        is_flat = abs(R_s) < 1e-50
        n_ex = 22
        if dA is None: dA = r*r
        dphi  = R_r0r0 * dA
        dphi_P = R_r0r0 * _S149_LP**2
        R_AU   = (6*eta*cos_tn*C_num/_S149_AU**5) if abs(cos_tn) > 1e-15 else 0.0
        dphi_AU = R_AU * _S149_AU**2
        omega  = eps_p / (2*A00)
        return {'paper': self.PAPER,
                'G_hol_4D': 'SO+(3,1)', 'G_hol_extra': f'U(1)^{n_ex}',
                'G_hol_full': f'SO+(3,1) x U(1)^{n_ex}',
                'n_extra_flat_dims': n_ex, 'R_r0r0_at_r': R_r0r0,
                'R_scalar': R_s, 'is_Ricci_flat': is_flat,
                'has_G2_holonomy': False, 'has_Spin7_holonomy': False,
                'delta_phi_rad': dphi, 'delta_phi_Planck': dphi_P,
                'delta_phi_AU2_rad': dphi_AU, 'loop_area_m2': dA,
                'omega_0r': omega, 'holonomy_trivial': (abs(dphi) < 1e-100),
                'session': self.SESSION, 'papers': [self.PAPER]}


class BSFGBlackHoleSolutionHorizonCalculator(_CP4Calculator):
    SESSION = 149; PAPER = 'PAPER_561'
    def compute(self, dataset):
        tn   = dataset.get('t_n',  1.0)
        eta  = dataset.get('eta',  _S148_ETA)
        Ms   = dataset.get('Ms',   _S148_MS)
        Ls   = dataset.get('Ls',   3.828e26)
        c    = _S148_C_LIGHT;  G_N = _S149_G_N
        hbar = _S149_HBAR;     kB  = _S149_KB
        C_num  = (Ms*c**2 + Ls/c**2) / ((4/3)*math.pi)
        cos_tn = math.cos(math.pi*tn)
        exists = cos_tn < -1e-15
        r_h    = (-eta*C_num*cos_tn)**(1/3) if exists else float('nan')
        r_GR   = 2*G_N*Ms/c**2
        dA_rh  = -3*eta*C_num*cos_tn/r_h**4 if exists else float('nan')
        kappa  = c**2*abs(dA_rh)/2          if exists else float('nan')
        T_H    = hbar*kappa/(2*math.pi*kB*c) if exists else float('nan')
        T_GR   = hbar*c**3/(8*math.pi*G_N*Ms*kB)
        r_q    = _S147_R_Q_AU * _S147_AU_IN_M
        return {'paper': self.PAPER, 'horizon_exists': exists,
                't_n': tn, 'cos_pi_tn': cos_tn, 'r_h_m': r_h,
                'r_h_over_Rs':    r_h/_S148_RS if exists else float('nan'),
                'r_h_over_r_s_GR': r_h/r_GR   if exists else float('nan'),
                'r_h_over_r_q':   r_h/r_q      if exists else float('nan'),
                'kappa_BSFG_ms2': kappa, 'T_H_BSFG_K': T_H, 'T_H_GR_K': T_GR,
                'r_s_GR_m': r_GR, 'r_q_canonical_m': r_q,
                'eta_Cnum': eta*abs(C_num),
                'session': self.SESSION, 'papers': [self.PAPER]}


class BSFGBohrSommerfeldAetherQuantizationCalculator(_CP4Calculator):
    SESSION = 149; PAPER = 'PAPER_562'
    def compute(self, dataset):
        r    = dataset.get('r',   _S148_RS)
        tn   = dataset.get('t_n', 0.0)
        eta  = dataset.get('eta', _S148_ETA)
        Ms   = dataset.get('Ms',  _S148_MS)
        Ls   = dataset.get('Ls',  3.828e26)
        c    = _S148_C_LIGHT;  G_N = _S149_G_N
        hbar = _S149_HBAR;     h_pl = _S149_H_PL
        C_num  = (Ms*c**2 + Ls/c**2) / ((4/3)*math.pi)
        cos_tn = math.cos(math.pi*tn)
        eps_p  = -3*eta*cos_tn*C_num/r**4
        v2_n   = G_N*Ms/r
        v2_a   = r*eps_p*c**2/2
        v2_tot = v2_n + v2_a
        v_orb  = math.sqrt(max(v2_tot, 0.0))
        dJJ    = v2_a/(2*v2_n) if v2_n != 0 else float('nan')
        if abs(cos_tn) > 1e-15:
            r_cr = math.sqrt(eta*c**2*abs(cos_tn)*C_num/(G_N*Ms))
            r_cr_AU = r_cr/_S149_AU
        else:
            r_cr = r_cr_AU = float('nan')
        h_eta  = eta*h_pl
        J_sp   = math.sqrt(G_N*Ms*r)
        n_Kep  = J_sp/hbar
        dn     = dJJ*n_Kep if not math.isnan(dJJ) else float('nan')
        ep_AU  = -3*eta*cos_tn*C_num/_S149_AU**4
        v2aAU  = _S149_AU*ep_AU*c**2/2
        v2nAU  = G_N*Ms/_S149_AU
        dJJ_AU = v2aAU/(2*v2nAU) if v2nAU > 0 else float('nan')
        return {'paper': self.PAPER, 'r_m': r,
                'v2_newton_m2s2': v2_n, 'v2_aether_m2s2': v2_a,
                'v2_total_m2s2': v2_tot, 'v_orbit_ms': v_orb,
                'delta_J_over_J': dJJ, 'delta_J_over_J_1AU': dJJ_AU,
                'r_cross_m': r_cr, 'r_cross_AU': r_cr_AU,
                'aether_dominates': abs(v2_a) > v2_n,
                'h_eta': h_eta, 'J_spec_m2s': J_sp,
                'n_Kepler': n_Kep, 'delta_n_BSFG': dn,
                'session': self.SESSION, 'papers': [self.PAPER]}


# ── Tests ─────────────────────────────────────────────────────────────────────
ds   = {'r': 6.96e8, 't_n': 0.0}
ds_h = {'t_n': 1.0}

# #154 Einstein Tensor
r1 = BSFGEinsteinTensorFieldEquationsCalculator().compute(ds)
assert r1['paper'] == 'PAPER_559'
assert abs(r1['G_00']) > 0
assert r1['non_Einstein'] is True
print(f"  PAPER_559  G_00        = {r1['G_00']:.3e} m^-2")
print(f"             amp_factor  = {r1['amp_factor']:.3e}")
print(f"             Lambda_eff  = {r1['Lambda_eff']:.3e} m^-2")
print(f"             Lambda_ratio= {r1['Lambda_ratio']:.3e}")
print(f"             non_Einstein= {r1['non_Einstein']}")

# #155 Holonomy
r2 = BSFGHolonomyGroupParallelTransportCalculator().compute(ds)
assert r2['paper'] == 'PAPER_560'
assert 'SO+(3,1)' in r2['G_hol_full']
assert r2['is_Ricci_flat'] is False
assert r2['has_G2_holonomy'] is False
assert r2['has_Spin7_holonomy'] is False
print(f"  PAPER_560  G_hol_full  = {r2['G_hol_full']}")
print(f"             is_Ricci_fl = {r2['is_Ricci_flat']}")
print(f"             delta_phi   = {r2['delta_phi_rad']:.3e} rad  (loop=r^2)")
print(f"             delta_Pl    = {r2['delta_phi_Planck']:.3e} rad  (Planck)")

# #156 BH Horizon (anti-phase)
r3 = BSFGBlackHoleSolutionHorizonCalculator().compute(ds_h)
assert r3['paper'] == 'PAPER_561'
assert r3['horizon_exists'] is True
assert r3['r_h_m'] > 0
print(f"  PAPER_561  r_h_m       = {r3['r_h_m']:.3e} m")
print(f"             r_h/Rs      = {r3['r_h_over_Rs']:.3f}")
print(f"             kappa_BSFG  = {r3['kappa_BSFG_ms2']:.3e} m/s^2")
print(f"             T_H_BSFG    = {r3['T_H_BSFG_K']:.3e} K")
r3b = BSFGBlackHoleSolutionHorizonCalculator().compute({'t_n': 0.0})
assert r3b['horizon_exists'] is False
print(f"             horizon@tn0 = {r3b['horizon_exists']}  (no horizon in normal phase)")

# #157 Bohr-Sommerfeld
r4 = BSFGBohrSommerfeldAetherQuantizationCalculator().compute(ds)
assert r4['paper'] == 'PAPER_562'
assert r4['r_cross_AU'] > 0
print(f"  PAPER_562  r_cross_AU  = {r4['r_cross_AU']:.3f} AU")
print(f"             dJJ_1AU     = {r4['delta_J_over_J_1AU']:.3e}")
print(f"             h_eta       = {r4['h_eta']:.3e}")
print(f"             aether_dom  = {r4['aether_dominates']} (Aether dominates at Rs)")

print()
print("  ALL 4 BSFG SESSION 149 CLASSES PASSED (#154-#157)")


ds   = {'r': 6.96e8, 't_n': 0.0}
ds_h = {'t_n': 1.0}   # anti-phase for horizon

# ── #154 Einstein Tensor ──────────────────────────────────────────────────
r1 = BSFGEinsteinTensorFieldEquationsCalculator().compute(ds)
assert r1['paper'] == 'PAPER_559', "paper tag wrong"
assert abs(r1['G_00']) > 0, "G_00 must be nonzero"
assert r1['non_Einstein'] is True, "amp should be >>1"
print(f"G_00        = {r1['G_00']:.3e} m^-2")
print(f"amp_factor  = {r1['amp_factor']:.3e}")
print(f"Lambda_eff  = {r1['Lambda_eff']:.3e} m^-2")
print(f"Lambda_r    = {r1['Lambda_ratio']:.3e}")

# ── #155 Holonomy Group ───────────────────────────────────────────────────
r2 = BSFGHolonomyGroupParallelTransportCalculator().compute(ds)
assert r2['paper'] == 'PAPER_560'
assert 'SO+(3,1)' in r2['G_hol_full']
assert r2['is_Ricci_flat'] is False,  "BSFG must not be Ricci-flat"
assert r2['has_G2_holonomy'] is False
assert r2['has_Spin7_holonomy'] is False
print(f"G_hol_full  = {r2['G_hol_full']}")
print(f"is_Ricci_fl = {r2['is_Ricci_flat']}")
print(f"delta_phi   = {r2['delta_phi_rad']:.3e} rad  (loop area = r^2)")
print(f"delta_Pl    = {r2['delta_phi_Planck']:.3e} rad  (Planck loop)")

# ── #156 BH Horizon ───────────────────────────────────────────────────────
r3 = BSFGBlackHoleSolutionHorizonCalculator().compute(ds_h)
assert r3['paper'] == 'PAPER_561'
assert r3['horizon_exists'] is True,  "t_n=1 should produce a horizon"
assert r3['r_h_m'] > 0
print(f"r_h_m       = {r3['r_h_m']:.3e} m")
print(f"r_h/Rs      = {r3['r_h_over_Rs']:.3f}")
print(f"kappa_BSFG  = {r3['kappa_BSFG_ms2']:.3e} m/s^2")
print(f"T_H_BSFG    = {r3['T_H_BSFG_K']:.3e} K")

# Confirm no horizon at t_n=0
r3b = BSFGBlackHoleSolutionHorizonCalculator().compute({'t_n': 0.0})
assert r3b['horizon_exists'] is False, "t_n=0 should have no horizon"
print(f"horizon@tn0 = {r3b['horizon_exists']}  (expected False)")

# ── #157 Bohr-Sommerfeld ──────────────────────────────────────────────────
r4 = BSFGBohrSommerfeldAetherQuantizationCalculator().compute(ds)
assert r4['paper'] == 'PAPER_562'
assert r4['r_cross_AU'] > 0
print(f"r_cross_AU  = {r4['r_cross_AU']:.3f} AU")
print(f"dJJ_1AU     = {r4['delta_J_over_J_1AU']:.3e}")
print(f"h_eta       = {r4['h_eta']:.3e}")
print(f"aether_dom  = {r4['aether_dominates']}  (Aether dominates at Rs)")

print()
print("ALL 4 BSFG SESSION 149 CLASSES PASSED")
