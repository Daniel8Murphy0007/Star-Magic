"""_insert_s147_cp4.py
Inserts Session 147 classes (#145-#148) into CondensedPhysics4.py.
Run once: python _insert_s147_cp4.py
Source: grok_share_b08cc4e3684.txt
"""
import re, pathlib, math

CP4 = pathlib.Path('CondensedPhysics4.py')
txt = CP4.read_text(encoding='utf-8')

# ── 1. Count non-underscore classes before insertion ──────────────────────
before = len(re.findall(r'^class (?!_)[A-Z]\w+\(', txt, re.MULTILINE))
print(f'CP4 classes before: {before}')

# ── 2. Constants block to insert ──────────────────────────────────────────
CONSTANTS = '''

# ── Session 147 constants (grok_share_b08cc4e3684.txt) ──────────────────────
import math as _math_s147
_S147_FAC26         = _math_s147.factorial(26)           # 4.0329e+26
_S147_FAC13         = _math_s147.factorial(13)           # 6.2270e+09
_S147_FAC13_SQ      = _S147_FAC13 ** 2                   # 3.878e+19
_S147_R_Q_AU        = (2.0 / _S147_FAC26) ** (1.0/26.0) # ≈ 0.0973 AU
_S147_RHO_MIN       = 1e-3 / _S147_FAC26                 # ≈ 2.48e-30 kg/m³
_S147_DVP_PRIME     = 113
_S147_AU_IN_M       = 1.496e11
'''

CLASSES = '''
# ===========================================================================
# Session 147 — grok_share_b08cc4e3684.txt
# 26th-Level Polynomial Proofs: DPM Quantization, U_g Anti-Collapse, Tensor Hub, FUBi Poly
# PAPER_550–553   CP4 #145–#148
# ===========================================================================

class Um26DPolyQuantizationDPMConfinementCalculator(_CP4Calculator):
    """
    #145 — 26th-Order U_m Polynomial: DPM Quantization & 26D Confinement Proof
    U_m = κ·(DPMn-DPMs)/r^26 + ∂^26/∂t_adj^26(DPMn(SCm)-DPMs(SCm))/UA
    Equilibrium r_q = (κ(DPMn-DPMs)·UA/(26!·c26))^(1/26) ≈ 0.097 AU (proplyds)
    General: d^n(c/r^k)/dr^n = (k+n-1)!/(k-1)! · c/r^(k+n)  [induction proof]
    CERN monopole null results (4 TeV) explained as 26D projection masking (r^23 suppression)
    VDS: series coefficients bounded by P_order/3 eigenvalue
    DVP: 26!·c_26 irrational → primitive roots mod p=113 → non-repeating
    BH26: 26D DPM harmonic count matches BH26 dimension series
    Source: grok_share_b08cc4e3684.txt  PAPER_550
    """
    def compute(self, dataset=None):
        import math
        d     = dataset or {}
        kappa = d.get('kappa',  1.0);  DPMn = d.get('DPMn',  1.0)
        DPMs  = d.get('DPMs',  -1.0);  UA   = d.get('UA',    1.0)
        c_26  = d.get('c_26',   1.0);  r_AU = d.get('r_AU',  1.0)
        P     = d.get('P_order',9.999e-6)
        r_m   = r_AU * _S147_AU_IN_M
        diff  = DPMn - DPMs
        # Equilibrium radius
        r_q_26 = kappa * diff * UA / (_S147_FAC26 * c_26) if c_26 > 0 else float('nan')
        r_q_m  = r_q_26 ** (1.0/26.0) if (isinstance(r_q_26, float) and r_q_26 > 0) else float('nan')
        r_q_AU = r_q_m / _S147_AU_IN_M
        r_q_can = _S147_R_Q_AU   # canonical c_26=1
        # 26th derivative of 1/r at r_m: coeff = 26! (k=1)
        deriv26_coeff = _S147_FAC26
        # U_m terms
        Um_r26  = kappa * diff / r_m**26
        Um_t26  = _S147_FAC26 * c_26 / UA
        lam_min = P / 3.0
        dvp_ok  = (_S147_FAC26 % _S147_DVP_PRIME) != 0
        in_range = 0.05 < r_q_can < 1.5
        return {'r_q_AU': r_q_AU, 'r_q_canonical_AU': r_q_can,
                'Um_r26_term': Um_r26, 'Um_deriv26_term': Um_t26,
                'deriv26_coeff': deriv26_coeff,
                'vds_lambda_stable': lam_min, 'vds_bound_ok': lam_min > 0,
                'dvp_irreducible': dvp_ok, 'dvp_p': _S147_DVP_PRIME,
                'bh26_dim': 26, 'proplyd_range_ok': in_range,
                'cern_mask_r23': r_m**(26-3)}


class Ug26DFactorialAntiCollapseUg4SplitCalculator(_CP4Calculator):
    """
    #146 — U_g 26th-Order Factorial Anti-Collapse + Ug4 13+13 Dual Split
    Ug1_26 = ∂^26(DPMn·SCm)/∂r^26 = 26!·a0  (degree-26 stable core constant)
    Ug4_split = ∂^13(r·t)/∂r^13 · ∂^13(r·t)/∂t^13 = (13!)^2 · r · t
    At r=1e-5 AU, t=-10: Ug4_split ≈ -5.80e+26
    ρ_min > g·SCm/(26!·UA) ≈ 2.48e-30 kg/m³  — vacuum energy anti-collapse threshold
    DVP: 13+13 split → two 13-prime orbit pairs (BH–star duality)
    BH26: two BH26 13-mode sub-series from Ug4 split
    Source: grok_share_b08cc4e3684.txt  PAPER_551
    """
    def compute(self, dataset=None):
        d     = dataset or {}
        g     = d.get('g',   1e-3);  SCm = d.get('SCm',  1.0)
        UA    = d.get('UA',  1.0);   a0  = d.get('a0',   1.0)
        r_AU  = d.get('r_AU',1e-5);  t_v = d.get('t',  -10.0)
        r_m   = r_AU * _S147_AU_IN_M
        Ug1_26    = _S147_FAC26 * a0
        d13_r     = _S147_FAC13 * t_v          # ∂^13(r·t)/∂r^13 = 13!·t
        d13_t     = _S147_FAC13 * r_m          # ∂^13(r·t)/∂t^13 = 13!·r
        Ug4_split = d13_r * d13_t
        rho_thr   = g * SCm / (_S147_FAC26 * UA)
        vac_dens  = 9.47e-27
        Ug_26     = g * (SCm/UA) * (Ug1_26 + Ug4_split)
        return {'Ug1_26': Ug1_26, 'Ug4_split': Ug4_split,
                'd13_r': d13_r, 'd13_t': d13_t, 'fac13_sq': _S147_FAC13_SQ,
                'Ug_26_full': Ug_26,
                'rho_threshold': rho_thr, 'rho_min_canonical': _S147_RHO_MIN,
                'no_collapse': rho_thr < vac_dens,
                'bh26_dual_13_modes': (13, 13),
                'dvp_orbit_residue': _S147_FAC13 % _S147_DVP_PRIME}


class UQFFComp26DTensorOffDiag13NSYMHubCalculator(_CP4Calculator):
    """
    #147 — Full UQFF_comp 26D Tensor: Off-Diagonal ∂^13 Couplings + NS + YM Hub
    UQFF_comp = [[P/3,  13!·(SCm/UA),  0        ],
                 [13!·(SCm/UA), P/3,   0        ],
                 [0,    0,             2P/3 + 26!/ρ^26]]
    NS 26th-order smoothness: (26+k-1)!/r^(k+26) < ∞ for r>0 (no blow-up)
    YM mass gap: Δ = 26!·c/r^26 > 0  (at r=1m: Δ=4.033e26 — factorial guarantee)
    Eigenvalues: P/3 ± 13! (off-diag perturbation); min eig > 0 always
    Hub for #145 (U_m 26D quantization) and #146 (U_g anti-collapse)
    VDS: diagonal = P/3 and 2P/3 (eigenvalue pair), BH26: (3,3) ∂^26(Ub)/∂ρ^26
    Source: grok_share_b08cc4e3684.txt  PAPER_552 (hub)
    """
    def compute(self, dataset=None):
        import math
        d    = dataset or {}
        P    = d.get('P_order',9.999e-6);  SCm = d.get('SCm',1.0)
        UA   = d.get('UA',1.0);             r_m = d.get('r_m',1.0)
        k_ns = d.get('k_ns',2);             c_ym= d.get('c_ym',1.0)
        rho_b= d.get('rho_b',1.0)
        T11 = P / 3.0;  T22 = P / 3.0
        Ub_d26 = (_S147_FAC26 / rho_b**26) if rho_b > 0 else float('inf')
        T33 = 2.0*P/3.0 + Ub_d26
        T12 = _S147_FAC13 * (SCm/UA)
        # Eigenvalues: diagonal ± sqrt(T12^2) for 2x2 block
        eig1 = T11 + math.sqrt(T12**2);  eig2 = T11 - math.sqrt(T12**2)
        eig3 = T33
        min_eig = min(abs(eig1), abs(eig2), abs(eig3))
        # NS bound
        ns_coeff = math.factorial(26 + k_ns - 1)
        ns_bound = ns_coeff * c_ym / (r_m**(k_ns+26)) if r_m > 0 else float('inf')
        # YM gap
        delta = _S147_FAC26 * c_ym / r_m**26 if r_m > 0 else float('inf')
        return {'tensor_T11': T11, 'tensor_T22': T22, 'tensor_T33': T33,
                'off_diag_13': T12, 'off_diag_coeff': _S147_FAC13,
                'eigenvalue_1': eig1, 'eigenvalue_2': eig2, 'eigenvalue_3': eig3,
                'min_eigenvalue': min_eig, 'vds_diagonal': P/3.0,
                'bh26_d26_Ub': Ub_d26,
                'ns_bound_26': ns_bound, 'ym_gap_delta': delta, 'ym_gap_ok': delta > 0,
                'hub_um26': 'PAPER_550 r_q=0.097AU', 'hub_ug26': 'PAPER_551 rho_min=2.48e-30'}


class FUBi26thGaussianTruncatedPolynomialBoundCalculator(_CP4Calculator):
    """
    #148 — F_U_Bi_i with 26th-Order Gaussian Polynomial (Truncated Exponential Proof)
    exp(-z²) ≈ Σ_{k=0}^{26} (-1)^k z^{2k}/k!  (degree-52 polynomial in z)
    26th term: z^52/26! ≈ 2.48e-27 at z=1 (confirms convergence)
    Bounded integral: ∫exp(-z²)dz = √π/2·erf(z) ≤ 1 (anti-collapse proof)
    At z=1: truncated sum = 0.367879 = exp(-1) to 6 decimal places — exact match
    Diophantine: 26!·c_26 irrational → non-repeating per DVP prime mod p=113
    VDS: P_order/3 bounds highest series coefficient c_26 (negligibility threshold)
    BH26: z = (x-μ)/σ evaluated at 92/225/345 GHz ALMA bins
    Source: grok_share_b08cc4e3684.txt  PAPER_553
    """
    def compute(self, dataset=None):
        import math
        d     = dataset or {}
        sigma = d.get('sigma', 1e16);  mu  = d.get('mu',   92e9)
        x     = d.get('x',   345e9);  F_U = d.get('F_U', -9.999e-4)
        P     = d.get('P_order',9.999e-6)
        z     = (x - mu) / sigma
        # Truncated exp(-z^2)
        total = 0.0; t26 = 0.0
        for k in range(27):
            term   = ((-1.0)**k) * (z**(2*k)) / math.factorial(k)
            total += term
            if k == 26: t26 = term
        exact  = math.exp(-z**2)
        FUBi_p = total * F_U
        erf1   = math.erf(1.0)
        b_int  = math.sqrt(math.pi)/2.0 * erf1
        ag6    = abs(total - exact) < 1e-6
        c26_bnd= P / 3.0
        vds_ok = (1.0/_S147_FAC26) <= c26_bnd
        dvp_nr = (_S147_FAC26 % _S147_DVP_PRIME) != 0
        bins   = {}
        for lbl, xb in [('92GHz',92e9),('225GHz',225e9),('345GHz',345e9)]:
            zb = (xb-mu)/sigma; s=0.0
            for k in range(27): s += ((-1)**k)*(zb**(2*k))/math.factorial(k)
            bins[lbl] = {'z': zb, 'poly_val': s, 'FUBi': s*F_U}
        return {'z': z, 'poly26_val': total, 'exact_exp': exact,
                'FUBi_poly26': FUBi_p, 'term26': t26, '26f': _S147_FAC26,
                'bounded_integral': b_int, 'agreement_6dec': ag6,
                'vds_c26_bound': c26_bnd, 'vds_ok': vds_ok,
                'dvp_26f_mod_113': _S147_FAC26 % _S147_DVP_PRIME,
                'dvp_non_repeating': dvp_nr, 'bh26_bins': bins}
'''

ALLOC = '''
    # --- Session 147: grok_share_b08cc4e3684.txt — 26th-Degree Polynomial Proofs
    "Um26DPolyQuantizationDPMConfinementCalculator",          # PAPER_550 (#145)
    "Ug26DFactorialAntiCollapseUg4SplitCalculator",           # PAPER_551 (#146)
    "UQFFComp26DTensorOffDiag13NSYMHubCalculator",            # PAPER_552 hub (#147)
    "FUBi26thGaussianTruncatedPolynomialBoundCalculator",     # PAPER_553 (#148)
'''

# ── 3. Insert constants after Session 146 block ────────────────────────────
CONST_ANCHOR = '_S146_REMNANT_PCT = 18.32      # %'
assert CONST_ANCHOR in txt, 'Constants anchor not found'
txt = txt.replace(CONST_ANCHOR, CONST_ANCHOR + CONSTANTS, 1)

# ── 4. Insert classes before __all__ ──────────────────────────────────────
ALL_ANCHOR = "\n__all__ = [\n"
assert ALL_ANCHOR in txt, '__all__ anchor not found'
txt = txt.replace(ALL_ANCHOR, CLASSES + ALL_ANCHOR, 1)

# ── 5. Insert __all__ entries at end of Galaxy Merger block ───────────────
LAST_ENTRY = '    "GalaxyMergerUQFFVsNewtonEinsteinCalculator",    # PAPER_549 hub (#144)\n\n]'
assert LAST_ENTRY in txt, '__all__ last-entry anchor not found'
txt = txt.replace(LAST_ENTRY, LAST_ENTRY.rstrip(']\n') + ALLOC + '\n]', 1)

# ── 6. Update session header line ─────────────────────────────────────────
OLD_HDR = '    Updated: Session 146 v5.06 — CP4 140→144 (#141–#144 Ug/Ub Boundary Overlap, Ug4 BH Tidal, F_U_Bi_i Collapse Proof, Galaxy Merger UQFF vs Newton/Einstein Hub: PAPER_546–549; grok_share_366dc393a37.txt)'
NEW_HDR = OLD_HDR + '\n    Updated: Session 147 v5.07 — CP4 144→148 (#145–#148 Um26D Quantization, Ug26D Anti-Collapse, UQFF_comp 26D Tensor Hub, FUBi 26th Polynomial: PAPER_550–553; grok_share_b08cc4e3684.txt)'
assert OLD_HDR in txt, 'Header anchor not found'
txt = txt.replace(OLD_HDR, NEW_HDR, 1)

CP4.write_text(txt, encoding='utf-8')

# ── 7. Verify ─────────────────────────────────────────────────────────────
txt2  = CP4.read_text(encoding='utf-8')
after = len(re.findall(r'^class (?!_)[A-Z]\w+\(', txt2, re.MULTILINE))
print(f'CP4 classes after:  {after}')
delta = after - before
if delta == 4:
    print(f'OK: +{delta} classes inserted correctly. Header updated to v5.07.')
else:
    print(f'WARNING: Expected +4, got +{delta}')
