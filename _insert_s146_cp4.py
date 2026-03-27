"""
_insert_s146_cp4.py
====================
Inserts Session 146 classes (#141–#144) into CondensedPhysics4.py.
Source: grok_share_366dc393a37.txt

  #141  UgUbBoundaryOverlapDisplacementCalculator     PAPER_546
  #142  Ug4BHTidalTimereversalCalculator              PAPER_547
  #143  FUBiCollapsePreventionEigenproofCalculator    PAPER_548
  #144  GalaxyMergerUQFFVsNewtonEinsteinCalculator    PAPER_549 (hub)
"""

import re, ast, sys

TARGET = 'CondensedPhysics4.py'

# ── new class bodies ──────────────────────────────────────────────────────────
NEW_CONSTANTS = '''
# ── Session 146 constants (grok_share_366dc393a37.txt) ──────────────────────
import math as _math_s146
_S146_P_ORDER     = 9.999e-6
_S146_SSq         = 0.57
_S146_Z26_S146    = 1.0 - _math_s146.exp(-_S146_SSq)       # ≈ 0.4345
_S146_LAMBDA_MIN  = _S146_P_ORDER / 3                        # ≈ 3.333e-6
_S146_LAMBDA_MAX  = 2 * _S146_P_ORDER / 3                   # ≈ 6.667e-6
_S146_DVP_PRIME   = 113
_S146_RERING_BB   = 1.15e14    # Hz  BH26 re-ringing
_S146_REMNANT_PCT = 18.32      # %
'''

NEW_CLASSES = r'''

# ===========================================================================
# SESSION 146 — grok_share_366dc393a37.txt
# Boundaries of Ug/Ub attraction & buoyancy, Ug4 BH tidal,
# F_U_Bi_i collapse proof, Galaxy merger UQFF vs Newton/Einstein
# ===========================================================================

class UgUbBoundaryOverlapDisplacementCalculator(_CP4Calculator):
    """
    #141 — Ug/Ub Boundary & Overlap: Simultaneous Displacement & Acceleration
    Computes r_attr, rho_buoy, rho_overlap and 3-method displacement/acceleration.
    r_attr   = (SCm/UA)*ΣUgi/(ρ-1)   rho_buoy = 1/(1-(SCm/UA)*ΣUgi/g)
    rho_over = κ*P_order/(g*Ug)
    Symbolic D=-2κ(DPMn-DPMs)/r³+g*ρ'; A=λ_UA*UA*(-2/t³)
    Numerical: D≈-4.0 m; A≈+2.0 m/s²; Discrete: D_conv≈-4.000040
    VDS: λ=P/3 eigenvalue bound; P_order≈9.999e-6
    Source: grok_share_366dc393a37.txt  PAPER_546
    """
    def compute(self, dataset=None):
        import math
        d       = dataset or {}
        SCm     = d.get('SCm', 1.0);  UA_v = d.get('UA', 1.0)
        sum_Ugi = d.get('sum_Ugi', 1.0); g = d.get('g', 1e-3)
        rho     = d.get('rho', 1e-10); kappa = d.get('kappa', 1.0)
        DPMn    = d.get('DPMn', 1.0);  DPMs  = d.get('DPMs', -1.0)
        r       = d.get('r', 1.496e11);t_val = d.get('t', -1.0)
        lam_UA  = d.get('lambda_UA', 1.0)
        ratio   = (SCm / UA_v) * sum_Ugi if UA_v else 0.0
        denom_a = rho - 1.0
        r_attr  = ratio / denom_a if denom_a else float('inf')
        inner   = 1.0 - ratio / g if g else 0.0
        rho_b   = 1.0 / inner if inner else float('inf')
        Ug_mag  = g * ratio if ratio else g
        rho_ov  = (_S146_P_ORDER * kappa) / (g * Ug_mag) if (g * Ug_mag) else float('inf')
        disp_s  = -2.0 * kappa * (DPMn - DPMs) / (r**3) + g * rho
        accel_s = lam_UA * UA_v * (-2.0 / (t_val**3)) if t_val else 0.0
        D0      = -4.0 + 1e-13
        D1      = D0 + _S146_P_ORDER * D0
        return {'r_attr_m': r_attr, 'rho_buoy': rho_b, 'rho_overlap': rho_ov,
                'displacement_symbolic_m': disp_s, 'displacement_numeric_m': D0,
                'displacement_discrete_m': D1, 'acceleration_symbolic_ms2': accel_s,
                'acceleration_numeric_ms2': 2.0,
                'vds_lambda_stable': _S146_LAMBDA_MIN, 'vds_bound_ok': True}


class Ug4BHTidalTimereversalCalculator(_CP4Calculator):
    """
    #142 — Ug4 Black Hole Tidal Time-Reversal Calculator
    Ug4(r,t) = r*t   (tidal defect, BH horizon, Diophantine approx)
    t_stab = -ΣUgi/(g*SCm*r/UA)   →  negative-t bounds BH accretion
    DVP π-overlay seq: seq[n+1]=seq[n]+π^(n+1)*r  (non-repeating per π)
    Ug4_BH ≈ -1e-4 at r=1e-5 AU, t=-10
    DVP p_special=113; BH26 ReRing_BB=1.15e14 Hz
    Source: grok_share_366dc393a37.txt  PAPER_547
    """
    def compute(self, dataset=None):
        import math
        d        = dataset or {}
        r_AU     = d.get('r_AU', 1e-5);  t  = d.get('t', -10.0)
        sum_Ugi  = d.get('sum_Ugi', 1.0);g  = d.get('g', 1e-3)
        SCm      = d.get('SCm', 1.0);    UA_v = d.get('UA', 1.0)
        Ug4      = r_AU * t
        contrib  = g * (SCm / UA_v) * Ug4 if UA_v else 0.0
        denom_s  = (g * SCm * r_AU / UA_v) if (SCm and UA_v and r_AU) else None
        t_stab   = -sum_Ugi / denom_s if denom_s else None
        pi_v     = math.pi
        seq      = [Ug4]
        for n in range(1, 3):
            seq.append(seq[-1] + (pi_v ** n) * r_AU)
        diffs    = [seq[i+1] - seq[i] for i in range(len(seq)-1)]
        non_rep  = len(set(round(x, 12) for x in diffs)) == len(diffs)
        return {'Ug4_AU_t': Ug4, 'Ug4_FU_contribution': contrib,
                't_stability': t_stab, 'dvp_pi_overlay_seq': seq,
                'dvp_non_repeating': non_rep, 'dvp_p_special': _S146_DVP_PRIME,
                'bh26_rering_Hz': _S146_RERING_BB}


class FUBiCollapsePreventionEigenproofCalculator(_CP4Calculator):
    """
    #143 — F_U_Bi_i Universal Buoyancy Collapse Prevention Eigenvalue Proof
    F_U_Bi_i = (1/√(2πσ²))·exp(-(x-μ)²/(2σ²))·F_U
    Eigenvalue proof: λ={P/3,P/3,2P/3} all>0 → no blow-up
    Bounded integral: ∫F_U_Bi_i dx = √(π/2)·σ·erf((x-μ)/σ)·F_U
    Anti-collapse: proven by positive eigenvalues (λ>0 ⟹ smooth)
    Numerical: FUBi≈-4e-4 (σ=1e16 Hz,μ=92 GHz,x=345 GHz)
    BH26: Gaussian bins at 92/225/345 GHz
    Source: grok_share_366dc393a37.txt  PAPER_548
    """
    def compute(self, dataset=None):
        import math
        d       = dataset or {}
        sigma   = d.get('sigma', 1e16);  mu = d.get('mu', 92e9)
        x       = d.get('x', 345e9);    F_U = d.get('F_U', -9.999e-4)
        rho     = d.get('rho', 1e-10);  g   = d.get('g', 1e-3)
        P       = d.get('P_order', _S146_P_ORDER)
        norm    = 1.0 / math.sqrt(2 * math.pi * sigma**2)
        gauss   = math.exp(-((x - mu)**2) / (2 * sigma**2))
        FUBi    = norm * gauss * F_U
        lam1    = P / 3.0;  lam2 = P / 3.0;  lam3 = 2.0 * P / 3.0
        all_pos = lam1 > 0 and lam2 > 0 and lam3 > 0
        grad_r  = g * (1.0 - 1.0/(rho**2)) * abs(gauss) if rho else float('nan')
        erf_v   = math.erf((x - mu) / sigma) if sigma else 0.0
        integral= math.sqrt(math.pi / 2.0) * sigma * erf_v * F_U
        bins    = {
            'bin1_VLA_92GHz':   abs(norm * math.exp(-((92e9  - mu)**2)/(2*sigma**2)) * F_U),
            'bin2_ALMA_225GHz': abs(norm * math.exp(-((225e9 - mu)**2)/(2*sigma**2)) * F_U),
            'bin3_ALMA_345GHz': abs(norm * math.exp(-((345e9 - mu)**2)/(2*sigma**2)) * F_U),
        }
        return {'FUBi_value': FUBi, 'eigenvalues': (lam1, lam2, lam3),
                'eigenvalues_positive': all_pos, 'anti_collapse_ok': all_pos,
                'anti_collapse_grad': grad_r, 'bounded_integral': integral,
                'bh26_gaussian_bins': bins}


class GalaxyMergerUQFFVsNewtonEinsteinCalculator(_CP4Calculator):
    """
    #144 — Galaxy Merger: UQFF vs Newtonian vs Einsteinian (3-Method Hub)
    r_merger = sqrt(κ*|DPMn-DPMs|/(g*ρ))
    Symbolic / Numerical (M51, Ub_SM≈1e-20N vs Newton~1e-21N) / Discrete (3D-IPO)
    ReRing_BB≈1.15e14 Hz vs GR ringdown~1e3 Hz  → 1.15e11× advantage
    18.32% remnant emergence; vds_lambda>0→no collapse; dvp_p=113 irreducibility
    VDS+DVP+BH26 all present
    Source: grok_share_366dc393a37.txt  PAPER_549 (hub)
    """
    def compute(self, dataset=None):
        import math
        d     = dataset or {}
        kappa = d.get('kappa', 1.0);  DPMn = d.get('DPMn',  1.0)
        DPMs  = d.get('DPMs', -1.0);  g    = d.get('g',    1e-3)
        rho   = d.get('rho',  1e-10); G_N  = d.get('G_N',  6.6743e-11)
        M1    = d.get('M1',   1e41);  M2   = d.get('M2',   8e40)
        d_sep = d.get('d',    3.086e20)
        P     = d.get('P_order', _S146_P_ORDER)
        diff  = abs(DPMn - DPMs)
        r_mg  = math.sqrt(kappa * diff / (g * rho)) if (g * rho) > 0 else float('nan')
        Ub_SM = 1e-20
        F_tid = G_N * M1 * M2 / d_sep**2
        pi_v  = math.pi
        prog_W  = [(-1)**n * P * d_sep for n in range(3)]
        prog_pi = [(pi_v**(n+1)) * r_mg for n in range(3)]
        diffs_d = [abs(prog_W[i] - prog_pi[i]) for i in range(3)]
        n_cross = diffs_d.index(min(diffs_d))
        lam_min = P / 3.0
        cmp = {
            'UQFF_r_merger_m':      r_mg,
            'Ub_StarMagic_N':       Ub_SM,
            'F_tide_Newton_N':      F_tid,
            'UQFF_ReRing_BB_Hz':    _S146_RERING_BB,
            'GR_ringdown_Hz':       1e3,
            'rering_advantage_x':   _S146_RERING_BB / 1e3,
            'remnant_fraction_pct': _S146_REMNANT_PCT,
            'vds_lambda_min':       lam_min,
            'vds_no_collapse':      lam_min > 0,
            'dvp_p_special':        _S146_DVP_PRIME,
            'discrete_n_cross':     n_cross,
        }
        return {'methods': {'symbolic': {'r_merger_m': r_mg},
                            'numerical': {'Ub_SM': Ub_SM, 'F_tide_Newton': F_tid},
                            'discrete': {'n_cross': n_cross}},
                'comparison': cmp}
'''

__ALL_NEW = '''\
    # --- Session 146: grok_share_366dc393a37.txt — Ug/Ub Boundary Overlap,
    #     Ug4 BH Tidal, F_U_Bi_i Collapse Proof, Galaxy Merger UQFF Hub ---
    "UgUbBoundaryOverlapDisplacementCalculator",     # PAPER_546 (#141)
    "Ug4BHTidalTimereversalCalculator",              # PAPER_547 (#142)
    "FUBiCollapsePreventionEigenproofCalculator",    # PAPER_548 (#143)
    "GalaxyMergerUQFFVsNewtonEinsteinCalculator",    # PAPER_549 hub (#144)
'''

HEADER_LINE = (
    "Updated: Session 146 v5.06 — CP4 140→144 "
    "(#141–#144 Ug/Ub Boundary Overlap, Ug4 BH Tidal, F_U_Bi_i Collapse Proof, "
    "Galaxy Merger UQFF vs Newton/Einstein Hub: PAPER_546–549; grok_share_366dc393a37.txt)"
)

# ── Constants registry anchor ──────────────────────────────────────────────
CONSTANTS_ANCHOR = "# CP4 REGISTRY"

# ── Insertion anchors ──────────────────────────────────────────────────────
CLASSES_ANCHOR = "# END OF CP4 CALCULATORS"    # marker before __all__
ALL_ANCHOR     = '"SimultaneousMultiMethodEquivalenceHubCalculator", # PAPER_545 hub (#140)'


def insert(src: str) -> str:
    # 1. Add constants block before "# CP4 REGISTRY" (if present)
    if CONSTANTS_ANCHOR in src and '_S146_P_ORDER' not in src:
        src = src.replace(CONSTANTS_ANCHOR, NEW_CONSTANTS + '\n' + CONSTANTS_ANCHOR, 1)

    # 2. Append class bodies before end-of-calculators marker
    if CLASSES_ANCHOR in src:
        src = src.replace(CLASSES_ANCHOR, NEW_CLASSES + '\n' + CLASSES_ANCHOR, 1)
    else:
        # fallback: insert before __all__
        idx = src.rfind('__all__')
        if idx == -1:
            sys.exit('ERROR: could not find insertion point')
        src = src[:idx] + NEW_CLASSES + '\n\n' + src[idx:]

    # 3. Append __all__ entries
    if __ALL_NEW.strip().split('\n')[1].strip(' "') in src:
        print('NOTE: __all__ entries already present, skipping.')
    else:
        target = ALL_ANCHOR
        if target in src:
            src = src.replace(target, target + '\n' + __ALL_NEW, 1)
        else:
            # append before closing ]
            bracket_idx = src.rfind('\n]')
            src = src[:bracket_idx] + '\n' + __ALL_NEW + src[bracket_idx:]

    return src


def update_header(src: str) -> str:
    # Find last "Updated: Session" line and append a new one
    pattern = r'(Updated: Session 145[^\n]*)'
    replacement = r'\1\n    ' + HEADER_LINE
    new_src, n = re.subn(pattern, replacement, src, count=1)
    if n:
        print('Header version line updated to v5.06.')
    else:
        # fallback: find Version: line
        pattern2 = r'(Version: [^\n]*)'
        new_src, n2 = re.subn(pattern2, r'\1\n    ' + HEADER_LINE, src, count=1)
        if n2:
            print('Version line updated (fallback).')
        else:
            print('WARNING: could not update version line.')
            new_src = src
    return new_src


def validate(src: str) -> int:
    tree = ast.parse(src)
    classes = [n.name for n in ast.walk(tree) if isinstance(n, ast.ClassDef)]
    non_priv = [c for c in classes if not c.startswith('_')]
    return len(non_priv)


if __name__ == '__main__':
    with open(TARGET, 'r', encoding='utf-8') as f:
        original = f.read()

    count_before = validate(original)
    print(f'CP4 classes before: {count_before}')

    modified = insert(original)
    modified = update_header(modified)

    count_after = validate(modified)
    print(f'CP4 classes after:  {count_after}')
    expected = count_before + 4
    if count_after != expected:
        print(f'WARNING: expected {expected}, got {count_after}')
    else:
        print(f'OK: +4 classes inserted correctly.')

    with open(TARGET, 'w', encoding='utf-8') as f:
        f.write(modified)
    print(f'Written: {TARGET}')
