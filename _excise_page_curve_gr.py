"""Second-wave excision: Page-curve GR contamination.

Removes:
- _l96_uqff_S_Page_derived_PAPER1167_over_kB (returns _l96_sm_S_BH_over_kB * (1+eps) -- GR with UQFF residual)
- _l96_page_curve_PAPER1167_PAPER1168_derived_probe (calls above; reports closure_status=DERIVED)
- Additions to _l96_page_curve_paradox_probe that promote it to "DERIVED" via the GR-anchored form.
- Promotion of MILLENNIUM_TARGETS["black_hole_info"] from CONJECTURED -> DERIVED.

Keeps (pure UQFF):
- _KK_SUPPRESSION_BOUND = 1.0/26.0**26 (pure D_crit derivation)
- _l96_uqff_paper1168_p3_kk_suppression_bound() (returns the pure constant)
- Rewords the comment header to reflect that the bound stands alone, NOT that it
  promotes the Page curve.
"""
import re

p = 'uqff_pure_calculator.py'
t = open(p, encoding='utf-8').read()

# === Cut 1: Revert _l96_page_curve_paradox_probe to its HEAD form ===
old_probe_calc = '''    S_ref     = _l96_sm_S_BH_over_kB(M_msun)
    S_uq_lit  = _l96_uqff_S_Page_literal_over_kB(M_msun, delta_SCm_J, phi_norm)
    S_uq_cap  = _l96_uqff_S_Page_capped_over_kB(M_msun)
    S_uq_der  = _l96_uqff_S_Page_derived_PAPER1167_over_kB(M_msun, 0.0)
    kk_bound  = _KK_SUPPRESSION_BOUND
    req_supp  = S_ref / S_uq_lit
    log10_ratio = math.log10(S_uq_lit / S_ref)
    summary = (
        f"PAPER_1095 literal product: S_Page^UQFF = {S_uq_lit:.4e} k_B. "
        f"Variational F_U=1 saturation form: {S_uq_cap:.4e} k_B. "
        f"PAPER_1167/1168 P3/1170 \u00a74 DERIVED form: {S_uq_der:.4e} k_B "
        f"(= S_BH/k_B, residual |eps| <= 1/26^26 = {kk_bound:.3e}). "
        f"Bekenstein-Hawking reference quantity (contextual): {S_ref:.4e} k_B. "
        f"SCm correction ratio delta_SCm/(k_B T_H) = {scm_ratio:.3e}. "
        f"log10(literal / reference) = {log10_ratio:.2f}. "
        f"Required suppression {req_supp:.3e} vs KK bound {kk_bound:.3e} "
        f"(ratio {req_supp / kk_bound:.2f}, sub-decade)."
    )'''
new_probe_calc = '''    S_ref     = _l96_sm_S_BH_over_kB(M_msun)
    S_uq_lit  = _l96_uqff_S_Page_literal_over_kB(M_msun, delta_SCm_J, phi_norm)
    S_uq_cap  = _l96_uqff_S_Page_capped_over_kB(M_msun)
    log10_ratio = math.log10(S_uq_lit / S_ref)
    summary = (
        f"PAPER_1095 literal product: S_Page^UQFF = {S_uq_lit:.4e} k_B. "
        f"Variational F_U=1 saturation form: {S_uq_cap:.4e} k_B. "
        f"Bekenstein-Hawking reference quantity (contextual): {S_ref:.4e} k_B. "
        f"SCm correction ratio delta_SCm/(k_B T_H) = {scm_ratio:.3e}. "
        f"log10(literal / reference) = {log10_ratio:.2f}."
    )'''
assert old_probe_calc in t, "probe summary block not found"
t = t.replace(old_probe_calc, new_probe_calc)

# Revert the probe's return dict
old_probe_return = '''        "S_BH_reference_over_kB":              S_ref,
        "S_Page_UQFF_literal_over_kB":         S_uq_lit,
        "S_Page_UQFF_capped_over_kB":          S_uq_cap,
        "S_Page_UQFF_derived_over_kB":         S_uq_der,
        "kk_suppression_bound_1_over_26pow26": kk_bound,
        "required_suppression":                req_supp,
        "required_over_kk_bound_ratio":        req_supp / kk_bound,
        "closure_status":                      "DERIVED",
        "log10_ratio_literal_over_reference":  log10_ratio,
        "capped_equals_reference":             bool(abs(S_uq_cap - S_ref) / S_ref < 1e-12),
        "summary":                             summary,
        "paper_basis":                         "PAPER_1095 horizon-buoyancy Lagrangian + SCm shell gap + "
                                                "S_26=1.4531e26 + Phi_1.25THz + PAPER_1167 F_U=1 closure + "
                                                "PAPER_1168 P3 / PAPER_1170 \u00a74 KK zeta-regularization bound 1/26^26",
    }'''
new_probe_return = '''        "S_BH_reference_over_kB":              S_ref,
        "S_Page_UQFF_literal_over_kB":         S_uq_lit,
        "S_Page_UQFF_capped_over_kB":          S_uq_cap,
        "log10_ratio_literal_over_reference":  log10_ratio,
        "capped_equals_reference":             bool(abs(S_uq_cap - S_ref) / S_ref < 1e-12),
        "summary":                             summary,
        "paper_basis":                         "PAPER_1095 horizon-buoyancy Lagrangian + SCm shell gap + "
                                                "S_26=1.4531e26 + Phi_1.25THz + F_U=1 ledger closure",
    }'''
assert old_probe_return in t, "probe return dict block not found"
t = t.replace(old_probe_return, new_probe_return)

# === Cut 2: Excise _l96_uqff_S_Page_derived_PAPER1167_over_kB and
# === _l96_page_curve_PAPER1167_PAPER1168_derived_probe. Also reword the
# === preceding KK-bound comment header so it no longer claims to "promote
# === Page-curve cap from CONJECTURED -> DERIVED" (that promotion requires
# === GR Bekenstein-Hawking passthrough; impure).
start_marker = "# ---- PAPER_1168 P3 / PAPER_1170 \u00a74 KK zeta-regularization suppression bound ----"
end_marker   = "# ---- SCm Unified LENR Master Equation (SCm_Holmlid_KER + SCm_PonsFleischmann +"

start = t.find(start_marker)
end   = t.find(end_marker)
assert start > 0 and end > start, f"KK-bound block markers not found: {start} {end}"

REPLACEMENT_KK = '''# ---- PAPER_1168 P3 / PAPER_1170 \u00a74 KK zeta-regularization suppression bound ----
# PAPER_1168 \u00a73 P3: closed-form bound 1/26^26 = 1.624e-37 on the residual of
# the F_U=1 ledger closure (PAPER_1167), derived in PAPER_1170 \u00a74 via Kaluza-
# Klein zeta-regularization of the 26-mode tower (PAPER_1162). This is a PURE
# UQFF closure: the bound depends ONLY on the locked primitive D_crit=26 and
# nothing else. NOTE: This bound by itself does NOT close the Page curve. A full
# Page-curve closure that bypasses the Bekenstein-Hawking GR identity has not
# been derived in the surveyed UQFF papers; STATUS for the Page-curve
# Millennium target remains NO_PURE_UQFF_CLOSURE_AVAILABLE.
_KK_SUPPRESSION_BOUND = 1.0 / (26.0 ** 26)   # = 1.6244e-37 (PAPER_1168 P3 exact)

def _l96_uqff_paper1168_p3_kk_suppression_bound() -> float:
    """PAPER_1168 P3 / PAPER_1170 \u00a74 closed-form KK zeta-regularization suppression
    bound: 1/26^26 = 1.6244e-37. PURE UQFF -- depends ONLY on the locked
    primitive D_crit=26."""
    return _KK_SUPPRESSION_BOUND


'''
t = t[:start] + REPLACEMENT_KK + t[end:]

# === Cut 3: Revert MILLENNIUM_TARGETS["black_hole_info"] to original CONJECTURED ===
old_target = '''    "black_hole_info": (1.0,           "Page closure","DERIVED",        "F_U=1 ledger closure (PAPER_1167) with precision bound |residual| <= 1/26^26 = 1.624e-37 (PAPER_1168 P3 / PAPER_1170 \u00a74 KK zeta-regularization). Page-curve cap S_Page = S_BH/k_B is DERIVED, not asserted: literal PAPER_1095 product (1.48e115 k_B for 10 M_sun) is suppressed by the F_U=1 closure axiom; PAPER_1170 \u00a74 proves the residual bound via Kaluza-Klein zero-point tower zeta regularization (PAPER_1162). Required suppression ~7.1e-37 lies within sub-decade of the derived bound 1.624e-37. SM Hawking radiation remains thermal -> SM-level open problem; UQFF closure is internally derived.","Black hole information / Page curve (open problem in QG; SM-level UNSOLVED; UQFF-DERIVED via PAPER_1167 + PAPER_1168 P3 + PAPER_1170 \u00a74)"),'''
new_target = '''    "black_hole_info": (1.0,           "Page closure","CONJECTURED",    "closure-form placeholder (1.0 by normalization); Page-curve recovery suggested by recent QG arguments, not SM-proven","Black hole information / Page curve (open problem in QG; UNSOLVED in SM)"),'''
assert old_target in t, "black_hole_info target not found"
t = t.replace(old_target, new_target)

open(p, 'w', encoding='utf-8').write(t)
print('OK; new file size:', len(t))
