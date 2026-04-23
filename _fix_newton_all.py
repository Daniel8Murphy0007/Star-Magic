#!/usr/bin/env python3
"""
_fix_newton_all.py  (v2 — comprehensive)
==========================================
Replace ALL Newton-as-foundation violations:
  g_base / g_grav / g_BH / Ug1 / ug1 / Ug_base / base_gravity = G * M / r**2
  → dpm_ug1_seed(M, r)   (NO G, NO r**2)

CANONICAL RULE (Star-Magic.txt, immutable):
  0_vacuum → grad(UA) → DPM_vortex → mu_s → Ug1[seed=DPM] → ...→ GM/r^2 [LAST]
  Ug1 seed = dpm_ug1_seed(M, r) = mu_s * M/r   NO G, NO r^2

SKIP (not violations):
  - g_projection / g_proj = _G * M / r**2   ← correct downstream label
  - G * M1 * M2 / r**2                       ← two-body force, not Ug1 seed
  - G * M / max(r**2, ...)                   ← non-standard denominator
  - return (G * M * M) / (C**4 * r) * ...   ← GW strain formula

New vs v1:
  Pattern D — dict-subscript M/r: G * p["M"] / p["r"]**2
  Pattern E — trailing multiplier: g_base = (G*M/(r**2)) * expansion * ...
  Pattern F — return G*M/(r**2) without outer parens
  + g_BH variable name and self.g_BH attribute assignment
  + CP4 added to file list (was incorrectly excluded)
"""

from __future__ import annotations
import re
import os
import ast

# ---------------------------------------------------------------------------
# Atom sub-patterns
# ---------------------------------------------------------------------------

G_ATOM = r'(?:self\.G|G_CONST|_G\b|G\b)'

# M / r expression: simple name, dotted attribute, or dict subscript
# Matches: M, self.M, params.M, p["M"], p['M_BH'], self.params['r_companion']
_SUBSCRIPT = r"""(?:\[["'][\w_]+["']\])"""
ATTR = r"""(?:\w+(?:\.\w+)*""" + _SUBSCRIPT + r"""?)"""

# Variable names that are forbidden Newton seeds (also self.g_BH assignments)
VAR_NAMES = (
    r'(?:'
    r'(?:self\.)?'
    r'(?:ug1|Ug1|Ug2|ug2|Ug_base|g_base|g_grav|g_BH|base_gravity)'
    r')'
)

# ---------------------------------------------------------------------------
# Helper: strip trailing close-paren and if/else guard from tail
# ---------------------------------------------------------------------------

_IF_GUARD = re.compile(
    r'\s*if\s+[\w.]+\s*(?:[><=!]+\s*[\d.e+\-]+)?\s*else\s+0(?:\.0?)?\s*$'
)

def clean_tail(raw: str) -> str:
    """Remove leading ) (outer-paren close) and trailing if-guard from tail."""
    t = raw.strip()
    if t.startswith(')'):
        t = t[1:].strip()
    t = _IF_GUARD.sub('', t).strip()
    return t

# ---------------------------------------------------------------------------
# Pattern A — assignment (all variants incl. dict subscripts + trailing mult.)
#   VAR = [( ] G * M [)] / [(] R**2 [)] [* tail | if ... else 0]
# ---------------------------------------------------------------------------
PAT_A = re.compile(
    r'^(?P<indent>[ \t]+)'
    r'(?P<var>' + VAR_NAMES + r')[ \t]*=[ \t]*'
    r'\(?[ \t]*'                                 # optional outer (
    r'' + G_ATOM + r'[ \t]*\*[ \t]*'            # G *
    r'(?P<M>' + ATTR + r')'                     # M expression (single mass)
    r'[ \t]*\)?[ \t]*/[ \t]*\(?[ \t]*'          # optional ) / optional (
    r'(?P<r>' + ATTR + r')[ \t]*\*\*[ \t]*2'   # r**2
    r'[ \t]*\)?'                                 # optional ) closing denominator (NO newline!)
    r'(?P<tail>[^\n]*)'                          # rest of line (tail or if-guard)
    r'[ \t]*$',
    re.MULTILINE
)

# ---------------------------------------------------------------------------
# Pattern B — k-scaled assignment: VAR = K * G * M / r**2 [* tail]
# ---------------------------------------------------------------------------
PAT_B = re.compile(
    r'^(?P<indent>[ \t]+)'
    r'(?P<var>(?:ug1|Ug1|Ug2|ug2|g_base|Ug_base))[ \t]*=[ \t]*'
    r'(?P<k>(?:self\.\w+(?:\[[\'"]\w+[\'"]\])?|[\w\.]+(?:\[[\'"]\w+[\'"]\])?))[ \t]*\*[ \t]*'
    r'' + G_ATOM + r'[ \t]*\*[ \t]*'
    r'(?P<M>' + ATTR + r')[ \t]*/[ \t]*'
    r'(?P<r>' + ATTR + r')[ \t]*\*\*[ \t]*2'
    r'(?P<tail>[^\n]*)'
    r'[ \t]*$',
    re.MULTILINE
)

# ---------------------------------------------------------------------------
# Pattern RETURN — return G*M/r**2 [* tail] (both plain and paren-wrapped)
# Does NOT match two-mass (M1 * M2 / r**2) or max() denominators.
# ---------------------------------------------------------------------------
PAT_RETURN = re.compile(
    r'^(?P<indent>[ \t]+return[ \t]+)'
    r'\(?[ \t]*'                                 # optional outer (
    r'' + G_ATOM + r'[ \t]*\*[ \t]*'            # G *
    r'(?P<M>' + ATTR + r')'                     # M expression (single mass)
    r'[ \t]*\)?[ \t]*/[ \t]*\(?[ \t]*'          # optional ) / optional (
    r'(?P<r>' + ATTR + r')[ \t]*\*\*[ \t]*2'   # r**2
    r'[ \t]*\)?'                                 # optional ) closing denominator (NO newline!)
    r'(?P<tail>[^\n]*)'                          # rest of line
    r'[ \t]*$',
    re.MULTILINE
)

# ---------------------------------------------------------------------------
# fix_content
# ---------------------------------------------------------------------------

def fix_content(content: str, fname: str) -> tuple[str, int]:
    """Apply all replacement patterns. Returns (new_content, count)."""
    count = 0

    # Pattern B first (most specific — k * G * M / r**2)
    def replace_b(m: re.Match) -> str:
        nonlocal count
        indent = m.group('indent')
        var = m.group('var')
        k = m.group('k')
        M = m.group('M')
        r = m.group('r')
        tail = clean_tail(m.group('tail'))
        count += 1
        base = f"{indent}{var} = {k} * dpm_ug1_seed({M}, {r})"
        if tail and not tail.startswith('#'):
            base += f" {tail}"
        return base

    content = PAT_B.sub(replace_b, content)

    # Pattern A — assignment: VAR = G*M/r**2 in all paren/subscript forms
    def replace_a(m: re.Match) -> str:
        nonlocal count
        indent = m.group('indent')
        var = m.group('var')
        M = m.group('M')
        r = m.group('r')
        tail = clean_tail(m.group('tail'))
        count += 1
        base = f"{indent}{var} = dpm_ug1_seed({M}, {r})"
        if tail and not tail.startswith('#'):
            base += f" {tail}"
        return base

    content = PAT_A.sub(replace_a, content)

    # Pattern RETURN — return G*M/r**2 [* tail]
    def replace_return(m: re.Match) -> str:
        nonlocal count
        indent = m.group('indent')   # includes "return "
        M = m.group('M')
        r = m.group('r')
        tail = clean_tail(m.group('tail'))
        count += 1
        base = f"{indent}dpm_ug1_seed({M}, {r})"
        if tail and not tail.startswith('#'):
            base += f" {tail}"
        return base

    content = PAT_RETURN.sub(replace_return, content)

    return content, count

# ---------------------------------------------------------------------------
# Syntax check
# ---------------------------------------------------------------------------

def verify_syntax(content: str, fname: str) -> bool:
    try:
        ast.parse(content)
        return True
    except SyntaxError as e:
        print(f"  SYNTAX ERROR in {fname}: {e.msg} at line {e.lineno}")
        ctx = content.splitlines()
        for i in range(max(0, e.lineno - 3), min(len(ctx), e.lineno + 2)):
            marker = " <--" if i + 1 == e.lineno else ""
            print(f"    {i+1}: {ctx[i][:100]}{marker}")
        return False

# ---------------------------------------------------------------------------
# Files to fix
# ---------------------------------------------------------------------------

FILES_TO_FIX = [
    'CondensedPhysics.py',
    'CondensedPhysics2.py',
    'CondensedPhysics3.py',
    'CondensedPhysics4.py',           # added: 48+ violations (g_grav dict, g_BH)
    'QCalc.py',
    'QCalc_API.py',
    'QCalc_Advanced.py',
    'QCalc_Performance.py',
    'QCalc_Wolfram_Extensions.py',
    'QCalc_Wolfram_Phase5.py',
    'QCalc_core_uqff.py',
    'QCalc_cpp_equations.py',
    'QCalc_extracted.py',
    'uqff_validation_test.py',
    'Phase6_Consolidated.py',
    'Phase7_Consolidated.py',
    '_fix_newton_compute.py',
    'add_uqff_to_8_models.py',
    'add_uqff_methods.py',
    'add_uqff_remaining.py',
    'add_uqff_v3.py',
    'standard_astrophysics_equations.py',
    'stellar_evolution_module.py',
    'triadic_validations_next.py',
    'grok_100_equations_module.py',
    'grok_100_equations_module_part2.py',
    'millennium_prize_uqff_calculator.py',
    'muge_cluster_3d_sim.py',
    'MUGE_equations_module.py',
    'alpha_clustering_lenr_module.py',
    'compact_objects_module.py',
    'smbh_binary_mergers.py',
    'updated_uqff_2025_module.py',
    'uqff_lagrangian_derivation.py',
    'PhysicsFramework.py',
    'grok_url_calculators.py',
    '99system_master_equation.py',
    'add_master_gravity.py',
    'source81_ngc346_extract.py',
]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

total_fixed = 0
failed = []
skipped = []

for fname in FILES_TO_FIX:
    if not os.path.exists(fname):
        skipped.append(fname)
        continue

    with open(fname, encoding='utf-8-sig', errors='ignore') as f:
        original = f.read()

    new_content, n = fix_content(original, fname)

    if n == 0:
        print(f"  SKIP (no matches): {fname}")
        continue

    if not verify_syntax(new_content, fname):
        print(f"  FAILED syntax check — NOT written: {fname}")
        failed.append(fname)
        continue

    with open(fname, 'w', encoding='utf-8', newline='\n') as f:
        f.write(new_content)

    total_fixed += n
    print(f"  FIXED {n:3d} violations: {fname}")

print(f"\nTotal replacements: {total_fixed}")
print(f"Failed (not written): {failed}")
print(f"Skipped (not found): {skipped}")
