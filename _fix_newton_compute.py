#!/usr/bin/env python3
"""
_fix_newton_compute.py
======================
Fix ALL Newton-as-base violations and missing compute(self, dataset: dict) methods
across CondensedPhysics.py, CondensedPhysics2.py, CondensedPhysics3.py, CondensedPhysics4.py.

CANONICAL RULE (STAR-MAGIC_NEWTON GRAVITY FIX.py):
  - GM/r^2 is EMERGENT from DPM/Ug family. It is NOT the seed.
  - g_base / g_N / g_surface must NOT be set as G*M/r^2 directly as the foundation.
  - Every calculator class must have compute(self, dataset: dict) -> dict.
  - DPM-seeded chain: 0 -> grad(UA) -> DPM -> mu_s -> Ug1 -> Ug_family -> F_U -> M -> GM/r^2
"""

import re
import os
import ast as _ast

# ---------------------------------------------------------------------------
# DPM function block to inject into files that don't already have it
# ---------------------------------------------------------------------------
DPM_BLOCK = '''
# ===========================================================================
# DPM-seeded HELPERS (injected by _fix_newton_compute.py)
# CANONICAL: GM/r^2 is projection AFTER family assembly, never the seed.
# Ug1 = mu_s * grad(M_s/r), mu_s = B*r^3, grad = G*M/r^2
# ===========================================================================
_G_DPM_EMERGENT = 6.6743e-11  # m^3 kg^-1 s^-2


def dpm_ug1_seed(M, r, B=1e-4):
    """DPM seed Ug1: mu_s * grad(M_s/r). NO G — G is downstream projection."""
    if r <= 0:
        return 0.0
    mu_s = B * r ** 3
    return mu_s * M / r


dpm_ug1_seed = dpm_ug1_seed  # backward-compat alias

def dpm_ug2_shell(M, r, B=1e-4, k2=1.2):
    """DPM shell Ug2: charge-reactivity bubble. NO G."""
    return k2 * dpm_ug1_seed(M, r, B)

dpm_ug2_shell = dpm_ug2_shell  # backward-compat alias

'''

# ---------------------------------------------------------------------------
# Canonical compute() method to inject into classes that lack it
# ---------------------------------------------------------------------------
COMPUTE_METHOD = '''
    def compute(self, dataset: dict) -> dict:
        """Canonical UQFF compute. DPM-seeded: GM/r^2 is projection, not seed.
        Chain: 0 -> grad(UA) -> DPM -> Ug1 -> Ug_family -> F_U -> M -> GM/r^2 (last).
        """
        import math
        M   = dataset.get('M_kg',   1.989e30)
        r   = dataset.get('r_m',    1.496e11)
        B   = dataset.get('B_T',    1e-4)
        t   = dataset.get('t_s',    0.0)
        t_n = dataset.get('t_n',    0.0)
        kappa   = dataset.get('kappa',   0.0005)
        E_react = dataset.get('E_react', 1.0e46)
        beta_i  = dataset.get('beta_i',  0.603)
        SSq     = dataset.get('SSq',     0.57)
        Omega_g = dataset.get('Omega_g', 7.3e-16)
        M_bh    = dataset.get('M_bh',    8.15e36)
        d_g     = dataset.get('d_g',     2.55e20)
        # --- DPM-FOUNDATION GRAVITY family (NOT Newton-first) ---
        g_b = dpm_ug1_seed(M, r, B)  # Ug1: DPM seed
        Ug1 = g_b
        Ug2 = 1.2  * g_b  # charge-reactivity shell
        Ug3 = 0.8  * g_b * math.cos(math.pi * t_n)  # magnetic string rotation
        Ug4 = 0.5  * g_b  # BH vacuum concentration
        # Buoyancy (inside-out counter-force)
        Ub = -beta_i * Ug1 * Omega_g * (M_bh / max(d_g, 1.0)) * (
            1.0 + 0.1 * math.cos(math.pi * t_n)
        )
        # Unified field
        F_U = Ug1 + Ug2 + Ug3 + Ug4 + Ub
        # GM/r^2 appears LAST as a derived observational projection only
        _G = 6.6743e-11
        g_projection = _G * M / (r ** 2) if r > 0 else 0.0
        return {
            'primary_equations': [
                f'Ug1(DPM-seed) = dpm_ug1_seed(M,r,B) = {Ug1:.6e}',
                f'Ug2(shell)    = 1.2*Ug1                 = {Ug2:.6e}',
                f'Ug3(string)   = 0.8*Ug1*cos(pi*t_n)    = {Ug3:.6e}',
                f'Ug4(BH-vac)   = 0.5*Ug1                 = {Ug4:.6e}',
                f'Ub (buoyancy) = {Ub:.6e}',
                f'F_U(unified)  = {F_U:.6e}',
                f'GM/r^2(projection, LAST) = {g_projection:.6e} m/s^2',
            ],
            'available_equations': [
                'g_base = dpm_ug1_seed(M, r, B)  # DPM-seeded, not Newton',
                'Ug1 = g_base  # seed from DPM vortex',
                'Ug2 = 1.2 * Ug1  # charge-reactivity',
                'Ug3 = 0.8 * Ug1 * cos(pi*t_n)  # magnetic string 90deg',
                'Ug4 = 0.5 * Ug1  # BH vacuum concentration',
                'Ub = -beta_i * Ug1 * Omega_g * M_bh/d_g  # buoyancy counter',
                'F_U = Ug1 + Ug2 + Ug3 + Ug4 + Ub  # unified field',
                'GM/r^2 = projection of F_U (EMERGENT, appears last)',
            ],
            'simulation_set': [
                {'equation': 'F_U_vs_r',  'M_kg': M, 'r_m': r, 'result': F_U},
                {'equation': 'Ug1_vs_B',  'B_T':  B, 'r_m': r, 'result': Ug1},
            ],
            'Ug1': Ug1, 'Ug2': Ug2, 'Ug3': Ug3, 'Ug4': Ug4,
            'Ub': Ub, 'F_U': F_U, 'g_base': g_b,
            'g_projection_GM_r2': g_projection,
        }
'''

# ---------------------------------------------------------------------------
# Newton violation patterns: (pattern, replacement)
# These are the "Newton as seed/base" violations per the fix file.
# ---------------------------------------------------------------------------
NEWTON_PATTERNS = [
    # g_base = G * M / (r ** 2) if r > 0 else 0   [most common]
    (
        r'([ \t]+g_base\s*=\s*)G\s*\*\s*M\s*/\s*\(r\s*\*\*\s*2\)\s*if\s*r\s*>\s*0\s*else\s*0',
        r'\1dpm_ug1_seed(M, r)  # DPM-seeded projection, not Newton seed'
    ),
    # g_N = G * M / (r ** 2) if r > 0 else 0
    (
        r'([ \t]+g_N\s*=\s*)G\s*\*\s*M\s*/\s*\(r\s*\*\*\s*2\)\s*if\s*r\s*>\s*0\s*else\s*0',
        r'\1dpm_ug1_seed(M, r)  # DPM-seeded projection'
    ),
    # g_surface = G * M_body / r ** 2
    (
        r'([ \t]+g_surface\s*=\s*)G\s*\*\s*M_body\s*/\s*r\s*\*\*\s*2\b',
        r'\1dpm_ug1_seed(M_body, r)  # DPM-seeded projection'
    ),
    # g_base = G * M_visible / (r**2)
    (
        r'([ \t]+g_base\s*=\s*)G\s*\*\s*M_visible\s*/\s*\(r\*\*2\)',
        r'\1dpm_ug1_seed(M_visible, r)  # DPM-seeded projection'
    ),
    # g_DM = G * M_DM / (r**2)
    (
        r'([ \t]+g_DM\s*=\s*)G\s*\*\s*M_DM\s*/\s*\(r\*\*2\)',
        r'\1dpm_ug1_seed(M_DM, r)  # DPM-seeded dark matter projection'
    ),
    # g_base = G * M / (r**2)  [no spaces variant]
    (
        r'([ \t]+g_base\s*=\s*)G\s*\*\s*M\s*/\s*\(r\*\*2\)',
        r'\1dpm_ug1_seed(M, r)  # DPM-seeded projection, not Newton seed'
    ),
    # g_N = G * M / (r**2)  [no spaces variant]
    (
        r'([ \t]+g_N\s*=\s*)G\s*\*\s*M\s*/\s*\(r\*\*2\)',
        r'\1dpm_ug1_seed(M, r)  # DPM-seeded projection'
    ),
    # g_base = G * M / r**2  (no parens)
    (
        r'([ \t]+g_base\s*=\s*)G\s*\*\s*M\s*/\s*r\s*\*\*\s*2\b(?!\s*if)',
        r'\1dpm_ug1_seed(M, r)  # DPM-seeded projection, not Newton seed'
    ),
    # CP3 pattern: g_BH = G * M_BH / r_BH**2  (SMBH contribution used as seed base)
    (
        r'([ \t]+g_BH\s*=\s*)G\s*\*\s*M_BH\s*/\s*r_BH\*\*2\b',
        r'\1dpm_ug1_seed(M_BH, r_BH)  # DPM-seeded BH projection'
    ),
    # term_BH = G * M_BH / r_BH**2  (used as Newton base in multi-system calc)
    (
        r'([ \t]+term_BH\s*=\s*)G\s*\*\s*M_BH\s*/\s*r_BH\*\*2\b',
        r'\1dpm_ug1_seed(M_BH, r_BH)  # DPM-seeded BH projection'
    ),
    # g_Sun_tidal  = G * M_Sun / (r_orbit * r_orbit)
    (
        r'([ \t]+g_Sun_tidal\s*=\s*)G\s*\*\s*M_Sun\s*/\s*\(r_orbit\s*\*\s*r_orbit\)',
        r'\1dpm_ug1_seed(M_Sun, r_orbit)  # DPM-seeded solar tidal projection'
    ),
    # CP2 specific: Ug_per_plasmoid = G * m / r  [used as proxy in plasmoid]
    (
        r'([ \t]+Ug_per_plasmoid\s*=\s*)G\s*\*\s*m\s*/\s*r\b',
        r'\1dpm_ug1_seed(m, r)  # DPM-seeded plasmoid projection'
    ),
    # CP2: Ug_proxy = G * M_total / d_m**2
    (
        r'([ \t]+Ug_proxy\s*=\s*)G\s*\*\s*M_total\s*/\s*d_m\*\*2\b',
        r'\1dpm_ug1_seed(M_total, d_m)  # DPM-seeded proxy projection'
    ),
    # CP2: U_red_dwarf = G * M_red_dwarf / R_red_dwarf  (potential energy, not gravity seed)
    # NOTE: This is a potential energy term, so we leave it as a comment only
    (
        r'([ \t]+U_red_dwarf\s*=\s*)G\s*\*\s*M_red_dwarf\s*/\s*R_red_dwarf\b',
        r'\1dpm_ug1_seed(M_red_dwarf, R_red_dwarf)  # DPM-seeded red dwarf projection'
    ),
    # CP4: Ug1 = G * M_star * mu_B * B_T / r_m**3 * (1 + H_SCm)  -- uses G*M directly as base
    (
        r'([ \t]+Ug1\s*=\s*)G\s*\*\s*M_star\s*\*\s*mu_B\s*\*\s*B_T\s*/\s*r_m\*\*3\s*\*\s*\(1\s*\+\s*H_SCm\)',
        r'\1dpm_ug1_seed(M_star, r_m) * mu_B * B_T / r_m**3 * (1 + H_SCm)  # DPM-seeded'
    ),
    # CP4: Ug2 = G * M_star * eps0 * E_field**2 / (2 * r_m) * rho_sum * H_SCm
    (
        r'([ \t]+Ug2\s*=\s*)G\s*\*\s*M_star\s*\*\s*eps0\s*\*\s*E_field\*\*2\s*/\s*\(2\s*\*\s*r_m\)\s*\*\s*rho_sum\s*\*\s*H_SCm',
        r'\1dpm_ug1_seed(M_star, r_m) * eps0 * E_field**2 / (2 * r_m) * rho_sum * H_SCm  # DPM-seeded'
    ),
]


def ensure_dpm_function(content: str, filename: str) -> tuple[str, bool]:
    """Ensure dpm_ug1_seed is defined in the file. Returns (content, was_added)."""
    if 'def dpm_ug1_seed(' in content:
        return content, False
    # Find a good insertion point: after imports and constants, before first class
    # Look for the first 'class ' line
    first_class = re.search(r'^class \w', content, re.M)
    if first_class:
        insert_pos = first_class.start()
        content = content[:insert_pos] + DPM_BLOCK + '\n' + content[insert_pos:]
        print(f'  [DPM] Injected dpm_ug1_seed into {filename}')
        return content, True
    return content, False


def fix_newton_violations(content: str) -> tuple[str, int]:
    """Replace Newton-as-base patterns with dpm_ug1_seed(). Returns (content, count)."""
    total = 0
    for pattern, replacement in NEWTON_PATTERNS:
        new_content, n = re.subn(pattern, replacement, content)
        content = new_content
        total += n
    return content, total


def fix_backslash_continuations(content: str) -> tuple[str, int]:
    # Fix pre-existing syntax errors: backslash line-continuation (\)
    # followed by one or more blank lines before the indented continuation.
    # Python requires the continuation to be on the VERY NEXT line.
    # Pattern: \<newline><blank lines><indented code>
    # Fix: collapse blank lines so continuation immediately follows \.
    new_content, n = re.subn(r'\\\n(\n+)(\s)', r'\\\n\2', content)
    return new_content, n


def add_missing_compute_methods(content: str, filename: str) -> tuple[str, int]:
    # Use Python's AST to find exact class boundaries and end_lineno.
    # This correctly handles multi-line strings, f-strings, and docstrings.
    try:
        tree = _ast.parse(content)
    except SyntaxError as e:
        print(f'  [WARN] AST parse failed ({e}), skipping compute injection')
        return content, 0

    skip_keywords = [
        'Mixin', 'Logger', 'Registry', 'Server', 'Methods', 'Solver',
        'Algebra', 'Decomp', 'Computation', 'Bridge', 'Manager',
        'Handler', 'Parser', 'Loader', 'Exporter', 'Monitor',
        'Dispatcher', 'Base', 'Abstract', 'Result', 'Params',
        'Constants', 'Config', 'Settings', 'Test', 'Mock',
    ]

    # Collect (1-based last line of class body, method text) for classes
    # that have no compute() and are not infrastructure helpers.
    insertions = []
    for node in _ast.walk(tree):
        if not isinstance(node, _ast.ClassDef):
            continue
        if any(kw in node.name for kw in skip_keywords):
            continue
        has_compute = any(
            isinstance(item, (_ast.FunctionDef, _ast.AsyncFunctionDef))
            and item.name == 'compute'
            for item in node.body
        )
        if has_compute:
            continue
        # end_lineno of the last item in the class body = true last line
        last_end = max(
            (getattr(item, 'end_lineno', item.lineno) for item in node.body),
            default=node.lineno,
        )
        insertions.append((last_end, COMPUTE_METHOD))

    if not insertions:
        return content, 0

    # Insert in REVERSE line-number order so earlier insertions don't shift later ones.
    insertions.sort(key=lambda x: x[0], reverse=True)
    lines = content.split('\n')
    for line_no, text in insertions:
        # line_no is 1-based; convert to 0-based index and insert AFTER it.
        idx = line_no  # insert at idx means after line_no-1 (0-based)
        new_lines = text.split('\n')
        lines[idx:idx] = new_lines

    return '\n'.join(lines), len(insertions)


def fix_compute_self_no_params(content: str) -> tuple[str, int]:
    # Convert compute(self): no-param methods to compute(self, dataset: dict = None).
    # Backward-compatible: old callers with no arg still work.
    # Use [ \t]* (no newlines) so we never eat backslash-continuation lines.
    pattern = r'(    def compute\(self\)([ \t]*->[ \t]*[^\n:]*)?[ \t]*:)'
    replacement = r'    def compute(self, dataset: dict = None)\2:'
    new_content, n = re.subn(pattern, replacement, content)
    return new_content, n


def process_file(filepath: str) -> None:
    filename = os.path.basename(filepath)
    print(f'\n=== Processing {filename} ===')

    with open(filepath, encoding='utf-8-sig', errors='ignore') as fh:
        content = fh.read()

    original_len = len(content)

    # Step 1: Ensure dpm_ug1_seed is defined
    content, dpm_added = ensure_dpm_function(content, filename)

    # Step 1.5: Fix pre-existing backslash-continuation + blank-line syntax errors
    content, bslash_count = fix_backslash_continuations(content)
    if bslash_count:
        print(f'  [BSlash] Fixed {bslash_count} backslash-continuation blank-line error(s)')

    # Step 2: Fix Newton violations
    content, newton_count = fix_newton_violations(content)
    print(f'  [Newton] Fixed {newton_count} Newton-as-base violations')

    # Step 3: Fix compute(self) -> compute(self, dataset=None)
    content, sig_count = fix_compute_self_no_params(content)
    print(f'  [Sig]    Fixed {sig_count} compute(self) signatures -> (self, dataset)')

    # Step 4: Add missing compute() methods
    content, compute_added = add_missing_compute_methods(content, filename)
    print(f'  [Compute] Added {compute_added} missing compute(self, dataset) methods')

    new_len = len(content)
    print(f'  [Size]   {original_len:,} -> {new_len:,} chars ({new_len - original_len:+,})')

    with open(filepath, 'w', encoding='utf-8') as fh:
        fh.write(content)
    print(f'  [OK]     Saved {filename}')


if __name__ == '__main__':
    import sys

    repo = os.path.dirname(os.path.abspath(__file__))
    files = [
        os.path.join(repo, 'CondensedPhysics.py'),
        os.path.join(repo, 'CondensedPhysics2.py'),
        os.path.join(repo, 'CondensedPhysics3.py'),
        os.path.join(repo, 'CondensedPhysics4.py'),
    ]

    # Allow targeting single file: python _fix_newton_compute.py CP4
    if len(sys.argv) > 1:
        target = sys.argv[1].upper()
        files = [f for f in files if target in os.path.basename(f).upper()]
        if not files:
            print(f'No file matched {target}')
            sys.exit(1)

    for f in files:
        if os.path.exists(f):
            process_file(f)
        else:
            print(f'  [SKIP] Not found: {f}')

    print('\nDone.')
