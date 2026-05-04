# -*- coding: utf-8 -*-
"""Apply Quantum Chain rules to entire codebase.
- Replaces 7.09e-37 / 7.09e-36 hardcoded literals with named constants derived
  from derive_from_quantum_chain() in each Python file.
- Adds derive_from_quantum_chain import block at top of each CP/QCalc file.
- Fixes index.js to use JS deriveFromQuantumChain() function.
- Fixes whitepaper calibration table unit labels kg/m3 -> J/m3 for rho_SCm/UA.
"""

from pathlib import Path
import re

# ── Canonical derived values (mirrors scm_vacuum_manifold.derive_from_quantum_chain) ──
E0 = 1e-20
f_SCm = 0.57
n_levels = 26
E_n = [E0 * 10**n for n in range(1, n_levels + 1)]
RHO_VAC_SCM_QC = sum(f_SCm * E for E in E_n)        # J/m³ — 6.333e+05
RHO_VAC_UA_QC  = sum(f_SCm * 10 * E for E in E_n)   # J/m³ — 6.333e+06

PREAMBLE_PY = '''\
# ── Quantum Chain Derived Constants (UQFF_THEORY.md) ──────────────────────────
# Vacuum density is emergent energy density J/m³, NOT kg/m³.
# SCm and UA are MASSLESS geometric substrates derived from 26-level H-atom geometry.
# All functions that use _RHO_VAC_SCM / _RHO_VAC_UA are automatically correct.
try:
    from scm_vacuum_manifold import derive_from_quantum_chain as _derive_qc
    _RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)   # J/m³ SCm energy density
    _RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)    # J/m³ UA  energy density (10x)
except ImportError:
    # Fallback: canonical numeric values if scm_vacuum_manifold not on path
    _RHO_VAC_SCM = {:.6e}   # J/m³ — SCm vacuum energy density (Quantum Chain)
    _RHO_VAC_UA  = {:.6e}   # J/m³ — UA  vacuum energy density (Quantum Chain)
# ─────────────────────────────────────────────────────────────────────────────
'''.format(RHO_VAC_SCM_QC, RHO_VAC_UA_QC)


def patch_python_file(path: Path) -> int:
    """Replace hardcoded 7.09e-37/7.09e-36 literals in a Python file.
    Returns count of replacements made."""
    txt = path.read_text('utf-8', errors='replace')

    # Check if already patched
    if '_RHO_VAC_SCM, _ = _derive_qc' in txt:
        print(f'  {path.name}: already patched — skipping')
        return 0

    # Count occurrences
    n37 = txt.count('7.09e-37')
    n36 = txt.count('7.09e-36')
    if n37 == 0 and n36 == 0:
        print(f'  {path.name}: no literals found — skipping')
        return 0

    # Replace literals with named constants
    txt2 = txt.replace('7.09e-37', '_RHO_VAC_SCM')
    txt2 = txt2.replace('7.09e-36', '_RHO_VAC_UA')

    # Insert preamble after any existing shebang/encoding line and imports block
    # Find the first non-comment, non-blank line that is an import or after imports
    lines = txt2.split('\n')
    insert_at = 0
    # Skip initial encoding comments / docstrings
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith('#') or stripped == '' or stripped.startswith('"""') or stripped.startswith("'''"):
            insert_at = i + 1
            continue
        # First real code line
        insert_at = i
        break

    lines.insert(insert_at, PREAMBLE_PY)
    txt_out = '\n'.join(lines)

    path.write_text(txt_out, 'utf-8')
    print(f'  {path.name}: patched {n37} × 7.09e-37  +  {n36} × 7.09e-36  (inserted at line {insert_at+1})')
    return n37 + n36


def patch_index_js(path: Path) -> int:
    """Add JS deriveFromQuantumChain() and replace hardcoded literals in index.js."""
    txt = path.read_text('utf-8', errors='replace')
    if 'deriveFromQuantumChain' in txt:
        print(f'  {path.name}: already patched — skipping')
        return 0

    n37 = txt.count('7.09e-37')
    n36 = txt.count('7.09e-36')

    js_preamble = '''\
// ── Quantum Chain Derived Constants (UQFF_THEORY.md) ─────────────────────────
// Vacuum density is emergent ENERGY density J/m³ (NOT kg/m³).
// SCm and UA are MASSLESS geometric substrates (26D H-atom geometry).
function deriveFromQuantumChain(nLevels = 26, fSCm = 0.57, V = 1.0) {
    const E0 = 1e-20; // J — base energy scale
    let rhoVacEnergy = 0;
    for (let n = 1; n <= nLevels; n++) {
        rhoVacEnergy += fSCm * (E0 * Math.pow(10, n));
    }
    rhoVacEnergy /= V; // J/m³
    const c = 2.99792458e8;
    return { rhoVacEnergy, rhoMassEq: rhoVacEnergy / (c * c) };
}
const { rhoVacEnergy: _RHO_VAC_SCM } = deriveFromQuantumChain(26, 0.57);  // J/m³
const { rhoVacEnergy: _RHO_VAC_UA  } = deriveFromQuantumChain(26, 5.7);   // J/m³ (10× SCm)
// ─────────────────────────────────────────────────────────────────────────────

'''
    txt2 = txt.replace('7.09e-37', '_RHO_VAC_SCM')
    txt2 = txt2.replace('7.09e-36', '_RHO_VAC_UA')

    # Insert after first comment block (before first 'const' or 'function' or 'var')
    lines = txt2.split('\n')
    insert_at = 0
    for i, line in enumerate(lines):
        s = line.strip()
        if s.startswith('//') or s == '':
            insert_at = i + 1
        else:
            break
    lines.insert(insert_at, js_preamble)
    path.write_text('\n'.join(lines), 'utf-8')
    print(f'  {path.name}: patched {n37} × 7.09e-37  +  {n36} × 7.09e-36')
    return n37 + n36


def fix_whitepaper_units(wp_dir: Path) -> int:
    """Replace kg/m³ with J/m³ for rho_SCm/rho_UA in calibration tables."""
    patterns = [
        # Match table rows: | rho_SCm | ... kg/m³ ... |
        (re.compile(r'(\\rho_\{?SCm\}?[^|]*?\|[^|]*?)kg/m\^\{?3\}?', re.IGNORECASE), r'\1J/m^{3}'),
        (re.compile(r'(\\rho_\{?UA\}?[^|]*?\|[^|]*?)kg/m\^\{?3\}?', re.IGNORECASE), r'\1J/m^{3}'),
        # Plain text in table cells
        (re.compile(r'(\|\s*\\rho_SCm\s*\|[^|]*?\|[^|]*?)kg/m\^3'), r'\1J/m^3'),
        (re.compile(r'(\|\s*\\rho_\{SCm\}\s*\|[^|]*?\|[^|]*?)kg/m\^3'), r'\1J/m^3'),
        # Calibration table value rows with 7.09e-37 kg/m³
        (re.compile(r'7\.09e-37\s+kg/m\^?\{?3\}?'), '7.09e-37 J/m^{3}'),
        (re.compile(r'7\.09e-36\s+kg/m\^?\{?3\}?'), '7.09e-36 J/m^{3}'),
        # | 7.09e-37 | kg/m³ | pattern
        (re.compile(r'(7\.09e-3[67][^|]*?\|[^|]*?)kg/m\^?\{?3\}?'), r'\1J/m^{3}'),
    ]

    fixed = 0
    for md in sorted(wp_dir.glob('*.md')):
        txt = md.read_text('utf-8', errors='replace')
        txt2 = txt
        for pat, repl in patterns:
            txt2 = pat.sub(repl, txt2)
        if txt2 != txt:
            md.write_text(txt2, 'utf-8')
            fixed += 1
    print(f'  Whitepapers fixed: {fixed}')
    return fixed


if __name__ == '__main__':
    print('=== QUANTUM CHAIN RULES — FULL CODEBASE APPLICATION ===\n')

    print('1. Patching Python CP/QCalc files...')
    py_files = ['CondensedPhysics.py', 'CondensedPhysics2.py', 'CondensedPhysics3.py',
                'CondensedPhysics4.py', 'QCalc.py']
    total_py = 0
    for f in py_files:
        p = Path(f)
        if p.exists():
            total_py += patch_python_file(p)
    print(f'  Total Python replacements: {total_py}\n')

    print('2. Patching index.js...')
    idx = Path('index.js')
    if idx.exists():
        patch_index_js(idx)
    print()

    print('3. Fixing whitepaper calibration unit labels...')
    wp_dir = Path('whitepapers')
    if wp_dir.exists():
        fix_whitepaper_units(wp_dir)
    print()

    print('=== DONE ===')
    print(f'  Derived _RHO_VAC_SCM = {RHO_VAC_SCM_QC:.6e} J/m³')
    print(f'  Derived _RHO_VAC_UA  = {RHO_VAC_UA_QC:.6e}  J/m³')
