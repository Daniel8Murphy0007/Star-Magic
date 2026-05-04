#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_fix_docstring_imports.py
Move the dpm_vacuum_manifold TOP block from inside docstrings to actual code.
"""
import re, shutil

IMPORT_BLOCK = (
    'from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc\n'
    '_RHO_VAC_SCM, _ = _derive_qc(n_levels=26, f_SCm=0.57)   '
    '# J/m\u00b3 SCm energy density\n'
    '_RHO_VAC_UA,  _ = _derive_qc(n_levels=26, f_SCm=5.7)    '
    '# J/m\u00b3 UA  energy density (10x)\n'
)

# The full embedded-in-docstring block to remove
EMBEDDED_PAT = re.compile(
    r'# \u2500\u2500 Quantum Chain Derived Constants \(UQFF_THEORY\.md\) '
    r'\u2500+\n'
    r'# Vacuum density is emergent energy density J/m\u00b3, NOT kg/m\u00b3\.\n'
    r'# SCm and UA are MASSLESS geometric substrates derived from '
    r'26-level H-atom geometry\.\n'
    r'# All functions that use _RHO_VAC_SCM / _RHO_VAC_UA are automatically correct\.\n'
    r'from dpm_vacuum_manifold import derive_from_quantum_chain as _derive_qc\n'
    r'_RHO_VAC_SCM, _ = _derive_qc\(n_levels=26, f_SCm=0\.57\).*?\n'
    r'_RHO_VAC_UA,  _ = _derive_qc\(n_levels=26, f_SCm=5\.7\).*?\n'
    r'# \u2500+\n',
    re.DOTALL
)


def fix_file(path):
    with open(path, encoding='utf-8', errors='replace') as f:
        content = f.read()

    # Step 1: Check if import block is already in real code (not in docstring)
    # Find position of first closing """ (the end of the module docstring)
    # The module docstring starts at the first """ that appears early in file.
    first_triple = content.find('"""')
    if first_triple == -1:
        print(f'  SKIP {path}: no triple-quote found')
        return

    second_triple = content.find('"""', first_triple + 3)
    if second_triple == -1:
        print(f'  SKIP {path}: no closing triple-quote found')
        return

    docstring_end = second_triple + 3  # position after closing """

    # Check if import block is before the docstring end (i.e., inside docstring)
    import_pos = content.find('from dpm_vacuum_manifold import derive_from_quantum_chain')
    if import_pos == -1:
        print(f'  SKIP {path}: no dpm import found')
        return

    if import_pos > docstring_end:
        print(f'  OK   {path}: import already outside docstring (pos {import_pos} > {docstring_end})')
        return

    print(f'  FIX  {path}: import at pos {import_pos} is inside docstring (ends at {docstring_end})')

    # Step 2: Remove embedded block from docstring
    new_content = EMBEDDED_PAT.sub('', content)
    if new_content == content:
        # Pattern didn't match, try simpler approach - find and remove just those lines
        print(f'    Pattern not matched - trying line-by-line removal')
        lines = content.split('\n')
        new_lines = []
        skip_next = False
        for i, line in enumerate(lines):
            if 'Quantum Chain Derived Constants (UQFF_THEORY.md)' in line:
                # Skip this and next 8 lines (the full embedded block)
                for j in range(9):
                    if i + j < len(lines):
                        pass  # mark for skip
                # Actually just do bulk removal
                break
            new_lines.append(line)
        # Fall back: direct string surgery
        marker = '# \u2500\u2500 Quantum Chain Derived Constants (UQFF_THEORY.md)'
        idx = content.find(marker)
        if idx != -1:
            # Find the line that ends the block (the ─────── line after UA import)
            # Look for the next non-comment, non-from/import line after the block
            end_marker = '# \u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500'
            end_idx = content.find(end_marker, idx)
            if end_idx != -1:
                end_idx = content.find('\n', end_idx) + 1  # include newline
                new_content = content[:idx] + content[end_idx:]
                print(f'    Removed block via string indices')
            else:
                print(f'    ERROR: could not find end of embedded block')
                return
        else:
            print(f'    ERROR: could not find start of embedded block')
            return

    # Step 3: Now find the recalculated docstring end and insert import block after it
    first_triple2 = new_content.find('"""')
    if first_triple2 == -1:
        print(f'    ERROR: no triple-quote in modified content')
        return
    second_triple2 = new_content.find('"""', first_triple2 + 3)
    if second_triple2 == -1:
        print(f'    ERROR: no closing triple-quote in modified content')
        return

    insert_pos = second_triple2 + 3  # right after closing """
    # Insert: \n\n{IMPORT_BLOCK}\n
    new_content = new_content[:insert_pos] + '\n\n' + IMPORT_BLOCK + new_content[insert_pos:]

    # Backup and write
    shutil.copy2(path, path + '.bak2')
    with open(path, 'w', encoding='utf-8') as f:
        f.write(new_content)
    print(f'    Written.')


files = [
    'CondensedPhysics.py',
    'CondensedPhysics2.py',
    'CondensedPhysics3.py',
    'CondensedPhysics4.py',
]

for f in files:
    fix_file(f)

print('\nDone. Running syntax check...')
import py_compile
for f in files:
    try:
        py_compile.compile(f, doraise=True)
        print(f'  OK  {f}')
    except py_compile.PyCompileError as e:
        print(f'  ERR {f}: {e}')
