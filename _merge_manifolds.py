#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
_merge_manifolds.py  -- one-time consolidation script

Merges scm_vacuum_manifold.py + ua_vacuum_manifold.py INTO dpm_vacuum_manifold.py.
Run once, then delete this script and the original scm/ua files.

Steps performed:
  1. Strip coding header, module docstring, from __future__ from each sub-file body.
  2. From ua body: remove 'from scm_vacuum_manifold import (...)' block.
  3. From ua body: remove duplicate RHO_VAC_SCM and RHO_VAC_UA definitions.
  4. From dpm body: remove 'from scm_vacuum_manifold import (...)' block.
  5. From dpm body: remove 'from ua_vacuum_manifold import (...)' block.
  6. From dpm body: remove the '# -- Import the two standalone...' comment line.
  7. Write new dpm_vacuum_manifold.py v3.0 as:
        [new header]
        SCm section separator + scm_body
        UA  section separator + ua_body
        DPM section separator + dpm_body
"""

import os
import sys

CWD = os.path.dirname(os.path.abspath(__file__))


def read_text(path):
    with open(path, 'r', encoding='utf-8') as f:
        return f.read()


def write_text(path, content):
    with open(path, 'w', encoding='utf-8') as f:
        f.write(content)


def strip_file_header(content):
    """Strip: # -*- coding -*- line, module docstring, from __future__ import.

    Returns the remaining content as a string (everything after those three
    header elements and any intervening blank lines).
    """
    lines = content.split('\n')
    i = 0
    n = len(lines)

    # 1. Strip coding declaration
    if i < n and lines[i].strip().startswith('# -*-'):
        i += 1

    # 2. Skip blank lines
    while i < n and not lines[i].strip():
        i += 1

    # 3. Strip module docstring (single-line or multi-line)
    if i < n and lines[i].strip().startswith('"""'):
        first = lines[i].strip()
        if first.count('"""') >= 2 and first != '"""':
            # Triple-quoted single-line: """text"""
            i += 1
        else:
            # Multi-line: opening """ is on its own line (or with text but only one """)
            i += 1
            while i < n and '"""' not in lines[i]:
                i += 1
            i += 1  # skip the closing """ line

    # 4. Skip blank lines
    while i < n and not lines[i].strip():
        i += 1

    # 5. Strip 'from __future__ import annotations' if present
    if i < n and lines[i].strip().startswith('from __future__'):
        i += 1

    # 6. Skip blank lines
    while i < n and not lines[i].strip():
        i += 1

    return '\n'.join(lines[i:])


def remove_import_block(content, module_name):
    """Remove a multi-line 'from <module_name> import (\\n...\\n)' block.

    Handles:
      from module import (
          name1,
          name2,
      )
    """
    lines = content.split('\n')
    out = []
    skipping = False
    for line in lines:
        stripped = line.strip()
        # Detect start of the import block
        if not skipping and stripped.startswith(f'from {module_name} import'):
            skipping = True
        if skipping:
            # End of block: a line that is just ')' or ends with ')' and has no comma
            # (catches both standalone ')' and inline single-line 'from x import (a)')
            if stripped == ')' or (stripped.endswith(')') and '(' in stripped):
                skipping = False
            continue  # skip this line
        out.append(line)
    return '\n'.join(out)


def remove_lines_containing(content, *patterns):
    """Remove any line that contains ALL patterns (if multiple given, each triggers removal)."""
    lines = content.split('\n')
    out = []
    for line in lines:
        keep = True
        for pat in patterns:
            if pat in line:
                keep = False
                break
        if keep:
            out.append(line)
    return '\n'.join(out)


# ---------------------------------------------------------------------------
# 1. Read source files
# ---------------------------------------------------------------------------
scm_path = os.path.join(CWD, 'scm_vacuum_manifold.py')
ua_path  = os.path.join(CWD, 'ua_vacuum_manifold.py')
dpm_path = os.path.join(CWD, 'dpm_vacuum_manifold.py')

for p in (scm_path, ua_path, dpm_path):
    if not os.path.exists(p):
        print(f"ERROR: {p} not found", file=sys.stderr)
        sys.exit(1)

scm_raw = read_text(scm_path)
ua_raw  = read_text(ua_path)
dpm_raw = read_text(dpm_path)

print(f"Read scm: {len(scm_raw.splitlines())} lines")
print(f"Read ua:  {len(ua_raw.splitlines())} lines")
print(f"Read dpm: {len(dpm_raw.splitlines())} lines")

# ---------------------------------------------------------------------------
# 2. Process scm body — strip header only
# ---------------------------------------------------------------------------
scm_body = strip_file_header(scm_raw)
print(f"scm_body after strip: {len(scm_body.splitlines())} lines")

# ---------------------------------------------------------------------------
# 3. Process ua body — strip header + remove scm import + remove duplicate constants
# ---------------------------------------------------------------------------
ua_body = strip_file_header(ua_raw)
ua_body = remove_import_block(ua_body, 'scm_vacuum_manifold')
# Remove the duplicate module-level definitions that scm already provides
ua_body = remove_lines_containing(ua_body, 'RHO_VAC_SCM = derive_from_quantum_chain()')
ua_body = remove_lines_containing(ua_body, 'RHO_VAC_UA = derive_from_quantum_chain(')
print(f"ua_body  after strip: {len(ua_body.splitlines())} lines")

# ---------------------------------------------------------------------------
# 4. Process dpm body — strip header + remove scm/ua imports + remove comment
# ---------------------------------------------------------------------------
dpm_body = strip_file_header(dpm_raw)
dpm_body = remove_import_block(dpm_body, 'scm_vacuum_manifold')
dpm_body = remove_import_block(dpm_body, 'ua_vacuum_manifold')
dpm_body = remove_lines_containing(dpm_body, '# -- Import the two standalone manifold layers')
print(f"dpm_body after strip: {len(dpm_body.splitlines())} lines")

# ---------------------------------------------------------------------------
# 5. Build consolidated file
# ---------------------------------------------------------------------------
NEW_HEADER = '''\
# -*- coding: utf-8 -*-
"""
dpm_vacuum_manifold.py  v3.0  (CONSOLIDATED)
Di-Pseudo-Monopole (DPM) Vacuum Calculator -- Quantum Chain Compliant

CONSOLIDATED FILE:
  Absorbs scm_vacuum_manifold.py (SCm base layer) and
           ua_vacuum_manifold.py (UA superstructure) into a single module.
  Legacy root modules scm_vacuum_manifold.py and ua_vacuum_manifold.py are
  restored for compatibility and direct SCm/UA inspection.
  Import dpm_vacuum_manifold.py as the canonical consolidated root.

THE QUANTUM CHAIN (canonical, Star-Magic.txt lines 11-22, IMMUTABLE):

  Step 0  0_vacuum   -> |grad(UA)|          vacuum tension differential
  Step 1  grad(UA)   -> DPM_vortex          a_DPM = F_DPM*f_DPM*E_vac/(c*V_sys)
  Step 2  DPM_vortex -> mu_s                mu_s = rho_A * V_DPM  (vortex volume)
  Step 3  mu_s       -> Ug1[seed=DPM]       Ug1 seeded from mu_s -- NOT from mass
  Step 4  Ug1        -> Ug_family           Ug2+Ug3+Ug4 simultaneously promoted
  Step 5  Ug_family  -> F_U                 + Um + FUBi + FUBii + UA_uv
  Step 6  F_U        -> crossing            FUBi(r) + FUBii(r) = 0  compaction
  Step 7  crossing   -> M_emergent          mass BORN at crossing, not before
  Step 8  M_emergent -> GM/r^2             LAST -- observational projection only

ARCHITECTURE (consolidated):
  Section A  SCm Vacuum Manifold  -- CW rotation, primordial vacuum base layer
  Section B  UA  Vacuum Manifold  -- CCW rotation, 4-layer UA superstructure
  Section C  DPM Quantum Chain    -- 8-step chain from vacuum to GM/r^2, all 118 elements

DPM FUNDAMENTAL EQUATIONS:
  DPM       = [UA\']/[SCm] = 10              (scale-invariant ratio)
  Grind_opp = omega_CW * SCm - omega_CCW * UA\'
  F_U       = Ug_sum - Ubi + Um

Author: Daniel T. Murphy  |  dpm_vacuum_manifold.py v3.0  |  May 2026
"""

from __future__ import annotations

'''

SCM_SEP = """\
# =============================================================================
# ==================== SECTION A: SCm VACUUM MANIFOLD ========================
# =============================================================================

"""

UA_SEP = """\


# =============================================================================
# ==================== SECTION B: UA VACUUM MANIFOLD =========================
# =============================================================================

"""

DPM_SEP = """\


# =============================================================================
# ==================== SECTION C: DPM QUANTUM CHAIN ASSEMBLY =================
# =============================================================================

"""

consolidated = (
    NEW_HEADER
    + SCM_SEP
    + scm_body
    + UA_SEP
    + ua_body
    + DPM_SEP
    + dpm_body
)

# ---------------------------------------------------------------------------
# 6. Write output
# ---------------------------------------------------------------------------
backup_path = dpm_path + '.bak'
import shutil
shutil.copy2(dpm_path, backup_path)
print(f"Backed up original dpm to: {backup_path}")

write_text(dpm_path, consolidated)
total = len(consolidated.splitlines())
print(f"Written {dpm_path}: {total} lines")
print("MERGE COMPLETE.")
