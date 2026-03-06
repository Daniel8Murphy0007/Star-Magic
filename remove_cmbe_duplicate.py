#!/usr/bin/env python3
"""
Remove the dead duplicate compute_master_buoyant_equation at line 93532
from TimeVaryingVacuumModel in CondensedPhysics.py.

Safety checks performed before writing:
  1. Confirms the target line contains the expected def signature.
  2. Confirms the block is inside TimeVaryingVacuumModel.
  3. Confirms the surviving definition at ~line 94141 is still present.
  4. Dry-run mode by default; pass --apply to write.
"""
import re, sys

DRY_RUN = '--apply' not in sys.argv
TARGET_LINE = 93532 - 1  # 0-based index

with open('CondensedPhysics.py', encoding='utf-8', errors='replace') as f:
    lines = f.readlines()

# ── Safety check 1: target line has the right def ──────────────────────────
assert re.match(r'\s+def compute_master_buoyant_equation', lines[TARGET_LINE]), \
    f'ERROR: line {TARGET_LINE+1} does not match expected def. Got: {lines[TARGET_LINE]!r}'

# ── Safety check 2: enclosing class is TimeVaryingVacuumModel ─────────────
enclosing = 'UNKNOWN'
for j in range(TARGET_LINE - 1, -1, -1):
    m = re.match(r'^class (\w+)', lines[j])
    if m:
        enclosing = m.group(1)
        break
assert enclosing == 'TimeVaryingVacuumModel', \
    f'ERROR: enclosing class is {enclosing!r}, expected TimeVaryingVacuumModel'

# ── Determine block extent (def + body until next def/class at same indent) ─
indent = len(lines[TARGET_LINE]) - len(lines[TARGET_LINE].lstrip())
end = TARGET_LINE + 1
while end < len(lines):
    l = lines[end]
    if l.strip() == '':
        end += 1
        continue
    li = len(l) - len(l.lstrip())
    if li <= indent and re.match(r'\s*(def |class )', l):
        break
    end += 1
# Also eat the blank separator line after the block if present
if end < len(lines) and lines[end].strip() == '':
    end += 1

print(f'Block to remove: lines {TARGET_LINE+1}–{end} (0-based {TARGET_LINE}–{end-1})')
print(f'Enclosing class: {enclosing}')
print('Block content:')
for ln in lines[TARGET_LINE:end]:
    print('  ' + ln, end='')
print()

# ── Safety check 3: surviving definition still present ─────────────────────
tail = ''.join(lines[end:])
assert 'def compute_master_buoyant_equation' in tail, \
    'ERROR: no surviving definition found after removal point!'
print('Surviving definition confirmed present.')

# ── Apply (or dry run) ──────────────────────────────────────────────────────
if DRY_RUN:
    print('\nDRY RUN — no changes written. Re-run with --apply to apply.')
else:
    # Remove ONLY the dead block; touch nothing else in the file.
    new_lines = lines[:TARGET_LINE] + lines[end:]
    # Detect original line endings (CRLF vs LF) to preserve them exactly.
    with open('CondensedPhysics.py', 'rb') as fb:
        raw = fb.read(4096)
    newline_seq = b'\r\n' if b'\r\n' in raw else b'\n'
    # Write back preserving original encoding and line endings.
    with open('CondensedPhysics.py', 'wb') as fb:
        for l in new_lines:
            # lines were read in text mode so they end with \n; re-encode with original endings
            lb = l.rstrip('\r\n').encode('utf-8', errors='replace')
            fb.write(lb + newline_seq)
    removed = len(lines) - len(new_lines)
    print(f'\nAPPLIED: removed {removed} lines. Total lines: {len(new_lines)}')

    # Post-check: confirm exactly one definition remains in TimeVaryingVacuumModel
    with open('CondensedPhysics.py', encoding='utf-8', errors='replace') as f:
        src = f.read()
    count = src.count('def compute_master_buoyant_equation')
    print(f'Remaining compute_master_buoyant_equation definitions: {count} (was 115, expected 114)')
