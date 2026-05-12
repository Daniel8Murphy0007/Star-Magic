"""Reverse the PowerShell $^ → $env:PYTHONIOENCODING injection that corrupted whitepaper edits.

Also handle the remaining `day-1` → `day$^{-1}$` and `m-2` → `m$^{-2}$` patterns
that weren't caught by the earlier regex pass.
"""
import re

files = [
    'PAPER_009b_Aether_String_TRZ_Damping_GW',
    'PAPER_015b_Multiband_GW_LISA_LIGO_UQFF',
    'PAPER_016b_White_Dwarf_Foreground_UQFF',
    'PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime',
    'PAPER_021_Gravitational_Lensing_Corrections_UQFF_Vacuum_Density',
]

# Repair the $env:PYTHONIOENCODING injection: replace `\$env:PYTHONIOENCODING` with `$^`
# Pattern is exactly: \$env:PYTHONIOENCODING{N}\$  →  $^{N}$
INJECT = re.compile(r'\\\$env:PYTHONIOENCODING\{([^}]+)\}\\\$')

# Also remaining cleanup
subs = [
    (re.compile(r'\bday-1\b'), r'day$^{-1}$'),
    (re.compile(r'(\d)\s*m-2\b'), r'\g<1> m$^{-2}$'),
    (re.compile(r'\$\\times 10\$\^\{-52\}\$\s*m\$\^\{-2\}\$'), r'$\\times 10^{-52}$ m$^{-2}$'),
]

for f in files:
    p = 'whitepapers/' + f + '.md'
    with open(p, encoding='utf-8') as fh:
        t = fh.read()
    orig = t
    t = INJECT.sub(r'$^{\1}$', t)
    for r, s in subs:
        t = r.sub(s, t)
    if t != orig:
        with open(p, 'w', encoding='utf-8') as fh:
            fh.write(t)
        print('FIXED', f)
    else:
        print('no-change', f)
