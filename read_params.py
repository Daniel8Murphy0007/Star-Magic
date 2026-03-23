"""Read parameters from extracted module .cpp files."""
import re, os

SYSTEM_MODULES = [
    'Abell2256UQFFModule', 'ASASSN14liUQFFModule', 'CentaurusAUQFFModule',
    'CrabNebulaUQFFModule', 'ElGordoUQFFModule', 'ESO137UQFFModule',
    'IC2163UQFFModule', 'J1610UQFFModule', 'JupiterAuroraeUQFFModule',
    'LagoonNebulaUQFFModule', 'M87JetUQFFModule', 'NGC1365UQFFModule',
    'NGC2207UQFFModule', 'RAquariiUQFFModule', 'SgrAStarUQFFModule',
    'SPTCLJ2215UQFFModule', 'StephanQuintetUQFFModule', 'VelaPulsarUQFFModule',
    'HydrogenResonanceUQFFModule', 'StarMagicUQFFModule', 'SMBHUQFFModule',
    'AstroSystemsUQFFModule', 'UQFFNebulaTriadicModule',
    'UQFFBuoyancyModule', 'UQFF8AstroSystemsModule', 'UQFF8AstroTriadicModule'
]

pat = re.compile(r'variables\["(\w+)"\]\s*=\s*\{([^}]+)\}')

print(f'{"Module":<35} {"M":>12} {"r":>12} {"omega0":>10} {"L_X":>12} {"t":>12} {"E_cm":>14}')
print('-'*100)
for mod in SYSTEM_MODULES:
    fname = mod + '.cpp'
    if not os.path.exists(fname):
        print(f'  MISSING: {fname}')
        continue
    with open(fname, encoding='utf-8', errors='replace') as f:
        text = f.read()
    params = {}
    for m in pat.finditer(text):
        key = m.group(1)
        val = m.group(2).split(',')[0].strip()
        if key not in params:
            params[key] = val
    M = params.get('M', '?')
    r = params.get('r', '?')
    om0 = params.get('omega0', '?')
    LX = params.get('L_X', '?')
    t = params.get('t', '?')
    Ecm = params.get('E_cm', '?')
    print(f'{mod:<35} {M:>12} {r:>12} {om0:>10} {LX:>12} {t:>12} {Ecm:>14}')

print('\n\n=== Unique comments per module ===')
for mod in SYSTEM_MODULES:
    fname = mod + '.cpp'
    if not os.path.exists(fname):
        continue
    with open(fname, encoding='utf-8', errors='replace') as f:
        text = f.read()
    # Show first comment block (watermark/params description)
    lines = text.splitlines()[:8]
    print(f'\n[{mod}]')
    for l in lines:
        if l.strip():
            print(f'  {l[:120]}')
