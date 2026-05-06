with open('The Quest for Unity.txt', encoding='utf-8', errors='replace') as f:
    content = f.read()

checks = [
    ('rho_A updated to 1.244e-23',    '1.244\u00d710\u207b\u00b2\u00b3 kg/m\u00b3'),
    ('DPM_RATIO derivation present',  'DPM_RATIO\u00b9\u00b3 / [SSq]'),
    ('P1-P6 block inserted',          'Falsifiable Predictions from UQFF'),
    ('E_crack 700 eV present',        '1.12e-19 J = 700 eV'),
    ('A_26 = 1,307,798,101 present',  '1,307,798,101'),
    ('rho_A dark matter text present','rho_A = rho_SCm * DPM_RATIO^13'),
    ('Old gm/cm3 density gone',       '10\u221223\u2009gm/cm3'),
    ('Old 10^-23 kg/m^3 gone',        '\u301610\u3017^(-23)'),
]
all_pass = True
for name, pattern in checks:
    found = pattern in content
    if name.startswith('Old'):
        ok = not found
        status = 'PASS' if ok else 'FAIL'
    else:
        ok = found
        status = 'PASS' if ok else 'FAIL'
    if not ok:
        all_pass = False
    print(f'  {status}  {name}')

# Count occurrences of the new density string
n = content.count('1.244\u00d710\u207b\u00b2\u00b3 kg/m\u00b3')
print(f'\n  Occurrences of new rho_A value: {n}')
print('\nALL PASS' if all_pass else '\nSOME FAILURES')
