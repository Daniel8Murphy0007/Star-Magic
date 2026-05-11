"""Fix \alphaalpha, \thetatheta, \phiphi, \mumu, \nunu etc → \alpha\alpha"""
import pathlib, re
ROOT = pathlib.Path('whitepapers')
GREEK = ['alpha','beta','gamma','delta','epsilon','zeta','eta','theta','iota','kappa','lambda','mu','nu','xi','rho','sigma','tau','upsilon','phi','chi','psi','omega']
n_total = 0
for p in ROOT.glob('PAPER_*.md'):
    s = p.read_text(encoding='utf-8')
    orig = s
    for g in GREEK:
        s = re.sub(r'\\' + g + g + r'\b', r'\\' + g + r'\\' + g, s)
    if s != orig:
        p.write_text(s, encoding='utf-8')
        n_total += 1
        print(f'{p.name}')
print(f'TOTAL: {n_total}')
