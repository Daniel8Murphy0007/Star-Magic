"""Repair mojibake in PAPER_172 by character substitution."""
import re
p = 'whitepapers/PAPER_172_FU_Complete_Unified_Field_Assembly.md'
t = open(p,'r',encoding='utf-8').read()

# Common mojibake substitutions (UTF-8 bytes interpreted as cp1252 then re-encoded UTF-8)
subs = {
    'Ã—': '*',          # ×
    'âˆ’': '-',          # −
    'â€"': '-',          # — em-dash
    'â€”': '-',          # — em-dash
    'â€“': '-',          # – en-dash
    'Î©': 'Omega',       # Ω
    'Î\x81': 'rho',     # corrupted ρ
    'Î\\mu': 'mu',
    'Î²': 'beta',
    'Î²^2': 'beta^2',
    'Î¼': 'mu',
    'Î¼_sw': 'mu_sw',
    'Î\u00b41/4': 'mu',  # mu in some forms
    'Î\xb01/4': 'mu',
    'Ï\x81': 'rho',
    'Ï\x80': 'pi',
    'Ï€': 'pi',
    'Ï\x84': 'tau',
    'Ï‰': 'omega',
    'â‚™': '_n',
    'â‚': '_',
    'Â³': '^3',
    'Â²': '^2',
    'Âµ': 'mu',
    'Â°': 'deg',
    'Î\u00b1': 'alpha',
    'Î±': 'alpha',
    'Î²': 'beta',
    'Î³': 'gamma',
    'Î´': 'delta',
    'Îµ': 'epsilon',
    'Î\xa3': 'Sigma',
    'Î£': 'Sigma',
    'Î›': 'Lambda',
    'Î¦': 'Phi',
    'Î¨': 'Psi',
    'Î\xb8': 'theta',
    'Î¸': 'theta',
    'Î\u00ba': 'kappa',
    'Îº': 'kappa',
    'Î»': 'lambda',
    'Î½': 'nu',
    'Î¾': 'xi',
    'Î\u0081': 'rho',  # explicit
    'Î': '',  # leftover - replace with empty (last resort)
    'Ï': '',
    'Ã': '',
    'â€': '',
    'â‚': '',
    'â': '',
}
# Apply longest-first to avoid clobbering prefixes
for k in sorted(subs.keys(), key=len, reverse=True):
    t = t.replace(k, subs[k])

open(p,'w',encoding='utf-8').write(t)
print('done')
