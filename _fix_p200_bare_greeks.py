import re, pathlib
fp = pathlib.Path('whitepapers/PAPER_200_Um_Universal_Magnetism_Taxonomy_Complete_Variant_Catalogue.md')
text = fp.read_text(encoding='utf-8')
lines = text.split('\n')
new_lines = []
in_code = False
in_block_math = False
GREEK = r'(mu|nabla|sim|Theta|rho|beta|Sigma|Delta|omega|gamma|phi|Phi|kappa|alpha|tau|sigma|epsilon|theta|lambda|Lambda|Omega|pi|chi|psi|Psi|eta|zeta|xi|Xi|approx|times|cdot|cos|sin|tan|exp|sqrt|log|ln|frac|leq|geq|neq|infty|propto|sum|int|partial|hbar|odot)'
for line in lines:
    stripped = line.strip()
    if stripped.startswith('```'):
        in_code = not in_code
        new_lines.append(line); continue
    if in_code:
        new_lines.append(line); continue
    dd = line.count('$$')
    if in_block_math:
        new_lines.append(line)
        if dd % 2 == 1:
            in_block_math = False
        continue
    else:
        if dd >= 2:
            new_lines.append(line); continue
        if dd == 1:
            new_lines.append(line); in_block_math = True; continue
    s = line
    out_chars = []
    j = 0
    inm = False
    while j < len(s):
        c = s[j]
        if c == '$':
            inm = not inm
            out_chars.append(c); j += 1; continue
        if not inm and c == '\\' and j+1 < len(s):
            m = re.match(r'\\' + GREEK + r'\b', s[j:])
            if m:
                token = m.group(0)
                out_chars.append(f'${token}$'); j += len(token); continue
        out_chars.append(c); j += 1
    new_lines.append(''.join(out_chars))
fp.write_text('\n'.join(new_lines), encoding='utf-8')
print('done')
