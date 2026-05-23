#!/usr/bin/env python3
"""Fix SCm_Rossi_ECat_Variants_Unified.md table"""

fpath = "whitepapers/SCm_Rossi_ECat_Variants_Unified.md"

with open(fpath, 'r', encoding='utf-8') as f:
    content = f.read()

# Replace the LaTeX table with Markdown table
old_table = r"""\begin{center}
\begin{tabular}{@{}lllll@{}}
\toprule
Variant & $T$ & Trigger & COP & SCm Mode \midrule
Early E-Cat (2011) & 200--400degree C & Electrolytic & 6--14 & Phonon + $F_{U_{Bi_i}}$ E-Cat X (2015) & $\sim$1400degree C & Gas heating & $>$20 & Enhanced $\Phi_T$ + EMF E-Cat SK/SK+ & $\sim$2600degree C & Plasma spark & $>$50 & $t_n$ coherence Star-Magic & ambient & Cold spark (pH $-37$) & 555:1 & Full $F_{U_{Bi_i}}$ + VDS \bottomrule
\end{tabular}
\end{center}"""

new_table = """| Variant | T | Trigger | COP | SCm Mode |
|---------|---|---------|-----|----------|
| Early E-Cat (2011) | 200--400°C | Electrolytic | 6--14 | Phonon + $F_{U_{Bi_i}}$ |
| E-Cat X (2015) | ~1400°C | Gas heating | >20 | Enhanced $\\Phi_T$ + EMF |
| E-Cat SK/SK+ | ~2600°C | Plasma spark | >50 | $t_n$ coherence |
| Star-Magic | ambient | Cold spark (pH -37) | 555:1 | Full $F_{U_{Bi_i}}$ + VDS |"""

if old_table in content:
    content = content.replace(old_table, new_table)
    print("[FIXED] Table replaced")
else:
    print("[ERROR] Table pattern not found, trying regex...")
    import re
    pattern = r'\\begin\{center\}.*?\\end\{center\}'
    if re.search(pattern, content, re.DOTALL):
        content = re.sub(pattern, new_table, content, flags=re.DOTALL)
        print("[FIXED] Table replaced via regex")
    else:
        print("[ERROR] Could not find table to replace")

with open(fpath, 'w', encoding='utf-8') as f:
    f.write(content)

print(f"[SAVED] {fpath}")
