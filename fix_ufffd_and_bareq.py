"""
fix_ufffd_and_bareq.py — Pattern-based restoration of corrupted characters.

UFFFD (U+FFFD) replacements  (confirmed high-confidence patterns only):
  1. digit[UFFFD]10 → digit×10  (multiplication before power of 10)
  2. [UFFFD]10       (space/start of number context)
  3. [UFFFD] = 5.0   → κ = 5.0    (UQFF kappa constant, extremely common)
  4. [UFFFD] = 0.0005 → κ = 0.0005
  5. [UFFFD]_i  /  [UFFFD]_SCm → κ_i / κ_SCm  (UQFF sub-constant names)
  6. [UFFFD] = [SSq]  → κ = [SSq]   (sometimes κ appears right before [SSq])
  7. ([UFFFD])  → (κ)  in UQFF constant context
  8. day[UFFFD][UFFFD] → day⁻¹   (common inverse-day unit)
  9. digit[UFFFD]space in unit cells: s⁻¹ / m⁻¹
 10. Space-[UFFFD]-capital → space–capital  (en-dash between proper nouns)
 11. digit[UFFFD]space-digit → digit × digit (multiplication)

Bare-? (ASCII 63) replacements (only in unambiguous math contexts):
  1. ×10?(\d)  →  ×10⁻\d     (scientific notation: ×10⁻N)
  2. ×10?(\d{2}) → ×10⁻\d\d  (two-digit exponents: ×10⁻NM)
  3. 10?(\d) space or | → 10⁻\d (plain 10^-N in tables)
  4. day?[UFFFD or ¹] → day⁻¹  (unit per day, UQFF constant denominator)
  5. s?¹  /  s?[UFFFD] → s⁻¹   (per second)
"""
import sys, re, glob, os, collections
sys.stdout.reconfigure(encoding='utf-8', errors='replace')

U = '\ufffd'

# ── UFFFD substitution rules (ordered by specificity) ─────────────────────
UFFFD_RULES = [
    # 1. digit / decimal ×10 pattern: "2.269[U]10" → "2.269×10"
    (re.compile(r'(\d[\d.]*)\ufffd(10)'), r'\1×\2'),
    # 2. " [U]10" (space before) → " ×10"
    (re.compile(r' \ufffd(10)'), r' ×\1'),
    # 3. "10[U]" might be 10× — handle "10[U]0." pattern (10×0.xxx)
    (re.compile(r'(10)\ufffd(0\.\d)'), r'\1×\2'),
    # 4. Kappa constant: "[U] = 5.0" or "[U] = 0.0005"
    (re.compile(r'\ufffd( = 5\.0)'), r'κ\1'),
    (re.compile(r'\ufffd( = 0\.0005)'), r'κ\1'),
    # 5. UQFF kappa sub-terms: "[U]_i", "[U]_SCm", "[U]_g", "[U]_0"
    (re.compile(r'\ufffd(_[A-Za-z\d_]+)'), r'κ\1'),
    # 6. Parenthesised: "(kappa = ..." from "([U] = 5.0"
    (re.compile(r'\((\ufffd) ='), r'(κ ='),
    # 7. day[U][U] or day[U]¹ → day⁻¹
    (re.compile(r'day\ufffd[\ufffd¹]'), 'day⁻¹'),
    (re.compile(r'day\ufffd\b'), 'day⁻'),
    # 8. " [U] " (isolated) between two word characters → em dash or ×
    #    if both sides are numbers → ×,  otherwise → –
    (re.compile(r'(\d) \ufffd (\d)'), r'\1 × \2'),
    (re.compile(r'([A-Za-z]) \ufffd ([A-Z])'), r'\1 – \2'),
    # 9. "[U]t_H" style expansion coupling: "[U]t" → "·t" (dot product)
    (re.compile(r'\ufffd(t_)'), r'·\1'),
    # 10. Ratio pattern "\ufffd 0." → "≈ 0." (near unity resonance)
    (re.compile(r'\ufffd (0\.\d)'), r'≈ \1'),
    # 11. Table "| [U] |" → "| – |"
    (re.compile(r'\| \ufffd \|'), '| – |'),
    # 12. End-of-line UFFFD after unit letter: "s[U]" or "m[U]" → "s⁻¹" if followed by space/end
    (re.compile(r'\b(s|m|J|K|T|W|kg|yr|Gyr)\ufffd\ufffd(\s|$|\|)'), r'\1⁻¹\2'),
    # 13. " [U][U] " (double) generic → ⁻¹ (most common double-UFFFD usage)
    # Only apply when surrounded by spaces (i.e., standalone)
    (re.compile(r'(?<=\s)\ufffd\ufffd(?=\s)'), '⁻¹'),
]

# ── Bare-? substitution rules (ASCII 63 in unambiguous math contexts) ──────
# These are applied ONLY inside patterns where ? = superscript minus is certain.
BAREQ_RULES = [
    # 1. ×10?(\d{1,2}) → ×10⁻\d\d  (scientific notation)
    (re.compile(r'(×10)\?(\d{1,2})'), lambda m: m.group(1) + '⁻' + m.group(2)),
    # 2. Plain " 10?(\d)" in table cells or space context
    (re.compile(r'(\s|^)(10)\?(\d{1,2})(\s|\||$)'), lambda m: m.group(1) + m.group(2) + '⁻' + m.group(3) + m.group(4)),
    # 3. "day?[UFFFD or ¹]" → "day⁻¹"  (per day)
    (re.compile(r'day\?[\ufffd¹]'), 'day⁻¹'),
    (re.compile(r'day\?\b'), 'day⁻'),
    # 4. "s?¹" or "s?[UFFFD]" → "s⁻¹"  (per second)
    (re.compile(r'\bs\?[¹\ufffd]'), 's⁻¹'),
    # 5. common "×10?4¹" or "×10?4²" patterns (already in ×10 rules above,
    #    but also handle "10?4¹" form)
    (re.compile(r'(10)\?(\d[¹²³⁴⁵⁶⁷⁸⁹⁰])'), lambda m: m.group(1) + '⁻' + m.group(2)),
    # 6. Multi-digit bare-? exponents that are still bare ASCII digits
    #    "10?47" → "10⁻⁴⁷"    (if 2 digits after ? both are exponent)
    # Only apply in table/equation context (pipe-delimited or space-terminated)
    (re.compile(r'(?<=[×\s])(10)\?(\d)(\d)(?=[\s|*\)])'), lambda m: m.group(1) + '⁻' + m.group(2) + m.group(3)),
]

SUPERSCRIPT_MAP = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')

def apply_rules(text, rules):
    for pat, repl in rules:
        if callable(repl):
            text = pat.sub(repl, text)
        else:
            text = pat.sub(repl, text)
    return text

changed_files = []
total_ufffd_fixed = 0
total_bareq_fixed = 0

for fpath in sorted(glob.glob('whitepapers/PAPER_*.md')):
    with open(fpath, encoding='utf-8', errors='replace') as fh:
        original = fh.read()

    text = original

    # Count UFFFD and bare-? before
    ufffd_before = text.count(U)
    # Bare ? in math context (rough count)
    bareq_before = len(re.findall(r'[×0-9]\?[\d\ufffd¹]|day\?|s\?[¹\ufffd]', text))

    # Apply UFFFD fixes
    text = apply_rules(text, UFFFD_RULES)

    # Apply bare-? fixes
    text = apply_rules(text, BAREQ_RULES)

    ufffd_after = text.count(U)
    bareq_after = len(re.findall(r'[×0-9]\?[\d\ufffd¹]|day\?|s\?[¹\ufffd]', text))

    ufffd_fixed = ufffd_before - ufffd_after
    bareq_fixed = bareq_before - bareq_after

    if ufffd_fixed > 0 or bareq_fixed > 0:
        total_ufffd_fixed += ufffd_fixed
        total_bareq_fixed += bareq_fixed
        changed_files.append((os.path.basename(fpath), ufffd_fixed, bareq_fixed))
        with open(fpath, 'w', encoding='utf-8') as fh:
            fh.write(text)

print(f'Files modified  : {len(changed_files)}')
print(f'UFFFD fixed     : {total_ufffd_fixed}')
print(f'Bare-? fixed    : {total_bareq_fixed}')
print()
print('Top 20 most-changed files:')
changed_files.sort(key=lambda x: x[1]+x[2], reverse=True)
for fname, u, q in changed_files[:20]:
    print(f'  {fname}: UFFFD={u} bareq={q}')

# Remaining UFFFD summary
remaining = 0
for fpath in glob.glob('whitepapers/PAPER_*.md'):
    with open(fpath, encoding='utf-8', errors='replace') as fh:
        remaining += fh.read().count(U)
print(f'\nRemaining UFFFD in corpus: {remaining}')
