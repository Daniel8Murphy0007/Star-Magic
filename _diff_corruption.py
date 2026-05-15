"""Compare pre/post 32e896fa to extract every (unicode_char -> ?) substitution
that the upgrade_latex.py script made for PAPER_020 as test case."""
import re
pre = open('_pre020.md', encoding='utf-8').read()
post = open('_post020.md', encoding='utf-8').read()

# Find every '?' in post and look at corresponding region in pre.
# Use difflib SequenceMatcher to align them.
from difflib import SequenceMatcher
sm = SequenceMatcher(None, pre, post, autojunk=False)
ops = sm.get_opcodes()

# Collect all (pre_substring, post_substring) for replace ops where post contains ?
subs = []
for tag, i1, i2, j1, j2 in ops:
    if tag == 'replace':
        a = pre[i1:i2]
        b = post[j1:j2]
        if '?' in b:
            subs.append((a, b))

# For each substitution print
print(f'TOTAL replacement ops with ?: {len(subs)}')
from collections import Counter
patterns = Counter()
for a, b in subs[:50]:
    print(f'  PRE: {a!r}')
    print(f'  POST:{b!r}')
    print()
