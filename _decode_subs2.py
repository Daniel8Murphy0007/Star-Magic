"""Build clean substitution map by computing the diff in pure Python and
writing UTF-8 bytes directly to file (no PowerShell pipe corruption)."""
import subprocess, re, unicodedata
from difflib import SequenceMatcher
from collections import Counter

PRE_SHA = '107906c7'
POST_SHA = 'd2f9bed6'
PATH = 'whitepapers/PAPER_020_Cosmic_Ray_Propagation_UQFF_Spacetime.md'

pre = subprocess.run(['git','show',f'{PRE_SHA}:{PATH}'], capture_output=True).stdout.decode('utf-8', errors='replace')
post = subprocess.run(['git','show',f'{POST_SHA}:{PATH}'], capture_output=True).stdout.decode('utf-8', errors='replace')

print(f'PRE: {len(pre)} chars, ?-count={pre.count("?")}')
print(f'POST: {len(post)} chars, ?-count={post.count("?")}')

sm = SequenceMatcher(None, pre, post, autojunk=False)
ops = sm.get_opcodes()

# Look for replacement ops with unicode -> ?
single_subs = Counter()  # (pre_char) -> count when pre is single non-ASCII and post is '?'
multi_subs = Counter()
for tag, i1, i2, j1, j2 in ops:
    if tag != 'replace':
        continue
    a = pre[i1:i2]
    b = post[j1:j2]
    if '?' not in b:
        continue
    if len(a) == 1 and b == '?':
        single_subs[a] += 1
    else:
        multi_subs[(a, b)] += 1

print('\n=== SINGLE-CHAR UNICODE -> ? SUBSTITUTIONS ===')
for c, cnt in single_subs.most_common():
    print(f'  {cnt:3}  U+{ord(c):04X}  {c!r}   {unicodedata.name(c,"?")}')

print('\n=== MULTI-CHAR SUBSTITUTIONS (pre -> post containing ?) ===')
for (a, b), cnt in multi_subs.most_common(40):
    a_names = ' / '.join(f'U+{ord(c):04X} {unicodedata.name(c,"?")}' for c in a)
    b_repr = b.replace('\n','\\n')
    print(f'  {cnt:3}  PRE={a!r:8}  POST={b_repr!r:8}    [{a_names}]')
