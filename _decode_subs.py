"""Build clean substitution map by reading the diff file with proper UTF-8."""
import re, unicodedata
from collections import Counter
text = open('_diff020.txt', encoding='utf-8').read()
# Parse pairs
pairs = re.findall(r"PRE: '(.*?)'\s+POST:'(.*?)'", text)
print(f'Pairs found: {len(pairs)}')
# Group by post containing single ?
sub_map = Counter()
for pre, post in pairs:
    # only consider when pre and post are short
    if len(pre) <= 3 and '?' in post:
        sub_map[(pre, post)] += 1

# Print unique substitutions sorted by frequency
print('\n=== UNIQUE SUBSTITUTIONS (pre -> post) ===')
for (pre, post), cnt in sub_map.most_common():
    pre_names = ' / '.join(f'U+{ord(c):04X} {unicodedata.name(c,"?")}' for c in pre)
    post_names = ' / '.join(f'U+{ord(c):04X} {unicodedata.name(c,"?")}' if c != '?' else 'U+003F QUESTION MARK' for c in post)
    print(f'  {cnt:3}  PRE={pre!r:12}  ->  POST={post!r:12}    [{pre_names}] -> [{post_names}]')
