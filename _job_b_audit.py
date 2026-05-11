import csv
from collections import Counter

rows = list(csv.DictReader(open('_job_b_categorization.csv', encoding='utf-8')))

for b in 'ABCDEFGIJZ':
    sub = [r for r in rows if r['bucket'] == b]
    print('=' * 70)
    print(f'BUCKET {b}  n={len(sub)}')
    kw = Counter(r['matched_keyword'] for r in sub)
    print('  top keywords:', kw.most_common(15))
    print('  samples:')
    for r in sub[:10]:
        print(f"    [{r['matched_keyword']:<30}] {r['filename']}")
