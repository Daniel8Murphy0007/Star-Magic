"""Corpus analysis: PAPER_001-900 size, structure, and coverage metrics."""
import os, re

WP_DIR = "whitepapers"
papers = {}
for f in os.listdir(WP_DIR):
    m = re.match(r'PAPER_(\d+).*\.md$', f)
    if m:
        num = int(m.group(1))
        path = os.path.join(WP_DIR, f)
        size = os.path.getsize(path)
        lines = sum(1 for _ in open(path, encoding='utf-8', errors='replace'))
        with open(path, encoding='utf-8', errors='replace') as fh:
            text = fh.read()
        has_abstract = '## Abstract' in text or '## abstract' in text.lower()
        has_secA = '§A.' in text or 'Cosmogenesis-Linked' in text
        has_secB = '§B.' in text or 'VDS/DVP/BSH' in text
        has_SM = '§SM' in text or 'SM Anchors' in text
        has_eqs = '```' in text or '$$' in text
        has_cvw = 'CVW' in text
        has_s204 = 'Session 204' in text
        has_s209 = 'Session 209' in text
        papers[num] = {
            'file': f, 'size': size, 'lines': lines,
            'abstract': has_abstract, 'secA': has_secA, 'secB': has_secB,
            'SM': has_SM, 'eqs': has_eqs, 'cvw': has_cvw,
            's204': has_s204, 's209': has_s209
        }

# Report
total = len(papers)
present = sorted(papers.keys())
missing = [n for n in range(1, 901) if n not in papers]

print(f"=== CORPUS ANALYSIS: PAPER_001-900 ===")
print(f"Total papers found: {total}")
print(f"Missing paper numbers: {len(missing)}")
if missing:
    print(f"  Missing: {missing[:30]}{'...' if len(missing)>30 else ''}")

# Size distribution
sizes = [papers[n]['lines'] for n in present]
stubs = [n for n in present if papers[n]['lines'] < 30]
small = [n for n in present if 30 <= papers[n]['lines'] < 100]
medium = [n for n in present if 100 <= papers[n]['lines'] < 300]
large = [n for n in present if papers[n]['lines'] >= 300]

print(f"\n=== SIZE DISTRIBUTION ===")
print(f"Stubs (<30 lines):   {len(stubs)}")
print(f"Small (30-99 lines): {len(small)}")
print(f"Medium (100-299):    {len(medium)}")
print(f"Large (300+):        {len(large)}")
print(f"Average lines:       {sum(sizes)/len(sizes):.0f}")
print(f"Min/Max:             {min(sizes)}/{max(sizes)}")

# Structure coverage
with_abstract = sum(1 for n in present if papers[n]['abstract'])
with_secA = sum(1 for n in present if papers[n]['secA'])
with_secB = sum(1 for n in present if papers[n]['secB'])
with_SM = sum(1 for n in present if papers[n]['SM'])
with_eqs = sum(1 for n in present if papers[n]['eqs'])
with_cvw = sum(1 for n in present if papers[n]['cvw'])
with_s204 = sum(1 for n in present if papers[n]['s204'])
with_s209 = sum(1 for n in present if papers[n]['s209'])

print(f"\n=== STRUCTURAL COVERAGE ===")
print(f"Has Abstract:        {with_abstract}/{total} ({100*with_abstract/total:.1f}%)")
print(f"Has §A Cosmogenesis: {with_secA}/{total} ({100*with_secA/total:.1f}%)")
print(f"Has §B VDS/DVP/BSH:  {with_secB}/{total} ({100*with_secB/total:.1f}%)")
print(f"Has §SM Anchors:     {with_SM}/{total} ({100*with_SM/total:.1f}%)")
print(f"Has equations:       {with_eqs}/{total} ({100*with_eqs/total:.1f}%)")
print(f"Has CVW ref:         {with_cvw}/{total} ({100*with_cvw/total:.1f}%)")
print(f"Has S204 ref:        {with_s204}/{total} ({100*with_s204/total:.1f}%)")
print(f"Has S209 ref:        {with_s209}/{total} ({100*with_s209/total:.1f}%)")

# Stubs list
if stubs:
    print(f"\n=== STUB PAPERS (<30 lines) ===")
    for n in stubs[:20]:
        print(f"  PAPER_{n:03d}: {papers[n]['lines']}L - {papers[n]['file']}")
    if len(stubs) > 20:
        print(f"  ... and {len(stubs)-20} more")

# Papers without key sections
no_abstract = [n for n in present if not papers[n]['abstract']]
no_eqs = [n for n in present if not papers[n]['eqs']]
print(f"\n=== PAPERS WITHOUT ABSTRACT ===")
print(f"Count: {len(no_abstract)}")
if no_abstract:
    print(f"  First 15: {no_abstract[:15]}")
print(f"\n=== PAPERS WITHOUT EQUATIONS ===")
print(f"Count: {len(no_eqs)}")
if no_eqs:
    print(f"  First 15: {no_eqs[:15]}")

# Session 209 papers
print(f"\n=== SESSION 209 PAPERS ===")
s209_papers = [n for n in present if papers[n]['s209']]
print(f"Count: {len(s209_papers)}")
print(f"Papers: {s209_papers}")
