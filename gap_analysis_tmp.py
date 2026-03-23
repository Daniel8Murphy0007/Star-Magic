import os

wp_dir = 'whitepapers'
papers = sorted(os.listdir(wp_dir))

targets = ['J1610', 'PLCK_G287', 'PLCK G287', 'PSZ2_G181', 'PSZ2 G181', 'ASKAP_J1832', 'ASKAP J1832', 'Sonification', 'Cosmic Neutrino Background', 'CNBModule', 'UQFFBuoyancy']
hits = {t: [] for t in targets}

for paper in papers:
    fpath = os.path.join(wp_dir, paper)
    with open(fpath, encoding='utf-8', errors='replace') as f:
        content = f.read()
    for t in targets:
        if t.lower() in content.lower():
            hits[t].append(paper)

print('=== Gap Analysis ===')
for t, papers_found in hits.items():
    status = papers_found if papers_found else ['NOT FOUND - NEW']
    print(f'{t}: {status}')

print()
print('=== Total papers in whitepapers/ ===')
print(len(papers))

# Show highest paper number
nums = []
for p in papers:
    m = __import__('re').search(r'PAPER_(\d+)', p)
    if m:
        nums.append(int(m.group(1)))
print('Max paper number:', max(nums) if nums else 'N/A')
print('Paper count:', len(nums))

print()
print('=== CNB setSystemParams ===')
with open('UQFFBuoyancyCNBModule.cpp', encoding='utf-8') as f:
    cnb = f.read()
sp = cnb.find('setSystemParams')
print(cnb[sp:sp+3000])
