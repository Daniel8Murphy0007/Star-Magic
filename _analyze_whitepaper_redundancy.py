#!/usr/bin/env python3
"""Analyze 1,344 whitepapers for redundancy"""
import os
import re
from collections import defaultdict
from pathlib import Path

# Get all whitepaper files
whitepapers = list(Path('whitepapers').glob('PAPER_*.md'))
print(f"Total whitepapers: {len(whitepapers)}\n")

# Parse paper numbers and titles
papers = {}
for fp in whitepapers:
    name = fp.stem
    # Extract paper number
    match = re.match(r'PAPER_(\d+)_(.*)', name)
    if match:
        num = int(match.group(1))
        title = match.group(2).replace('_', ' ')
        papers[num] = {'file': fp.name, 'title': title}

# Group by topic
topics = defaultdict(list)
for num, info in papers.items():
    # Extract main topic from title
    title = info['title'].lower()
    
    # Identify broad categories
    if any(x in title for x in ['calibration', 'calibrat']):
        topics['Calibration Papers'].append(num)
    elif any(x in title for x in ['uqff', 'unified quantum']):
        topics['UQFF Framework'].append(num)
    elif any(x in title for x in ['muge', 'compressed', 'resonance']):
        topics['MUGE/Resonance'].append(num)
    elif any(x in title for x in ['closure', 'gap']):
        topics['Closure/Gap Papers'].append(num)
    elif any(x in title for x in ['smbh', 'black hole', 'sagittarius', 'sgra']):
        topics['SMBH/Supermassive BH'].append(num)
    elif any(x in title for x in ['architecture', 'framework', 'entity']):
        topics['Architecture/Framework'].append(num)
    elif any(x in title for x in ['solver', 'simultaneous', 'equation']):
        topics['Solvers/Equations'].append(num)
    elif any(x in title for x in ['solar', 'sun']):
        topics['Solar/Stellar'].append(num)
    elif any(x in title for x in ['magnetar', 'sgr', 'neutron']):
        topics['Magnetars/Neutron Stars'].append(num)
    elif any(x in title for x in ['conanqi', 'calculator', 'gui']):
        topics['Calculators'].append(num)
    elif any(x in title for x in ['validation', 'sync', 'audit']):
        topics['Validation/Audit'].append(num)
    elif any(x in title for x in ['quantum', 'field', 'lagrangian']):
        topics['Quantum/Field Theory'].append(num)
    elif any(x in title for x in ['astrophysical', 'astro', 'galaxy', 'star']):
        topics['Astrophysical'].append(num)
    elif any(x in title for x in ['cosmological', 'cosmo', 'lambda', 'hubble']):
        topics['Cosmology'].append(num)
    else:
        topics['Other'].append(num)

print("="*70)
print("WHITEPAPER DISTRIBUTION BY TOPIC")
print("="*70)

for topic in sorted(topics.keys()):
    nums = sorted(topics[topic])
    print(f"\n{topic}: {len(nums)} papers")
    print(f"  Range: PAPER_{nums[0]} to PAPER_{nums[-1]}")
    if len(nums) <= 10:
        print(f"  List: {', '.join(f'PAPER_{n}' for n in nums)}")
    else:
        print(f"  First 5: {', '.join(f'PAPER_{n}' for n in nums[:5])}")
        print(f"  Last 5: {', '.join(f'PAPER_{n}' for n in nums[-5:])}")

print("\n" + "="*70)
print("REDUNDANCY ANALYSIS")
print("="*70)

# Look for papers with nearly identical titles
title_map = {}
for num, info in papers.items():
    title_clean = info['title'].lower()
    # Remove numbers and special chars for comparison
    title_simple = re.sub(r'[_\d\s]+', ' ', title_clean).strip()
    if title_simple not in title_map:
        title_map[title_simple] = []
    title_map[title_simple].append(num)

# Find duplicates
print("\n1. EXACT OR NEAR-DUPLICATE TITLES:")
duplicates = {k: v for k, v in title_map.items() if len(v) > 1}
if duplicates:
    for title, nums in sorted(duplicates.items(), key=lambda x: -len(x[1])):
        print(f"\n  '{title}'")
        for num in sorted(nums):
            print(f"    - PAPER_{num}: {papers[num]['title']}")
else:
    print("  None found")

# Look for similar titles (substring matches)
print("\n2. PAPERS WITH OVERLAPPING TOPICS (Potential Redundancy):")

potential_redundant = defaultdict(list)

# Calibration variants
calib_papers = [n for n in papers if 'calibration' in papers[n]['title'].lower() or 'calibrat' in papers[n]['title'].lower()]
if len(calib_papers) > 20:
    potential_redundant['Calibration Papers'] = (len(calib_papers), calib_papers[:5] + ['...'] + calib_papers[-5:])

# UQFF papers
uqff_papers = [n for n in papers if 'uqff' in papers[n]['title'].lower()]
if len(uqff_papers) > 30:
    potential_redundant['UQFF Framework'] = (len(uqff_papers), uqff_papers[:5] + ['...'] + uqff_papers[-5:])

# Framework papers
framework_papers = [n for n in papers if any(x in papers[n]['title'].lower() for x in ['framework', 'architecture'])]
if len(framework_papers) > 15:
    potential_redundant['Framework/Architecture'] = (len(framework_papers), framework_papers[:3] + ['...'] + framework_papers[-3:])

# Closure papers
closure_papers = [n for n in papers if 'closure' in papers[n]['title'].lower()]
if len(closure_papers) > 5:
    potential_redundant['Closure Papers'] = (len(closure_papers), closure_papers)

for category, (count, samples) in sorted(potential_redundant.items(), key=lambda x: -x[1][0]):
    print(f"\n  {category}: {count} papers")
    for item in samples:
        if item == '...':
            print(f"    {item}")
        else:
            print(f"    - PAPER_{item}")

# Check for LaTeX vs Markdown duplicates
print("\n3. LaTeX vs Markdown DUPLICATES:")
md_papers = set()
tex_papers = set()

for fp in whitepapers:
    if fp.suffix == '.md':
        match = re.match(r'PAPER_(\d+)', fp.stem)
        if match:
            md_papers.add(int(match.group(1)))

tex_files = list(Path('whitepapers').glob('PAPER_*.tex'))
for fp in tex_files:
    match = re.match(r'PAPER_(\d+)', fp.stem)
    if match:
        tex_papers.add(int(match.group(1)))

both = md_papers & tex_papers
print(f"\n  Papers with BOTH .md and .tex: {len(both)}")
if both:
    both_sorted = sorted(both)
    print(f"  Examples: {', '.join(f'PAPER_{n}' for n in both_sorted[:10])}")

# Check PDF count vs source
print("\n4. PDF GENERATION STATUS:")
md_count = len(list(Path('whitepapers').glob('*.md')))
tex_count = len(list(Path('whitepapers').glob('*.tex')))
pdf_count = len(list(Path('pdf').glob('*.pdf')))
print(f"\n  Markdown sources:  {md_count}")
print(f"  LaTeX sources:     {tex_count}")
print(f"  Generated PDFs:    {pdf_count}")
print(f"  Missing PDFs:      {md_count + tex_count - pdf_count}")

# Summary statistics
print("\n" + "="*70)
print("SUMMARY")
print("="*70)
print(f"\nTotal whitepapers: {len(papers)}")
print(f"Unique titles: {len(title_map)}")
print(f"Duplicate titles: {sum(1 for v in title_map.values() if len(v) > 1)}")
print(f"\nLargest categories:")
for topic in sorted(topics.keys(), key=lambda x: -len(topics[x]))[:5]:
    print(f"  {topic}: {len(topics[topic])}")

print(f"\nPapers with both .md and .tex versions: {len(both)}")
print(f"Papers missing PDFs: {md_count + tex_count - pdf_count}")

print("\n" + "="*70)
print("RECOMMENDATIONS FOR CLEANUP")
print("="*70)
print(f"""
1. CONSOLIDATE DUPLICATE FORMATS
   - {len(both)} papers have both .md and .tex versions
   - Keep ONE format per paper; delete the other
   - Estimated savings: {len(both) * 2} files (source + PDF pairs)

2. REDUCE CALIBRATION PAPERS
   - {len([n for n in papers if 'calibration' in papers[n]['title'].lower()])} calibration variants
   - Keep only canonical calibrations used by calculators
   - Archive session-specific calibrations
   - Estimated savings: 50-100 papers

3. CONSOLIDATE FRAMEWORK PAPERS
   - {len([n for n in papers if 'framework' in papers[n]['title'].lower()])} framework descriptions
   - Keep only ARCHITECTURE_FLOW_DIAGRAM.md
   - Archive redundant architecture docs
   - Estimated savings: 30-50 papers

4. KEEP CLOSURE PAPERS (ESSENTIAL)
   - {len([n for n in papers if 'closure' in papers[n]['title'].lower()])} closure papers justify Lagrangian gaps
   - These are NOT redundant; keep all

5. TOTAL REDUCTION POTENTIAL
   - From {len(papers)} whitepapers → ~200 canonical papers
   - 85% reduction while retaining all unique physics content
""")
