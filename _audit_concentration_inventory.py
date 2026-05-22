"""Find EVERYTHING - not just what's in master_closures.csv.
Answers: real vs fantasy, what's missing, what target files exist for concentration.
"""
import csv, re, os
from pathlib import Path
from collections import Counter, defaultdict

ROOT = Path('.')

# 1. ALL _session*.py at root
all_sessions = sorted(ROOT.glob('_session*.py'))
all_session_names = {p.name for p in all_sessions}

# 2. Sessions REFERENCED in master_closures.csv
referenced = set()
with open('master_closures.csv', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        s = (row.get('script') or '').strip()
        if s:
            referenced.add(Path(s).name)

# 3. Sessions on disk but NOT in CSV
orphans = sorted(all_session_names - referenced)
# 4. Sessions in CSV but missing on disk
ghosts = sorted(referenced - all_session_names)

# 5. Classify orphans by size + keyword
KEYWORDS = {
    'variational': re.compile(r'variational|sustain|action_principle|euler.lagrange|stationary', re.I),
    'lagrangian':  re.compile(r'lagrangian|hamiltonian|noether|symmetry', re.I),
    'manifold':    re.compile(r'manifold|metric|connection|christoffel|geodesic', re.I),
    'quantum':     re.compile(r'quantum|wavefunction|operator|commutator|schrodinger', re.I),
    'cosmology':   re.compile(r'cosmolog|hubble|friedmann|flrw|big.bang|inflation', re.I),
    'gravity':     re.compile(r'gravity|gravitation|einstein|riemann|curvature|black.hole', re.I),
    'thermo':      re.compile(r'thermo|entropy|temperature|partition|boltzmann', re.I),
    'em':          re.compile(r'maxwell|electromag|gauge|yang.mills', re.I),
    'fluid':       re.compile(r'fluid|navier|stokes|vortex|turbulence', re.I),
    'numeric':     re.compile(r'verify|check|test|trace|debug|inspect', re.I),
}

def classify(p: Path):
    sz = p.stat().st_size
    try:
        text = p.read_text(encoding='utf-8', errors='ignore')[:5000]
    except Exception:
        text = ''
    tags = [k for k, rx in KEYWORDS.items() if rx.search(text) or rx.search(p.name)]
    has_print = 'print(' in text
    has_main = '__main__' in text or 'if __name__' in text
    return sz, tags, has_print, has_main

# 6. Also count OTHER python derivation-flavor files NOT named _session
other_py_patterns = ['*_derivation*.py', '*derive*.py', '*_proof*.py', '*_chain*.py', '*qcg*.py',
                     '*lagrangian*.py', '*variational*.py', '*sustainarity*.py', '*sustainability*.py']
other_files = set()
for pat in other_py_patterns:
    for p in ROOT.glob(pat):
        if not p.name.startswith('_session') and p.is_file():
            other_files.add(p)
other_files = sorted(other_files)

# 7. EXISTING concentration target files
TARGET_CANDIDATES = [
    'CondensedPhysics.py', 'CondensedPhysics2.py', 'CondensedPhysics3.py',
    'CondensedPhysics_OutputData.py',
    'MAIN_1_CoAnQi.cpp', 'index.js', 'MAIN_1.cpp',
    'QCalc.py', 'QCalcGeom.py',
    'source2.cpp', 'physics_backend.cpp',
    'derivations/uqff_all_derivations.py',
]
existing_targets = []
for t in TARGET_CANDIDATES:
    p = ROOT / t
    if p.exists():
        existing_targets.append((t, p.stat().st_size, sum(1 for _ in open(p, encoding='utf-8', errors='ignore'))))

# 8. Look specifically for "variational sustainarity"
variational_hits = []
for p in all_sessions + other_files:
    try:
        txt = p.read_text(encoding='utf-8', errors='ignore')
    except Exception:
        continue
    if re.search(r'variational|sustainar|sustainabilit', txt, re.I):
        variational_hits.append((p.name, p.stat().st_size))

# ============= REPORT =============
print(f"=== TOTAL INVENTORY ===")
print(f"_session*.py on disk:              {len(all_sessions)}")
print(f"_session*.py referenced in CSV:    {len(referenced)}")
print(f"ORPHANS (on disk, not in CSV):     {len(orphans)}")
print(f"GHOSTS  (in CSV, not on disk):     {len(ghosts)}")
print(f"OTHER derivation .py files:        {len(other_files)}")
print(f"GRAND TOTAL derivation-flavor:     {len(all_sessions) + len(other_files)}")
print()

print(f"=== ORPHANED SESSIONS (not in master_closures.csv) ===")
if orphans:
    for n in orphans[:50]:
        sz = (ROOT / n).stat().st_size
        print(f"  {n}  ({sz} B)")
    if len(orphans) > 50:
        print(f"  ... +{len(orphans)-50} more")
else:
    print("  (none - every disk session is in CSV)")
print()

print(f"=== GHOSTS (CSV references file missing) ===")
if ghosts:
    for n in ghosts[:30]:
        print(f"  {n}")
    if len(ghosts) > 30:
        print(f"  ... +{len(ghosts)-30} more")
else:
    print("  (none)")
print()

print(f"=== OTHER DERIVATION FILES (non-_session) ===")
for p in other_files:
    sz = p.stat().st_size
    print(f"  {p.name}  ({sz} B)")
print()

print(f"=== 'VARIATIONAL / SUSTAINARITY' KEYWORD HITS ===")
if variational_hits:
    for n, sz in variational_hits:
        print(f"  {n}  ({sz} B)")
else:
    print("  NONE FOUND - this is fantasy or named differently")
print()

print(f"=== EXISTING CONCENTRATION TARGETS (no new files needed) ===")
for name, sz, lines in existing_targets:
    print(f"  {name:<45} {sz:>12} B  {lines:>8} lines")
print()

# Size tiers of all sessions => real vs stub
tiers = Counter()
substantial = []
for p in all_sessions:
    sz = p.stat().st_size
    if sz < 500:
        tiers['stub_<500B'] += 1
    elif sz < 2000:
        tiers['small_500-2K'] += 1
    elif sz < 10000:
        tiers['medium_2-10K'] += 1
    else:
        tiers['LARGE_>=10K'] += 1
        substantial.append((p.name, sz))

print(f"=== SESSION SIZE TIERS (real vs stub) ===")
for k in ('stub_<500B', 'small_500-2K', 'medium_2-10K', 'LARGE_>=10K'):
    print(f"  {k:<20} {tiers[k]}")
print()
print(f"=== LARGEST 30 SESSIONS (most likely REAL derivations) ===")
substantial.sort(key=lambda x: -x[1])
for n, sz in substantial[:30]:
    print(f"  {n}  ({sz} B)")

# Save full orphan + ghost lists
with open('_audit_concentration_inventory.csv', 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f)
    w.writerow(['file', 'kind', 'size_bytes', 'in_master_csv'])
    for p in all_sessions:
        in_csv = p.name in referenced
        w.writerow([p.name, 'session', p.stat().st_size, in_csv])
    for p in other_files:
        w.writerow([p.name, 'other_derivation', p.stat().st_size, False])
print()
print(f"Full per-file table -> _audit_concentration_inventory.csv")
