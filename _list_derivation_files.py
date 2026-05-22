"""List ALL files with 'derivation', 'derive', 'derivative', 'audit', or
'closure' in the filename. Mark each one as either CALLED or NOT CALLED by
_uqff_program.py --audit."""
import sys, re
from pathlib import Path
sys.path.insert(0, '.')
from _uqff_program import ROOT, SESSION_RE, EXTRA_DERIVATION_FILES

# What _uqff_program.py calls
called = set()
for p in ROOT.glob('_session*.py'):
    called.add(p.name)
for f in EXTRA_DERIVATION_FILES:
    called.add(f)

KEYWORDS = re.compile(r'(derivation|derive|derivative|audit|closure)', re.I)

matches = []
for p in sorted(ROOT.rglob('*.py')):
    # Skip virtualenvs and node_modules
    s = str(p).lower()
    if '.venv' in s or 'node_modules' in s or 'site-packages' in s:
        continue
    if KEYWORDS.search(p.name):
        rel = p.relative_to(ROOT)
        matches.append((p.name, str(rel), p.stat().st_size))

print(f"Total files matching derivation/derive/audit/closure keywords: {len(matches)}\n")

called_list = []
notcalled_list = []
for name, rel, size in matches:
    if name in called:
        called_list.append((name, rel, size))
    else:
        notcalled_list.append((name, rel, size))

print(f"=== CALLED BY _uqff_program.py ({len(called_list)}) ===")
for name, rel, size in called_list:
    print(f"  {size:>10,} bytes  {rel}")

print(f"\n=== NOT CALLED BY _uqff_program.py ({len(notcalled_list)}) ===")
for name, rel, size in notcalled_list:
    print(f"  {size:>10,} bytes  {rel}")
