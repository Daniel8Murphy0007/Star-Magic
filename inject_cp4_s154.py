"""inject_cp4_s154.py — inject Session 154 classes into CondensedPhysics4.py"""
import re, sys

# ── 1. Read files ─────────────────────────────────────────────────────────────
with open('CondensedPhysics4.py', encoding='utf-8') as f:
    cp4 = f.read()

with open('session_154_physics_registry.py', encoding='utf-8') as f:
    reg = f.read()

# ── 2. Extract the five class blocks from the registry ────────────────────────
#  Pattern: separator line + comment header + class body, until next separator or if __name__
pat = re.compile(
    r'(# -{74,}\n# #16[1-5].+?\n# -{74,}\n'   # header block
    r'class \w+.*?'                              # class definition
    r'(?=\n# -{74,}\n# #16|\nif __name__))',     # stop at next class header or __main__
    re.DOTALL
)
blocks = pat.findall(reg)
if len(blocks) != 5:
    sys.exit(f'ERROR: expected 5 blocks, got {len(blocks)}')

# ── 3. Build the constants block (module-level constants needed by the classes)
#  These are defined at the top of the registry and must precede the classes in CP4.
const_header = "\n# =============================================================================\n"
const_header += "# Session 154 constants (Big Bang Hypergraph / Nuclear Epoch)\n"
const_header += "# =============================================================================\n"

consts_pat = re.compile(
    r'(_FAC26\s*=.+?_DVP_P\s*=.+?\n)',
    re.DOTALL
)
# Pull the constants section from registry (after imports, before epoch table)
const_match = re.search(
    r'(import math\n\n# -{74,}\n# Constants\n# -{74,}\n)(.*?)(\n# Mayan epoch)',
    reg, re.DOTALL
)
if const_match:
    const_body = const_match.group(2)
else:
    # Fallback: grab everything from _FAC26 through _EPOCHS definition
    const_body = re.search(r'(_FAC26 = math\.factorial.*?\n\n)', reg, re.DOTALL)
    const_body = const_body.group(1) if const_body else ''

# Grab the full constants + helper functions block
# From first _FAC26 line to just before first '# ---' class separator
const_section_match = re.search(
    r'(import math\n\n)(# -{74,}\n# Constants.*?(?=\n# -{74,}\n# #161))',
    reg, re.DOTALL
)
if const_section_match:
    const_section = const_section_match.group(2)
else:
    # Manual extraction: everything from _FAC26 to first class separator
    const_section = re.search(
        r'(_FAC26 = math\.factorial.*?)(?=\n# -{74,}\n# #161)',
        reg, re.DOTALL
    )
    const_section = const_section.group(1) if const_section else ''

# ── 4. Build insertion text ───────────────────────────────────────────────────
insert_text = (
    "\n\n"
    "# =============================================================================\n"
    "# Session 154: Universal Epoch / Periodic Table UQFF\n"
    "# Source: grok_share_efc8a971378f.txt\n"
    "# Classes #161-165  PAPER_573 / PAPER_575-578\n"
    "# =============================================================================\n"
)
insert_text += const_section + "\n\n"
for b in blocks:
    insert_text += b + "\n\n"

# ── 5. Insert before __all__ ──────────────────────────────────────────────────
MARKER = "\n__all__ = ["
if MARKER not in cp4:
    sys.exit('ERROR: __all__ marker not found in CP4')

cp4_new = cp4.replace(MARKER, insert_text + MARKER, 1)

# ── 6. Update __all__ list ────────────────────────────────────────────────────
new_entries = (
    "    # --- Session 154: Universal Epoch / Periodic Table UQFF ---\n"
    '    "UniversalEpoch3DIPONuclearConvergenceCalculator",      # PAPER_573 (#161)\n'
    '    "DPMPyramidSumNuclearBindingPeriodicTableCalculator",   # PAPER_575 (#162)\n'
    '    "UQFFAtomicMassStandardModelErrorFactorCalculator",     # PAPER_576 (#163)\n'
    '    "IslandOfStability5thEpochSuperheavyElementsCalculator",# PAPER_577 (#164)\n'
    '    "UQFFCompEigenvalueQuantumGravityLinkageCalculator",    # PAPER_578 (#165)\n'
)
OLD_TAIL = (
    '    "AldersOlbersBSFGMetricGapAnalysisCalculator",            # PAPER_566 (#160)\n'
    "\n"
    "]"
)
NEW_TAIL = (
    '    "AldersOlbersBSFGMetricGapAnalysisCalculator",            # PAPER_566 (#160)\n'
    "\n"
    + new_entries
    + "]"
)
if OLD_TAIL not in cp4_new:
    sys.exit('ERROR: __all__ tail anchor not found')

cp4_new = cp4_new.replace(OLD_TAIL, NEW_TAIL, 1)

# ── 7. Write ──────────────────────────────────────────────────────────────────
with open('CondensedPhysics4.py', 'w', encoding='utf-8') as f:
    f.write(cp4_new)

print(f"Done. Lines: {cp4_new.count(chr(10))}")
