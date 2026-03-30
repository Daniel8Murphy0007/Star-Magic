"""Inject olbers_insert.txt before __all__ in CondensedPhysics4.py."""
import pathlib

cp4_path   = pathlib.Path('CondensedPhysics4.py')
insert_path = pathlib.Path('olbers_insert.txt')

content = cp4_path.read_text(encoding='utf-8')
insert  = insert_path.read_text(encoding='utf-8')

MARKER = '\n__all__ = ['
pos = content.rfind(MARKER)   # rfind = last occurrence (only one, but safe)
if pos == -1:
    raise RuntimeError('__all__ not found in CondensedPhysics4.py')

# Insert the new classes+constants right before __all__
new_content = content[:pos] + '\n\n' + insert + MARKER + content[pos + len(MARKER):]

# Also add the 3 class names near the end of __all__
OLD_LAST = '    "BSFGBohrSommerfeldAetherQuantizationCalculator",         # PAPER_562 (#157)\n\n]'
NEW_LAST = (
    '    "BSFGBohrSommerfeldAetherQuantizationCalculator",         # PAPER_562 (#157)\n'
    '    # --- Session 153: Alders/Olbers Paradox Resolution ---\n'
    '    "AldersOlbersParadoxDPMShellFluxCalculator",              # PAPER_564 (#158)\n'
    '    "AldersOlbersVDSNumberSystemResolutionCalculator",         # PAPER_565 (#159)\n'
    '    "AldersOlbersBSFGMetricGapAnalysisCalculator",            # PAPER_566 (#160)\n'
    '\n]'
)
if OLD_LAST not in new_content:
    raise RuntimeError('Could not find last __all__ entry to append')

new_content = new_content.replace(OLD_LAST, NEW_LAST, 1)

cp4_path.write_text(new_content, encoding='utf-8')
print('Done. Lines:', new_content.count('\n') + 1)
