"""Fix PAPER numbering in CondensedPhysics4.py: 563/564/565 -> 564/565/566"""
content = open('CondensedPhysics4.py', encoding='utf-8-sig').read()

# Fix class 1 docstring + PAPER value: 563 -> 564
content = content.replace('CP4 #158 \u2014 PAPER_563:', 'CP4 #158 \u2014 PAPER_564:')
content = content.replace('codebase resolution  PAPER_563', 'codebase resolution  PAPER_564')

# Fix class 2 docstring: 564 -> 565
content = content.replace('CP4 #159 \u2014 PAPER_564:', 'CP4 #159 \u2014 PAPER_565:')
content = content.replace('codebase resolution  PAPER_564', 'codebase resolution  PAPER_565')

# Fix class 3 docstring: was not changed yet (still says 565)
content = content.replace('CP4 #160 \u2014 PAPER_565:', 'CP4 #160 \u2014 PAPER_566:')
# codebase resolution PAPER_565 was already overwritten to 566 by powershell... check
# The docstring line for class 3: "codebase resolution  PAPER_5XX" 
# After PS it will be PAPER_566 (3rd pass changed 565->566 correctly)
# So leave PAPER_566 for class 3 as-is

# Fix the PAPER = '...' values
# After PS all three are PAPER_566. Fix first two:
content = content.replace("PAPER   = 'PAPER_566'", 'PAPER_564_SENTINEL', 1)
content = content.replace("PAPER   = 'PAPER_566'", 'PAPER_565_SENTINEL', 1)
# Third instance remains as PAPER_566 (correct for class 3)
content = content.replace('PAPER_564_SENTINEL', "PAPER   = 'PAPER_564'")
content = content.replace('PAPER_565_SENTINEL', "PAPER   = 'PAPER_565'")

# Cross-references in class 3 docstring
content = content.replace(
    'AldersOlbersParadoxDPMShellFluxCalculator (#158, PAPER_563)',
    'AldersOlbersParadoxDPMShellFluxCalculator (#158, PAPER_564)'
)
content = content.replace(
    'AldersOlbersVDSNumberSystemResolutionCalculator (#159, PAPER_564)',
    'AldersOlbersVDSNumberSystemResolutionCalculator (#159, PAPER_565)'
)

# gap_analysis dict cross-refs
content = content.replace('CP4 #158 PAPER_563', 'CP4 #158 PAPER_564')
content = content.replace('CP4 #159 PAPER_564', 'CP4 #159 PAPER_565')

# __all__ comments
content = content.replace('# PAPER_563 (#158)', '# PAPER_564 (#158)')
content = content.replace('# PAPER_564 (#159)', '# PAPER_565 (#159)')
content = content.replace('# PAPER_565 (#160)', '# PAPER_566 (#160)')
content = content.replace('PAPER_563\u2013565', 'PAPER_564-566')
content = content.replace('PAPER_563-565', 'PAPER_564-566')

# VMI2 reference in gap_analysis
content = content.replace(
    "'Olbers_DPM_shells':     'AldersOlbersParadoxDPMShellFluxCalculator CP4 #158 PAPER_563",
    "'Olbers_DPM_shells':     'AldersOlbersParadoxDPMShellFluxCalculator CP4 #158 PAPER_564"
)
content = content.replace(
    "'Olbers_VDS_number_sys': 'AldersOlbersVDSNumberSystemResolutionCalculator CP4 #159 PAPER_564",
    "'Olbers_VDS_number_sys': 'AldersOlbersVDSNumberSystemResolutionCalculator CP4 #159 PAPER_565"
)

# __all__ session comment
content = content.replace(
    '# --- Session 153: Alders/Olbers Paradox Resolution \u2014 PAPER_563\u2013565 ---',
    '# --- Session 153: Alders/Olbers Paradox Resolution \u2014 PAPER_564-566 ---'
)

open('CondensedPhysics4.py', 'w', encoding='utf-8').write(content)
print('Renumbering done')
