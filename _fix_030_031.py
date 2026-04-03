import re

def fix_file(path, header):
    content = open(path, 'r', encoding='utf-8').read()

    abstracts = [m.start() for m in re.finditer(r'^## Abstract$', content, re.MULTILINE)]
    if len(abstracts) <= 1:
        print(f'Already clean: {path}')
        return

    # Take content from first Abstract to start of second Abstract
    first_copy = content[abstracts[0] : abstracts[1]]

    # Strip garbled artifacts
    first_copy = re.sub(r'\.Groups\[1\]\.Value\s*\n', '', first_copy)
    first_copy = re.sub(r'"PAPER_\{0:D3\}"\s*-f\s*\$n\s*\n', '', first_copy)
    first_copy = re.sub(r'\s{4}"PAPER_\{0:D3\}"\s*-f\s*\$n\s*\n', '', first_copy)
    first_copy = first_copy.rstrip()

    # Remove trailing `---` separator that was between copies
    first_copy = re.sub(r'\n\n---\s*$', '', first_copy)

    # Check if UQFF Production Framework appendix is already in first copy
    has_prod_appendix = '## Appendix: UQFF Production Framework Reference' in first_copy

    # Extract it from the full document (last occurrence)
    prod_appendix = ''
    if not has_prod_appendix:
        m = re.search(r'(## Appendix: UQFF Production Framework Reference.*)', content, re.DOTALL)
        if m:
            prod_appendix = '\n\n---\n\n' + m.group(1).rstrip()

    result = header.rstrip() + '\n\n' + first_copy + prod_appendix + '\n'
    open(path, 'w', encoding='utf-8').write(result)
    print(f'Fixed {path}: {len(content)} -> {len(result)} chars, {len(abstracts)} abstracts -> 1')

header_030 = """# PAPER #30b \u2014 Dark Sector Mediators in UQFF

**Title:** Dark Sector Mediator Constraints from LFV B\u2070 \u2192 K*\u2070 \u03c4\u00b1e\u2213 Searches via the Unified Quantum Field Framework

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (\u03ba = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**arXiv Reference:** 2506.15347 (LFV B\u2070 \u2192 K*\u2070 \u03c4\u00b1e\u2213, LHCb 5.4 fb\u207b\u00b9)  
**Validator:** `bsm_physics_validation.py` \u2014 PASSED  
**Index Slot:** \u00a71.4 BSM Physics, Paper #30  

---"""

header_031 = """# PAPER #31b \u2014 Flavor Anomalies Resolution via UQFF

**Title:** Resolution of B-Physics Flavor Anomalies at Future e\u207ae\u207b Factories: UQFF Predictions for the ECFA Higgs Factory Program

**Author:** Daniel T. Murphy  
**Framework:** UQFF Star-Magic (\u03ba = 0.0005/day, [SSq] = 0.57)  
**Date:** March 7, 2026  
**arXiv Reference:** 2506.15390 (ECFA Higgs factory study, e\u207ae\u207b colliders)  
**Supporting Data:** 2506.15256 (Belle II |V_cb|, LFU ratio), 2506.15347 (LFV limits)  
**Validator:** `bsm_physics_validation.py` \u2014 PASSED  
**Index Slot:** \u00a71.4 BSM Physics, Paper #31  

---"""

fix_file('whitepapers/PAPER_030_Dark_Sector_Mediators_UQFF.md', header_030)
fix_file('whitepapers/PAPER_031_Flavor_Anomalies_Resolution_UQFF.md', header_031)
