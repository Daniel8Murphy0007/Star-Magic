#!/usr/bin/env python3
"""Generate uqff_derivation_output.txt from Trinity files"""

import sys
sys.path.insert(0, '.')

from _uqff_derivation_equations import DerivationRegistry
from _uqff_reference_documentation import ReferenceDocumentationRegistry

output = []
output.append('='*80)
output.append('UQFF DERIVATION EQUATIONS OUTPUT REGISTRY')
output.append('Generated from: _uqff_derivation_equations.py + _uqff_reference_documentation.py')
output.append('='*80)
output.append('')

for i, eq in enumerate(DerivationRegistry.DERIVATIONS, 1):
    ref = ReferenceDocumentationRegistry.PAPERS[i-1] if i <= len(ReferenceDocumentationRegistry.PAPERS) else None
    output.append(f'[{i:03d}] {eq.name}')
    output.append(f'     Equation: {eq.formula}')
    if ref:
        output.append(f'     Sources: {" | ".join(ref.papers)}')
    output.append(f'     CVW: {eq.cvw}')
    output.append('')

output.append('')
output.append('='*80)
output.append(f'Total Derivations: {len(DerivationRegistry.DERIVATIONS)}')
output.append('='*80)

with open('uqff_derivation_output.txt', 'w', encoding='utf-8') as f:
    f.write('\n'.join(output))

print(f'✓ Generated uqff_derivation_output.txt ({len(DerivationRegistry.DERIVATIONS)} equations)')
