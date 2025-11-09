#!/usr/bin/env python3
"""
Fix MAIN_1.cpp encoding and syntax issues
- Replace Unicode characters with ASCII equivalents
- Comment out markdown-style documentation
- Fix unescaped quotes in comments
"""

import re

# Read the file
with open('MAIN_1.cpp', 'r', encoding='utf-8', errors='replace') as f:
    content = f.read()

# Unicode character replacements
replacements = {
    '×': '*',  # multiplication sign
    '─': '-',  # em dash  
    '—': '-',  # em dash variant
    '–': '-',  # en dash
    'Δ': 'Delta',  # delta
    'δ': 'delta',  # small delta
    'ρ': 'rho',  # rho
    'ω': 'omega',  # omega
    'Ω': 'Omega',  # capital omega
    'Λ': 'Lambda',  # lambda
    'λ': 'lambda',  # small lambda
    'π': 'pi',  # pi
    'α': 'alpha',  # alpha
    'β': 'beta',  # beta
    'γ': 'gamma',  # gamma
    'μ': 'mu',  # mu
    'ν': 'nu',  # nu
    'τ': 'tau',  # tau
    'σ': 'sigma',  # sigma
    'φ': 'phi',  # phi
    'ψ': 'psi',  # psi
    'θ': 'theta',  # theta
    'ε': 'epsilon',  # epsilon
    'η': 'eta',  # eta
    '∫': 'integral',  # integral
    '∝': 'proportional_to',  # proportional to
    '∈': 'in',  # element of
    '√': 'sqrt',  # square root
    '∞': 'infinity',  # infinity
    '≈': 'approx',  # approximately equal
    '≠': '!=',  # not equal
    '≤': '<=',  # less than or equal
    '≥': '>=',  # greater than or equal
    '±': '+/-',  # plus-minus
    '→': '->',  # arrow
    '←': '<-',  # left arrow
    '⊗': 'x',  # tensor product
    '∇': 'nabla',  # nabla
    '∂': 'd',  # partial derivative
    '∑': 'sum',  # summation
    '∏': 'product',  # product
    '⁰': '^0', '¹': '^1', '²': '^2', '³': '^3',
    '⁴': '^4', '⁵': '^5', '⁶': '^6', '⁷': '^7',
    '⁸': '^8', '⁹': '^9',
    '₀': '_0', '₁': '_1', '₂': '_2', '₃': '_3',
    '₄': '_4', '₅': '_5', '₆': '_6', '₇': '_7',
    '₈': '_8', '₉': '_9',
    '"': '"',  # smart quotes
    '"': '"',  # smart quotes
    ''': "'",  # smart apostrophe
    ''': "'",  # smart apostrophe
    '…': '...',  # ellipsis
}

# Apply replacements
for old, new in replacements.items():
    content = content.replace(old, new)

# Fix lines that start with numbers and markdown (e.g., "1. **UQFF Core**:")
# These should be wrapped in comments
lines = content.split('\n')
fixed_lines = []
in_uncommented_block = False

for i, line in enumerate(lines):
    # Check if line starts with number followed by period and markdown
    if re.match(r'^\d+\.\s+\*\*', line):
        if not in_uncommented_block:
            fixed_lines.append('/*')
            in_uncommented_block = True
        # Remove the markdown formatting
        fixed_line = re.sub(r'\*\*(.*?)\*\*', r'\1', line)
        fixed_lines.append(fixed_line)
    elif re.match(r'^#{1,4}\s+', line):  # Markdown headers
        if not in_uncommented_block:
            fixed_lines.append('/*')
            in_uncommented_block = True
        # Remove the markdown header
        fixed_line = re.sub(r'^#{1,4}\s+', '', line)
        fixed_lines.append(fixed_line)
    else:
        # Check if we should close the comment block
        if in_uncommented_block and line.strip().startswith('//'):
            fixed_lines.append('*/')
            in_uncommented_block = False
            fixed_lines.append(line)
        else:
            fixed_lines.append(line)

# Close any open comment block
if in_uncommented_block:
    fixed_lines.append('*/')

content = '\n'.join(fixed_lines)

# Add main() function at the end if not present
if 'int main()' not in content and 'int main(' not in content:
    main_function = '''

// Main function for testing
int main() {
    std::cout << "MAIN_1 UQFF Calculator\\n";
    std::cout << "Compressed UQFF equation: g(r,t) = sum_{i=1 to 26} (Ug1_i + Ug2_i + Ug3_i + Ug4i_i)\\n";
    
    // Example calculation parameters
    double r = 1e10;  // 10 billion meters
    double M = 1e30;  // Solar mass in kg
    double omega0 = 1e-6;  // Angular frequency
    double t = 0.0;  // Time
    
    std::cout << "Sample parameters:\\n";
    std::cout << "  r = " << r << " m\\n";
    std::cout << "  M = " << M << " kg\\n";
    std::cout << "  omega_0 = " << omega0 << " s^-1\\n";
    std::cout << "  t = " << t << " s\\n";
    
    // TODO: Implement actual g(r,t) calculation
    std::cout << "\\nCalculation framework initialized.\\n";
    
    return 0;
}
'''
    content += main_function

# Write the fixed file
with open('MAIN_1_fixed.cpp', 'w', encoding='utf-8') as f:
    f.write(content)

print("✓ Created MAIN_1_fixed.cpp with fixes:")
print("  - Replaced Unicode characters with ASCII")
print("  - Commented out markdown-style documentation")
print("  - Added main() function if missing")
print("\nRun: g++ -std=c++17 MAIN_1_fixed.cpp -o MAIN_1")
