#!/usr/bin/env python3
"""
Comprehensive MAIN_1.cpp Repair Script
Fixes: Unicode characters, uncommented documentation, missing main()
"""

import re
import sys

print("MAIN_1.cpp Comprehensive Repair")
print("=" * 50)

# Read file
try:
    with open('MAIN_1.cpp', 'r', encoding='utf-8', errors='replace') as f:
        content = f.read()
except FileNotFoundError:
    print("ERROR: MAIN_1.cpp not found!")
    sys.exit(1)

print("✓ File loaded")

# Step 1: Replace Unicode characters (including problematic apostrophes)
# First pass: handle multi-byte characters
content = content.replace('\u2018', "'").replace('\u2019', "'")  # Smart single quotes
content = content.replace('\u201c', '"').replace('\u201d', '"')  # Smart double quotes  
content = content.replace('\u2013', '-').replace('\u2014', '-')  # En/em dashes

unicode_map = {
    '×': '*', '–': '-', '—': '-', '−': '-',
    'Δ': 'Delta', 'δ': 'delta', 'ρ': 'rho', 'ω': 'omega',
    'Ω': 'Omega', 'Λ': 'Lambda', 'λ': 'lambda', 'π': 'pi',
    'α': 'alpha', 'β': 'beta', 'γ': 'gamma', 'μ': 'mu',
    'ν': 'nu', 'τ': 'tau', 'σ': 'sigma', 'φ': 'phi',
    'ψ': 'psi', 'θ': 'theta', 'ε': 'epsilon', 'η': 'eta',
    '∫': 'integral', '∝': 'proportional_to', '√': 'sqrt',
    '∞': 'infinity', '≈': 'approx', '±': '+/-', '…': '...',
    '³': '^3', '²': '^2', '¹': '^1', '⁰': '^0',
}

unicode_count = 0
for old, new in unicode_map.items():
    count = content.count(old)
    if count > 0:
        content = content.replace(old, new)
        unicode_count += count

print(f"✓ Replaced {unicode_count} Unicode characters")

# Step 2: Fix uncommented documentation
lines = content.split('\n')
fixed_lines = []
in_comment_block = False
blocks_commented = 0

for i, line in enumerate(lines):
    # Detect uncommented documentation patterns
    is_doc_line = (
        re.match(r'^\d+\.\s+\*\*', line) or  # "1. **Title**:"
        re.match(r'^#{1,4}\s+', line) or      # "### Header"
        (not line.strip().startswith('//') and 
         not line.strip().startswith('/*') and 
         not line.strip().startswith('*') and
         not line.strip().startswith('#include') and
         not line.strip().startswith('#define') and
         not line.strip().startswith('using') and
         len(line.strip()) > 100 and
         ('This comment' in line or 'Elaboration' in line or 
          'Overview of' in line or 'Derivation from' in line or
          'Variables and' in line or 'Long-Form Calculations' in line or
          'Significance and' in line))
    )
    
    if is_doc_line:
        if not in_comment_block:
            fixed_lines.append('/*')
            in_comment_block = True
            blocks_commented += 1
        # Clean markdown
        clean_line = re.sub(r'\*\*(.*?)\*\*', r'\1', line)
        clean_line = re.sub(r'^#{1,4}\s+', '', clean_line)
        fixed_lines.append(clean_line)
    else:
        # Close block if we hit a comment or code line
        if in_comment_block and (line.strip().startswith('//') or 
                                  line.strip().startswith('#') or
                                  line.strip() == ''):
            fixed_lines.append('*/')
            in_comment_block = False
        fixed_lines.append(line)

if in_comment_block:
    fixed_lines.append('*/')

content = '\n'.join(fixed_lines)
print(f"✓ Commented {blocks_commented} documentation blocks")

# Step 3: Add main() if missing
if 'int main(' not in content:
    main_func = '''

// ============================================================================
// MAIN FUNCTION - UQFF 26-Layer Gravity Calculator
// ============================================================================

int main() {
    std::cout << "========================================\\n";
    std::cout << "  MAIN_1 UQFF Calculator\\n";
    std::cout << "  26-Layer Gravity Framework\\n";
    std::cout << "========================================\\n\\n";
    
    std::cout << "Compressed UQFF equation:\\n";
    std::cout << "g(r,t) = sum_{i=1 to 26} (Ug1_i + Ug2_i + Ug3_i + Ug4i_i)\\n\\n";
    
    // Physical constants
    const double G = 6.6743e-11;      // Gravitational constant (m^3 kg^-1 s^-2)
    const double h_bar = 1.0546e-34;  // Reduced Planck constant (J s)
    const double c = 2.99792458e8;    // Speed of light (m/s)
    const double rho_vac_UA = 7.09e-36; // Universal Aether vacuum density (J/m^3)
    
    // Sample system parameters (Solar-like)
    double r = 1.0e10;        // 10 billion meters
    double M = 1.989e30;      // Solar mass (kg)
    double omega0 = 1.0e-6;   // Angular frequency (s^-1)
    double t = 0.0;           // Time (s)
    
    std::cout << "Sample Parameters:\\n";
    std::cout << "  r = " << r << " m\\n";
    std::cout << "  M = " << M << " kg\\n";
    std::cout << "  omega_0 = " << omega0 << " s^-1\\n";
    std::cout << "  t = " << t << " s\\n\\n";
    
    // Calculate 26-layer gravity
    double g_total = 0.0;
    std::cout << "Calculating 26 layers...\\n";
    
    for (int i = 1; i <= 26; i++) {
        // Layer-specific parameters
        double r_i = r / static_cast<double>(i);
        double Q_i = static_cast<double>(i);
        double SCm_i = static_cast<double>(i * i);
        double f_TRZ_i = 1.0 / static_cast<double>(i);
        double f_Um_i = static_cast<double>(i);
        double alpha_i = 0.01;  // DPM stability factor
        
        // E_DPM,i = (h_bar * c / r_i^2) * Q_i * SCm_i
        double E_DPM_i = (h_bar * c / (r_i * r_i)) * Q_i * SCm_i;
        
        // Ug1_i = E_DPM,i / r_i^2 * [UA]_i * f_TRZ_i
        double Ug1_i = (E_DPM_i / (r_i * r_i)) * rho_vac_UA * f_TRZ_i;
        
        // Ug2_i = E_DPM,i / r_i^2 * [SCm]_i * f_Um_i
        double Ug2_i = (E_DPM_i / (r_i * r_i)) * SCm_i * f_Um_i;
        
        // Ug3_i = (h_bar * omega_i / 2) * Q_i * cos(2 * PI * f_i * t) / r_i
        double omega_i = omega0;  // Simplified: same for all layers
        double f_i = omega0 / (2.0 * M_PI);
        double Ug3_i = (h_bar * omega_i / 2.0) * Q_i * std::cos(2.0 * M_PI * f_i * t) / r_i;
        
        // Ug4i_i = (G * M_i / r_i^2) * (1 + alpha_i) * [SCm]_i
        double M_i = M / static_cast<double>(i);
        double Ug4i_i = (G * M_i / (r_i * r_i)) * (1.0 + alpha_i) * SCm_i;
        
        // Sum all Ug terms for this layer
        double g_layer = Ug1_i + Ug2_i + Ug3_i + Ug4i_i;
        g_total += g_layer;
        
        if (i <= 3 || i == 26) {
            std::cout << "  Layer " << i << ": " << g_layer << " m/s^2\\n";
        } else if (i == 4) {
            std::cout << "  ...\\n";
        }
    }
    
    std::cout << "\\nResults:\\n";
    std::cout << "  g_total = " << g_total << " m/s^2\\n";
    std::cout << "  For comparison: Earth surface g = 9.81 m/s^2\\n";
    std::cout << "\\n========================================\\n";
    std::cout << "Calculation complete!\\n";
    std::cout << "========================================\\n";
    
    return 0;
}
'''
    content += main_func
    print("✓ Added main() function")
else:
    print("✓ main() function already present")

# Write repaired file
with open('MAIN_1_repaired.cpp', 'w', encoding='utf-8') as f:
    f.write(content)

print("\n" + "=" * 50)
print("SUCCESS: Created MAIN_1_repaired.cpp")
print("=" * 50)
print("\nNext steps:")
print("  1. Test compile:")
print("     g++ -std=c++17 -c MAIN_1_repaired.cpp")
print("  2. Build executable:")
print("     g++ -std=c++17 MAIN_1_repaired.cpp -o MAIN_1.exe")
print("  3. Run:")
print("     ./MAIN_1.exe")
print("  4. If successful, replace original:")
print("     mv MAIN_1_repaired.cpp MAIN_1.cpp")
