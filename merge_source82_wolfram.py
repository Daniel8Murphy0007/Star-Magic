#!/usr/bin/env python3
"""
Merge SMBH and Virgo Cluster source82_wolfram.cpp versions
C++20/MSVC 19.44+ compliant
Combines 15 SMBH classes + 9 unique Virgo classes = 24 total
"""

import re
from datetime import datetime

def read_file(filename):
    """Read file with utf-8 encoding"""
    with open(filename, 'r', encoding='utf-8', errors='ignore') as f:
        return f.read()

def extract_classes(content, start_marker, end_marker):
    """Extract class definitions between markers"""
    classes = []
    class_pattern = r'class\s+(\w+)'
    
    # Find all class definitions
    for match in re.finditer(class_pattern, content):
        class_name = match.group(1)
        if class_name == 'PhysicsTerm':
            continue  # Skip base class
        
        # Find class start - go back to find any comment block
        class_line_start = match.start()
        
        # Look backwards for comment block (// ========)
        comment_start = class_line_start
        lines_before = content[:class_line_start].split('\n')
        for i in range(len(lines_before) - 1, max(0, len(lines_before) - 20), -1):
            if '========' in lines_before[i] or 'CLASS' in lines_before[i]:
                comment_start = len('\n'.join(lines_before[:i]))
                break
        
        class_start = comment_start
        
        # Find class end (search for matching braces)
        brace_count = 0
        in_class = False
        class_end = class_start
        
        for i in range(match.end(), len(content)):
            if content[i] == '{':
                brace_count += 1
                in_class = True
            elif content[i] == '}':
                brace_count -= 1
                if in_class and brace_count == 0:
                    # Look for closing semicolon
                    class_end = i + 1
                    while class_end < len(content) and content[class_end] in ' \t\n;':
                        if content[class_end] == ';':
                            class_end += 1
                            break
                        class_end += 1
                    break
        
        class_code = content[class_start:class_end]
        classes.append((class_name, class_code.strip()))
    
    return classes

def generate_merged_file(smbh_classes, virgo_classes):
    """Generate unified C++20-compliant source82_wolfram.cpp"""
    
    # Filter out duplicate SMBHMSigmaRelationTerm from Virgo (keep SMBH version)
    virgo_unique = [(name, code) for name, code in virgo_classes if name != 'SMBHMSigmaRelationTerm']
    
    header = f"""// ============================================================================
// UNIFIED SOURCE82_WOLFRAM.CPP - DUAL-SCALE PHYSICS
// ============================================================================
// Merged: December 4, 2025
// Part 1: SMBH Local-Scale Physics (15 classes, 10^6-10^9 M_sun)
// Part 2: Virgo Cluster Cosmological-Scale Physics (9 unique classes, 10^15 M_sun)
// Total Classes: 24 (15 SMBH + 9 Virgo, 1 duplicate removed)
// Compiler: C++20 / MSVC 19.44.35207+
// Copyright: Daniel T. Murphy
// ============================================================================

#include <cmath>
#include <string>
#include <map>
#include <memory>
#include <vector>
#include <complex>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Physical constants (C++20 constexpr)
constexpr double G_CONST = 6.6743e-11;      // Gravitational constant (m³/kg·s²)
constexpr double M_SUN = 1.989e30;          // Solar mass (kg)
constexpr double MPC_TO_M = 3.086e22;       // Megaparsec to meters
constexpr double KPC_TO_M = 3.086e19;       // Kiloparsec to meters
constexpr double KEV_TO_J = 1.602e-16;      // keV to Joules
constexpr double K_BOLTZ = 1.381e-23;       // Boltzmann constant (J/K)
constexpr double M_PROTON = 1.673e-27;      // Proton mass (kg)
constexpr double C_LIGHT = 2.998e8;         // Speed of light (m/s)
constexpr double YEAR_TO_S = 3.156e7;       // Year to seconds

// ============================================================================
// BASE PHYSICS TERM INTERFACE
// ============================================================================

class PhysicsTerm {{
public:
    virtual ~PhysicsTerm() = default;  // C++20 defaulted virtual destructor
    virtual double compute(double t, const std::map<std::string, double>& params) const = 0;
    virtual std::string getName() const = 0;
    virtual std::string getDescription() const = 0;
    virtual bool validate(const std::map<std::string, double>& params) const {{ return true; }}
}};

// ============================================================================
// PART 1: SMBH LOCAL-SCALE PHYSICS (15 CLASSES)
// Scale: 10^6 - 10^9 M_sun
// Physics: Black hole M-σ relation, vacuum energy, quantum coupling, UQFF layers
// ============================================================================

"""
    
    # Add SMBH classes
    smbh_section = "\n\n".join([code for name, code in smbh_classes if name != 'PhysicsTerm'])
    
    virgo_header = """

// ============================================================================
// PART 2: VIRGO CLUSTER COSMOLOGICAL-SCALE PHYSICS (9 UNIQUE CLASSES)
// Scale: 10^15 M_sun cluster mass
// Physics: Galaxy cluster dynamics, ICM, dark matter, X-ray emission
// Distance: 16.5 Mpc, Virial radius: 2.2 Mpc, Velocity dispersion: 700 km/s
// ============================================================================

"""
    
    # Add unique Virgo classes
    virgo_section = "\n\n".join([code for name, code in virgo_unique])
    
    # Generate registration function (C++20 std::make_unique)
    all_class_names = [name for name, code in smbh_classes if name != 'PhysicsTerm'] + \
                      [name for name, code in virgo_unique]
    
    registration_func = f"""

// ============================================================================
// UNIFIED REGISTRATION FUNCTION (C++20 std::make_unique)
// ============================================================================

/*
void registerWolframTerms_source82(PhysicsTermRegistry& registry) {{
    // Part 1: SMBH Local-Scale Terms (15)
"""
    
    for name in [n for n, c in smbh_classes if n != 'PhysicsTerm']:
        reg_name = name.replace('Term', '').replace('SMBH', 'SMBH_')
        registration_func += f'    registry.registerPhysicsTerm("{reg_name}", std::make_unique<{name}>(), "wolfram-smbh");\n'
    
    registration_func += '\n    // Part 2: Virgo Cluster Cosmological-Scale Terms (9)\n'
    
    for name in [n for n, c in virgo_unique]:
        reg_name = name.replace('Term', '').replace('VirgoCluster', 'Virgo_')
        registration_func += f'    registry.registerPhysicsTerm("{reg_name}", std::make_unique<{name}>(), "wolfram-virgo");\n'
    
    registration_func += f"""    
    std::cout << "Registered 24 Wolfram terms from source82_wolfram.cpp (15 SMBH + 9 Virgo)" << std::endl;
}}
*/

// ============================================================================
// SUMMARY: 24 UNIFIED PHYSICS TERM CLASSES
// ============================================================================
// SMBH LOCAL-SCALE (15 classes):
"""
    
    for i, (name, code) in enumerate([c for c in smbh_classes if c[0] != 'PhysicsTerm'], 1):
        registration_func += f'//   {i}. {name}\n'
    
    registration_func += '\n// VIRGO CLUSTER COSMOLOGICAL-SCALE (9 classes):\n'
    
    for i, (name, code) in enumerate(virgo_unique, 16):
        registration_func += f'//   {i}. {name}\n'
    
    registration_func += f"""//
// Duplicate removed: SMBHMSigmaRelationTerm (kept SMBH version from workspace)
// ============================================================================
// Integration: Include in MAIN_1_CoAnQi.cpp and uncomment registration function
// Build: MSVC 19.44+ with /std:c++20
// ============================================================================
"""
    
    # Combine all sections
    merged_content = header + smbh_section + virgo_header + virgo_section + registration_func
    
    return merged_content

def main():
    print("=" * 80)
    print("SOURCE82_WOLFRAM.CPP MERGE UTILITY")
    print("C++20/MSVC 19.44+ Compliant")
    print("=" * 80)
    print()
    
    # Read both versions
    print("[1/5] Reading workspace SMBH version...")
    smbh_content = read_file('source82_wolfram.cpp')
    print(f"      ✓ Loaded {len(smbh_content)} bytes")
    
    print("[2/5] Reading Virgo Cluster version...")
    virgo_content = read_file('source82_wolfram_VIRGO_EXTRACTION.cpp')
    print(f"      ✓ Loaded {len(virgo_content)} bytes")
    
    # Extract classes
    print("[3/5] Extracting class definitions...")
    smbh_classes = extract_classes(smbh_content, "SOURCE82-SPECIFIC", "")
    virgo_classes = extract_classes(virgo_content, "CLASS 820", "SUMMARY")
    
    print(f"      ✓ SMBH classes: {len(smbh_classes)}")
    print(f"      ✓ Virgo classes: {len(virgo_classes)}")
    
    # Check for duplicates
    smbh_names = [name for name, code in smbh_classes]
    virgo_names = [name for name, code in virgo_classes]
    duplicates = set(smbh_names) & set(virgo_names)
    
    if duplicates:
        print(f"      ⚠ Duplicate classes found: {', '.join(duplicates)}")
        print(f"      → Keeping SMBH version, removing from Virgo")
    
    # Generate merged file
    print("[4/5] Generating unified C++20-compliant file...")
    merged_content = generate_merged_file(smbh_classes, virgo_classes)
    
    # Write output
    print("[5/5] Writing merged source82_wolfram.cpp...")
    with open('source82_wolfram_MERGED.cpp', 'w', encoding='utf-8') as f:
        f.write(merged_content)
    
    print(f"      ✓ Written {len(merged_content)} bytes")
    print()
    print("=" * 80)
    print("MERGE COMPLETE")
    print("=" * 80)
    print(f"Total classes: {len(smbh_classes) + len(virgo_classes) - len(duplicates)}")
    print(f"  - SMBH local-scale: {len(smbh_classes)}")
    print(f"  - Virgo cosmological-scale: {len(virgo_classes) - len(duplicates)}")
    print(f"  - Duplicates removed: {len(duplicates)}")
    print()
    print("Output: source82_wolfram_MERGED.cpp")
    print("Next: Review merged file, then rename to source82_wolfram.cpp")
    print()

if __name__ == '__main__':
    main()
