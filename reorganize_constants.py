#!/usr/bin/env python3
"""
Reorganize CONSTANTS in QCalc.py - distribute extracted constants into proper sections.
"""
import re

# Read QCalc.py
with open('QCalc.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Find and remove the "EXTRACTED FROM source*.js" block
extracted_block_pattern = re.compile(
    r"\n    # ═+\n    # EXTRACTED FROM source\*\.js.*?(?=\n\n\n# ═|\n}\n)",
    re.DOTALL
)

# Find the extracted block
match = extracted_block_pattern.search(content)
if match:
    extracted_block = match.group(0)
    print(f"Found extracted block: {len(extracted_block)} chars")
    
    # Parse constants from extracted block
    const_pattern = re.compile(r"'([^']+)':\s*([^,\n]+),")
    extracted_constants = dict(const_pattern.findall(extracted_block))
    print(f"Parsed {len(extracted_constants)} constants")
    
    # Remove the extracted block
    content = content[:match.start()] + content[match.end():]
    print("Removed extracted block")
else:
    print("ERROR: Could not find extracted block")
    exit(1)

# Define category mappings for the 256 constants
# Each will be added to appropriate existing section or new section

# Categorize constants
categories = {
    'MAGNETIC_FIELD_CONSTANTS': [],
    'MASS_REFERENCE_VALUES': [],
    'DISTANCE_SCALE_REFERENCES': [],
    'TIMESCALE_REFERENCES': [],
    'VELOCITY_REFERENCES': [],
    'ENERGY_LUMINOSITY_POWER': [],
    'FREQUENCY_OSCILLATION': [],
    'COUPLING_EXTENDED': [],
    'QUANTUM_UNCERTAINTY': [],
    'VACUUM_DENSITY_EXTENDED': [],
    'COSMOLOGICAL_EXTENDED': [],
    'MISCELLANEOUS_PHYSICS': []
}

for key, val in extracted_constants.items():
    # Magnetic
    if key.startswith(('B', 'magnetic', 'num_magnetic')):
        categories['MAGNETIC_FIELD_CONSTANTS'].append((key, val))
    # Mass
    elif key.startswith(('M_', 'Mbh', 'mass', 'ns_mass', 'ejecta', 'gas_mass', 'trap_mass', 'primary_mass', 'secondary_mass', 'progenitor', 'proton_mass')):
        categories['MASS_REFERENCE_VALUES'].append((key, val))
    # Distance
    elif key.startswith(('d_', 'dg', 'r_', 'radius', 'distance', 'separation', 'tidal_radius', 'length', 'size', 'L_jet')):
        categories['DISTANCE_SCALE_REFERENCES'].append((key, val))
    # Time
    elif key.startswith(('t_', 'tau_', 'T_', 'dt', 'period', 'age', 'evolution_timescale')):
        categories['TIMESCALE_REFERENCES'].append((key, val))
    # Velocity
    elif key.startswith(('v_', 'V_', 'relative_velocity', 'gas_v')):
        categories['VELOCITY_REFERENCES'].append((key, val))
    # Energy/Luminosity/Power
    elif key.startswith(('E_', 'L_', 'P_', 'explosion_energy')):
        categories['ENERGY_LUMINOSITY_POWER'].append((key, val))
    # Frequency
    elif key.startswith(('f', 'omega', 'fosc', 'fov', 'force')):
        categories['FREQUENCY_OSCILLATION'].append((key, val))
    # Coupling
    elif key.startswith(('k_', 'k1', 'k2', 'k3', 'k4', 'kappa')):
        categories['COUPLING_EXTENDED'].append((key, val))
    # Quantum
    elif key.startswith(('Delta', 'delta', 'epsilon', 'integral', 'sigma')):
        categories['QUANTUM_UNCERTAINTY'].append((key, val))
    # Vacuum/Density
    elif key.startswith(('rho', 'Evac', 'C_vac')):
        categories['VACUUM_DENSITY_EXTENDED'].append((key, val))
    # Cosmological
    elif key.startswith(('H', 'Lambda', 'Omega', 'z_')):
        categories['COSMOLOGICAL_EXTENDED'].append((key, val))
    else:
        categories['MISCELLANEOUS_PHYSICS'].append((key, val))

# Generate new section blocks
new_sections = ""

for section_name, constants in categories.items():
    if constants:
        new_sections += f"""
    # ═══════════════════════════════════════════════════════════════════════════
    # {section_name} (Extracted from source*.js - 163 files)
    # ═══════════════════════════════════════════════════════════════════════════
"""
        for key, val in sorted(constants):
            new_sections += f"    '{key}': {val},\n"

# Find insertion point - just before closing brace of CONSTANTS
# Look for the line with just "}"
closing_pattern = re.search(r"\n}\n\n\n# ═+\n# UQFF SCALE SYSTEM", content)
if closing_pattern:
    insert_pos = closing_pattern.start()
    content = content[:insert_pos] + new_sections + content[insert_pos:]
    print(f"Inserted {len(new_sections)} chars of organized constants")
else:
    print("ERROR: Could not find insertion point")
    exit(1)

# Write back
with open('QCalc.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("\n✓ Reorganized constants into category sections:")
for cat, items in sorted(categories.items(), key=lambda x: -len(x[1])):
    if items:
        print(f"  {cat}: {len(items)} constants")
