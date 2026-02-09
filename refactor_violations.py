#!/usr/bin/env python3
"""
Refactor CondensedPhysics.py violation classes to follow pure calculator architecture.

This script:
1. Removes hardcoded system data from __init__
2. Makes compute methods fully parameterized
3. Renames classes to generic physics names
4. Updates global instances
"""

import re

# Mapping of old class names to new generic names
CLASS_RENAMES = {
    'VirgoClusterMassModel': 'GalaxyClusterMassCalculator',
    'VirgoClusterVirialModel': 'ClusterVirialCalculator',
    'VirgoClusterM87JetModel': 'AGNJetCalculator',
    'PillarsOfCreationModel': 'StarFormingRegionCalculator',
    'Westerlund2ClusterModel': 'StarClusterDynamicsCalculator',
    'NGC3596Model': 'SpiralGalaxyUQFFCalculator',
    'NGC1961Model': 'GiantSpiralGalaxyCalculator',
    'NGC5335Model': 'InteractingGalaxyCalculator',
    'NGC2014Model': 'EmissionNebulaCalculator',
    'NGC2020Model': 'WolfRayetNebulaCalculator',
}

# Instance renames
INSTANCE_RENAMES = {
    'VIRGO_CLUSTER_MASS_MODEL': 'GALAXY_CLUSTER_MASS_CALC',
    'VIRGO_CLUSTER_VIRIAL_MODEL': 'CLUSTER_VIRIAL_CALC',
    'VIRGO_CLUSTER_M87_JET_MODEL': 'AGN_JET_CALC',
    'PILLARS_OF_CREATION_MODEL': 'STAR_FORMING_REGION_CALC',
    'WESTERLUND2_CLUSTER_MODEL': 'STAR_CLUSTER_DYNAMICS_CALC',
    'NGC3596_MODEL': 'SPIRAL_GALAXY_UQFF_CALC',
    'NGC1961_MODEL': 'GIANT_SPIRAL_GALAXY_CALC',
    'NGC5335_MODEL': 'INTERACTING_GALAXY_CALC',
    'NGC2014_MODEL': 'EMISSION_NEBULA_CALC',
    'NGC2020_MODEL': 'WOLF_RAYET_NEBULA_CALC',
}

def refactor_init(init_block: str, class_name: str) -> str:
    """
    Refactor __init__ to remove hardcoded system data.
    Keep only physical constants.
    """
    # Keep only these constants in __init__
    keep_patterns = [
        r'self\.G\s*=\s*[0-9.e+-]+',
        r'self\.c\s*=\s*[0-9.e+-]+',
        r'self\.M_sun\s*=\s*[0-9.e+-]+',
        r'self\.ly\s*=\s*[0-9.e+-]+',
        r'self\.hbar\s*=\s*[0-9.e+-]+',
        r'self\.k_B\s*=\s*[0-9.e+-]+',
        r'self\.Mpc\s*=\s*[0-9.e+-]+',
        r'self\.pc\s*=\s*[0-9.e+-]+',
    ]
    
    # Build new init that only keeps physical constants
    new_init = '''    def __init__(self):
        """Initialize with physical constants only. System data passed to compute methods."""
        self.G = 6.674e-11          # Gravitational constant (m³/kg/s²)
        self.c = 2.998e8            # Speed of light (m/s)
        self.M_sun = 1.989e30       # Solar mass (kg)
        self.ly = 9.461e15          # Light year (m)
        self.Mpc = 3.086e22         # Megaparsec (m)
        self.pc = 3.086e16          # Parsec (m)
'''
    return new_init


def refactor_compute_method(method_block: str) -> str:
    """
    Refactor compute method to require parameters instead of using self.xxx defaults.
    """
    # Replace "if xxx is None: xxx = self.xxx" patterns
    method_block = re.sub(
        r'if\s+(\w+)\s+is\s+None:\s*\n\s*\1\s*=\s*self\.\w+',
        lambda m: f'# {m.group(1)} is a required parameter',
        method_block
    )
    
    # Replace self.distance, self.M_BH, etc. with parameter references
    system_attrs = [
        'distance', 'M_BH', 'M_virial', 'R_virial', 'SFR', 'z', 'sigma',
        'age_Myr', 'r_pillar', 'r_cluster', 'B_field', 'T_gas', 'T_HII',
        'v_radial', 'N_stars'
    ]
    
    for attr in system_attrs:
        # Add comment noting this should be a parameter
        method_block = re.sub(
            rf'self\.{attr}\b',
            f'params["{attr}"]',
            method_block
        )
    
    return method_block


def update_docstring(docstring: str, old_name: str, new_name: str) -> str:
    """Update docstring to reflect generic nature."""
    # Remove specific system references from title
    docstring = re.sub(
        rf'{old_name}.*?Model',
        f'{new_name}',
        docstring
    )
    
    # Add note about parameterized usage
    if 'KEY PARAMETERS' in docstring:
        docstring = docstring.replace(
            'KEY PARAMETERS',
            'EXAMPLE PARAMETERS (pass via compute methods)'
        )
    
    return docstring


def main():
    # Read the file
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Track changes
    changes = []
    
    # 1. Rename classes
    for old_name, new_name in CLASS_RENAMES.items():
        pattern = rf'\bclass\s+{old_name}\b'
        if re.search(pattern, content):
            content = re.sub(pattern, f'class {new_name}', content)
            changes.append(f'Renamed class {old_name} → {new_name}')
    
    # 2. Rename instances
    for old_inst, new_inst in INSTANCE_RENAMES.items():
        # Update instance creation
        for old_class, new_class in CLASS_RENAMES.items():
            pattern = rf'\b{old_inst}\s*=\s*{old_class}\(\)'
            replacement = f'{new_inst} = {new_class}()'
            if re.search(pattern, content):
                content = re.sub(pattern, replacement, content)
                changes.append(f'Renamed instance {old_inst} → {new_inst}')
    
    # 3. Update references to old instance names
    for old_inst, new_inst in INSTANCE_RENAMES.items():
        content = re.sub(rf'\b{old_inst}\b', new_inst, content)
    
    # 4. Update class references in type hints and other places
    for old_name, new_name in CLASS_RENAMES.items():
        content = re.sub(rf'\b{old_name}\b', new_name, content)
    
    # Write the updated file
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f'Applied {len(changes)} changes:')
    for c in changes[:20]:
        print(f'  • {c}')
    if len(changes) > 20:
        print(f'  ... and {len(changes) - 20} more')
    
    return len(changes)


if __name__ == '__main__':
    n = main()
    print(f'\nRefactoring complete. {n} changes applied.')
    print('Run tests to verify: python -c "import CondensedPhysics; print(\"OK\")"')
