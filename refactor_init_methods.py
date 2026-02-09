#!/usr/bin/env python3
"""
Phase 2 Refactoring: Remove hardcoded system data from __init__ methods.
"""

import re

def main():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    changes = 0
    
    # Define the classes and their specific hardcoded patterns to remove
    # Format: (class_name, list of self.xxx patterns that are system-specific)
    classes_to_fix = [
        ('SpiralGalaxyUQFFCalculator', [
            'self.distance', 'self.r =', 'self.SFR', 'self.B_field', 'self.z =',
            'self.M_BH', 'self.sigma', 'self.age_Myr'
        ]),
        ('GiantSpiralGalaxyCalculator', [
            'self.distance', 'self.r =', 'self.SFR', 'self.B_field', 'self.z =',
            'self.M_BH', 'self.sigma', 'self.age_Myr'
        ]),
        ('InteractingGalaxyCalculator', [
            'self.distance', 'self.r =', 'self.SFR', 'self.B_field', 'self.z =',
            'self.M_BH', 'self.sigma', 'self.age_Myr'
        ]),
        ('EmissionNebulaCalculator', [
            'self.distance', 'self.r =', 'self.SFR', 'self.B_field', 'self.z =',
            'self.age_Myr', 'self.T_gas', 'self.v_radial'
        ]),
        ('WolfRayetNebulaCalculator', [
            'self.distance', 'self.r =', 'self.SFR', 'self.B_field', 'self.z =',
            'self.age_Myr', 'self.T_gas', 'self.v_radial'
        ]),
        ('StarFormingRegionCalculator', [
            'self.distance', 'self.r_pillar', 'self.SFR', 'self.B_field', 'self.z =',
            'self.age_Myr', 'self.T_gas', 'self.v_radial'
        ]),
        ('StarClusterDynamicsCalculator', [
            'self.distance', 'self.r_cluster', 'self.SFR', 'self.B_field', 'self.z =',
            'self.age_Myr', 'self.T_HII', 'self.v_radial', 'self.N_stars'
        ]),
        ('ClusterVirialCalculator', [
            'self.M_virial', 'self.R_virial', 'self.sigma'
        ]),
        ('AGNJetCalculator', [
            'self.M_BH', 'self.v_jet', 'self.L_jet', 'self.theta_jet'
        ]),
    ]
    
    for class_name, _ in classes_to_fix:
        # Find class and its __init__
        pattern = rf'(class {class_name}:.*?)(def __init__\(self\):)(.*?)(def \w+\(self)'
        match = re.search(pattern, content, re.DOTALL)
        
        if match:
            docstring = match.group(1)
            init_header = match.group(2)
            init_body = match.group(3)
            next_method = match.group(4)
            
            # Update docstring to be generic
            docstring = re.sub(r'NGC\s*\d+|Virgo|M87|Pillars|Westerlund|Carina', 
                              'generic system', docstring, flags=re.IGNORECASE)
            docstring = re.sub(r'KEY PARAMETERS \(from.*?\):',
                              'EXAMPLE PARAMETERS (pass via compute methods):',
                              docstring)
            docstring = re.sub(r'TRIADIC SOLUTIONS.*?©',
                              '©', docstring, flags=re.DOTALL)
            
            # Create minimal __init__ with only physical constants
            new_init_body = '''
        """Initialize with physical constants only. System data passed to compute methods."""
        # Physical constants
        self.G = 6.674e-11          # Gravitational constant (m³/kg/s²)
        self.c = 2.998e8            # Speed of light (m/s)
        self.M_sun = 1.989e30       # Solar mass (kg)
        self.ly = 9.461e15          # Light year (m)
        self.Mpc = 3.086e22         # Megaparsec (m)
        self.pc = 3.086e16          # Parsec (m)
        
        # UQFF calibration constants (universal, not system-specific)
        self.f_UA = 0.999           # Aether fraction
        self.f_SCm = 0.001          # Superconductive medium fraction
        self.R_EB = 1.0             # Energy balance ratio
        self.rho_vac_SCm = 7.09e-37 # Vacuum density [SCm] (J/m³)
        self.rho_vac_UA = 7.09e-36  # Vacuum density [UA] (J/m³)
        self.kappa = 0.0005         # day⁻¹ decay rate (universal)
        self.k_eta = 2.75e8         # LENR calibration (universal)
        self.Delta_k_eta = 7.25e8   # Buoyancy coefficient (universal)
        self.V_ratio = 1/33         # Boyle's Law V_little/V_big (universal)
        
    '''
            
            # Rebuild the class start
            new_class_start = docstring + init_header + new_init_body + next_method
            
            # Replace in content
            old_class_start = match.group(0)
            content = content.replace(old_class_start, new_class_start, 1)
            changes += 1
            print(f'Refactored: {class_name}')
    
    # Write the updated file
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f'\nTotal classes refactored: {changes}')
    return changes


if __name__ == '__main__':
    main()
