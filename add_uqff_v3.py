#!/usr/bin/env python3
"""Add 8 UQFF Master Equation methods to all Model classes using line-based insertion."""

UQFF_BLOCK = '''
    # === 8 UQFF MASTER EQUATIONS (Auto-added) ===
    
    def compute_UQFF_base(self, r=None, t=0):
        """UQFF Base: F_U = Ug1 + Ug2 + Ug3 + Ug4 - Ub + Ui"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        Ug = G * M / r**2
        return {'Ug': Ug, 'F_U': Ug, 'equation': f"F_U = G*M/r^2 = {Ug:.3e} m/s^2"}
    
    def compute_compressed_equation(self, r=None):
        """Compressed MUGE gravity equation"""
        import numpy as np
        G, H_0 = 6.674e-11, 2.268e-18
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_N = G * M / r**2
        g_H = H_0**2 * r
        g = g_N - g_H
        return {'g_compressed': g, 'g_N': g_N, 'g_Hubble': g_H}
    
    def compute_resonance_equation(self, r=None, t=0):
        """Resonance MUGE with frequency modes"""
        import numpy as np
        G, H_0, c = 6.674e-11, 2.268e-18, 2.998e8
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        aDPM = G * M / r**2
        aExpFreq = H_0 * c * np.sin(H_0 * t)
        return {'g_resonance': aDPM + aExpFreq, 'aDPM': aDPM, 'aExpFreq': aExpFreq}
    
    def compute_superconductive_equation(self, T=None, T_c=None):
        """Superconductive coefficient H_SCm"""
        T = T if T else getattr(self, 'T', 4.2)
        T_c = T_c if T_c else getattr(self, 'T_c', 9.2)
        H_SCm = max(0, 1 - (T / T_c)**4) if T < T_c else 0
        return {'H_SCm': H_SCm, 'T': T, 'T_c': T_c}
    
    def compute_buoyant_equation(self, r=None, V=None):
        """Buoyancy component Ub"""
        import numpy as np
        G, rho_vac = 6.674e-11, 7.09e-36
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        if V is None: V = (4/3) * np.pi * r**3
        g = G * M / r**2
        Ub = (rho_vac * V * g) / M
        return {'Ub': Ub, 'g': g, 'V': V}
    
    def compute_master_buoyant_equation(self, r=None, v=0):
        """Master Buoyant: F_U_Bi_i = M * (Ug_i - Ub_i + Ui_i)"""
        import numpy as np
        G, rho_vac, k_eta, c = 6.674e-11, 7.09e-36, 1e-113, 2.998e8
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        V = (4/3) * np.pi * r**3
        Ug = G * M / r**2
        Ub = (rho_vac * V * Ug) / M
        Ui = -k_eta * rho_vac * V * (v / c)**2
        a_net = Ug - Ub + Ui
        return {'F_U_Bi_i': M * a_net, 'Ug': Ug, 'Ub': Ub, 'Ui': Ui, 'a_net': a_net}
    
    def compute_triadic_equation(self, r=None):
        """Triadic gravity: g_matter + g_field + g_vacuum"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_matter = G * M / r**2
        g_field, g_vacuum = 1e-15, 1e-20
        return {'g_triadic': g_matter + g_field + g_vacuum, 'g_matter': g_matter}
    
    def compute_quadratic_equation(self, r=None):
        """Quadratic polynomial fit: g = a*r^2 + b*r + c"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_base = G * M / r**2
        a, b, c = -1e-30, 1e-20, g_base
        return {'g_quad': a * r**2 + b * r + c, 'a': a, 'b': b, 'c': c}
'''

def add_uqff():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Find all Model class definitions and their line numbers
    model_classes = []
    for i, line in enumerate(lines):
        if line.startswith('class ') and 'Model' in line and ('(' in line or ':' in line):
            model_classes.append((i, line.strip()))
    
    print(f"Found {len(model_classes)} Model classes")
    
    # Find which ones already have compute_UQFF_base
    has_uqff = set()
    for i, line in enumerate(lines):
        if 'def compute_UQFF_base' in line:
            # Find the class this method belongs to
            for j in range(i, -1, -1):
                if lines[j].startswith('class '):
                    class_name = lines[j].split('(')[0].split(':')[0].replace('class ', '').strip()
                    has_uqff.add(class_name)
                    break
    
    print(f"{len(has_uqff)} already have UQFF")
    
    # Process each class that needs UQFF
    insertion_points = []
    for idx, (line_num, class_def) in enumerate(model_classes):
        class_name = class_def.split('(')[0].split(':')[0].replace('class ', '').strip()
        if class_name in has_uqff:
            continue
        
        # Find the end of this class (next class definition or end of file)
        if idx + 1 < len(model_classes):
            next_class_line = model_classes[idx + 1][0]
        else:
            next_class_line = len(lines)
        
        # Find the last non-empty, non-comment line before next class
        # that has proper class indentation (4 spaces = method level)
        insert_line = None
        for j in range(next_class_line - 1, line_num, -1):
            stripped = lines[j].strip()
            if stripped and not stripped.startswith('#'):
                # Check if it's inside the class (has some indentation)
                if lines[j].startswith('    ') or lines[j].startswith('\t'):
                    insert_line = j + 1  # Insert after this line
                    break
        
        if insert_line:
            insertion_points.append((insert_line, class_name))
    
    print(f"Will add UQFF to {len(insertion_points)} models")
    
    # Insert UQFF blocks from bottom to top (to preserve line numbers)
    insertion_points.sort(reverse=True)
    
    uqff_lines = UQFF_BLOCK.split('\n')
    
    for insert_line, class_name in insertion_points:
        # Insert the UQFF block
        for i, uqff_line in enumerate(reversed(uqff_lines)):
            lines.insert(insert_line, uqff_line + '\n')
        print(f"  Added UQFF to {class_name} at line {insert_line}")
    
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.writelines(lines)
    
    print(f"\nModified {len(insertion_points)} models")
    print("Done!")

if __name__ == '__main__':
    add_uqff()
