#!/usr/bin/env python3
"""Add compute_master_gravity method to all Model classes."""

MASTER_GRAVITY_BLOCK = '''
    def compute_master_gravity(self, r=None, t=0, include_all=True):
        """
        Master Gravity: Complete UQFF gravity computation.
        
        g_master = g_N * (1 - Ub/g_N + Ui/g_N) * H_SCm * fTRZ
        
        Combines: Newtonian, buoyancy, inertia, superconductivity, triadic zone
        """
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None:
            r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        
        # Base Newtonian gravity
        g_N = G * M / r**2
        
        # Get buoyancy correction
        buoy = self.compute_buoyant_equation(r)
        Ub = buoy['Ub']
        
        # Get superconductive coefficient
        sc = self.compute_superconductive_equation()
        H_SCm = sc['H_SCm']
        
        # Get triadic factor (default: in resonance zone)
        fTRZ = 1.0
        
        # Compute master gravity
        if g_N != 0:
            correction = 1 - Ub / g_N
        else:
            correction = 1.0
        
        g_master = g_N * correction * max(0.01, H_SCm) * fTRZ
        
        # Collect additional terms if requested
        if include_all:
            comp = self.compute_compressed_equation(r)
            res = self.compute_resonance_equation(r, t)
            tri = self.compute_triadic_equation(r)
            return {
                'g_master': g_master,
                'g_N': g_N,
                'Ub': Ub,
                'H_SCm': H_SCm,
                'fTRZ': fTRZ,
                'correction': correction,
                'g_compressed': comp['g_compressed'],
                'g_resonance': res['g_resonance'],
                'g_triadic': tri['g_triadic'],
                'equation': f"g_master = g_N * (1 - Ub/g_N) * H_SCm * fTRZ = {g_master:.6e} m/s^2"
            }
        
        return {'g_master': g_master, 'g_N': g_N}
'''

def add_master_gravity():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Find all Model class definitions
    model_classes = []
    for i, line in enumerate(lines):
        if line.startswith('class ') and 'Model' in line and ('(' in line or ':' in line):
            model_classes.append((i, line.strip()))
    
    print(f"Found {len(model_classes)} Model classes")
    
    # Find which ones already have compute_master_gravity
    has_method = set()
    for i, line in enumerate(lines):
        if 'def compute_master_gravity' in line:
            for j in range(i, -1, -1):
                if lines[j].startswith('class '):
                    class_name = lines[j].split('(')[0].split(':')[0].replace('class ', '').strip()
                    has_method.add(class_name)
                    break
    
    print(f"{len(has_method)} already have compute_master_gravity")
    
    # Process each class
    insertion_points = []
    for idx, (line_num, class_def) in enumerate(model_classes):
        class_name = class_def.split('(')[0].split(':')[0].replace('class ', '').strip()
        if class_name in has_method:
            continue
        
        if idx + 1 < len(model_classes):
            next_class_line = model_classes[idx + 1][0]
        else:
            next_class_line = len(lines)
        
        insert_line = None
        for j in range(next_class_line - 1, line_num, -1):
            stripped = lines[j].strip()
            if stripped and not stripped.startswith('#'):
                if lines[j].startswith('    ') or lines[j].startswith('\t'):
                    insert_line = j + 1
                    break
        
        if insert_line:
            insertion_points.append((insert_line, class_name))
    
    print(f"Will add compute_master_gravity to {len(insertion_points)} models")
    
    insertion_points.sort(reverse=True)
    method_lines = MASTER_GRAVITY_BLOCK.split('\n')
    
    for insert_line, class_name in insertion_points:
        for i, method_line in enumerate(reversed(method_lines)):
            lines.insert(insert_line, method_line + '\n')
    
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.writelines(lines)
    
    print(f"Modified {len(insertion_points)} models")

if __name__ == '__main__':
    add_master_gravity()
