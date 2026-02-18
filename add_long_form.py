#!/usr/bin/env python3
"""Add long_form_solution method to all Model classes that don't have it."""

LONG_FORM_BLOCK = '''
    def long_form_solution(self, r=None, t=0):
        """
        Long-form solution with complete equation display.
        Returns detailed step-by-step calculation for all UQFF equations.
        """
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None:
            r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        
        # Collect all equation results
        results = {}
        steps = []
        
        # 1. UQFF Base
        uqff = self.compute_UQFF_base(r, t)
        results['UQFF_base'] = uqff
        steps.append(f"UQFF Base: F_U = G*M/r^2 = {uqff['Ug']:.6e} m/s^2")
        
        # 2. Compressed equation
        comp = self.compute_compressed_equation(r)
        results['compressed'] = comp
        steps.append(f"Compressed: g = g_N - g_H = {comp['g_compressed']:.6e} m/s^2")
        
        # 3. Resonance equation
        res = self.compute_resonance_equation(r, t)
        results['resonance'] = res
        steps.append(f"Resonance: g_res = aDPM + aExpFreq = {res['g_resonance']:.6e} m/s^2")
        
        # 4. Superconductive
        sc = self.compute_superconductive_equation()
        results['superconductive'] = sc
        steps.append(f"Superconductive: H_SCm = {sc['H_SCm']:.6f}")
        
        # 5. Buoyant equation
        buoy = self.compute_buoyant_equation(r)
        results['buoyant'] = buoy
        steps.append(f"Buoyancy: Ub = {buoy['Ub']:.6e} m/s^2")
        
        # 6. Master buoyant
        master = self.compute_master_buoyant_equation(r)
        results['master_buoyant'] = master
        steps.append(f"Master Buoyant: F_U_Bi_i = {master['F_U_Bi_i']:.6e} N")
        
        # 7. Triadic
        tri = self.compute_triadic_equation(r)
        results['triadic'] = tri
        steps.append(f"Triadic: g_triadic = {tri['g_triadic']:.6e} m/s^2")
        
        # 8. Quadratic
        quad = self.compute_quadratic_equation(r)
        results['quadratic'] = quad
        steps.append(f"Quadratic: g_quad = {quad['g_quad']:.6e} m/s^2")
        
        return {
            'model': self.__class__.__name__,
            'parameters': {'M': M, 'r': r, 't': t},
            'results': results,
            'long_form': '\\n'.join(steps)
        }
'''

def add_long_form():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # Find all Model class definitions
    model_classes = []
    for i, line in enumerate(lines):
        if line.startswith('class ') and 'Model' in line and ('(' in line or ':' in line):
            model_classes.append((i, line.strip()))
    
    print(f"Found {len(model_classes)} Model classes")
    
    # Find which ones already have long_form_solution
    has_long_form = set()
    for i, line in enumerate(lines):
        if 'def long_form_solution' in line:
            for j in range(i, -1, -1):
                if lines[j].startswith('class '):
                    class_name = lines[j].split('(')[0].split(':')[0].replace('class ', '').strip()
                    has_long_form.add(class_name)
                    break
    
    print(f"{len(has_long_form)} already have long_form_solution")
    
    # Process each class that needs long_form_solution
    insertion_points = []
    for idx, (line_num, class_def) in enumerate(model_classes):
        class_name = class_def.split('(')[0].split(':')[0].replace('class ', '').strip()
        if class_name in has_long_form:
            continue
        
        # Find the end of this class
        if idx + 1 < len(model_classes):
            next_class_line = model_classes[idx + 1][0]
        else:
            next_class_line = len(lines)
        
        # Find last method line before next class
        insert_line = None
        for j in range(next_class_line - 1, line_num, -1):
            stripped = lines[j].strip()
            if stripped and not stripped.startswith('#'):
                if lines[j].startswith('    ') or lines[j].startswith('\t'):
                    insert_line = j + 1
                    break
        
        if insert_line:
            insertion_points.append((insert_line, class_name))
    
    print(f"Will add long_form_solution to {len(insertion_points)} models")
    
    # Insert from bottom to top
    insertion_points.sort(reverse=True)
    long_form_lines = LONG_FORM_BLOCK.split('\n')
    
    for insert_line, class_name in insertion_points:
        for i, lf_line in enumerate(reversed(long_form_lines)):
            lines.insert(insert_line, lf_line + '\n')
        print(f"  Added long_form_solution to {class_name}")
    
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.writelines(lines)
    
    print(f"\nModified {len(insertion_points)} models")

if __name__ == '__main__':
    add_long_form()
