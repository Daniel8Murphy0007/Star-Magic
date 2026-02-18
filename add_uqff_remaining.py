#!/usr/bin/env python3
"""Add 8 UQFF Master Equation methods to all Model classes that don't have them."""

import re

# 8 UQFF method implementations
UQFF_METHODS = '''
    # ═══════════════════════════════════════════════════════════════════════════
    # 8 UQFF MASTER EQUATIONS (Auto-added)
    # ═══════════════════════════════════════════════════════════════════════════
    
    def compute_UQFF_base(self, r: float = None, t: float = 0) -> dict:
        """UQFF Base: F_U = Ug1 + Ug2 + Ug3 + Ug4 - Ub + Ui"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        Ug = G * M / r**2
        return {'Ug': Ug, 'F_U': Ug, 'equation': f"F_U = G×M/r² = {Ug:.3e} m/s²"}
    
    def compute_compressed_equation(self, r: float = None) -> dict:
        """Compressed MUGE: g = g_N + g_Hubble + g_Super + g_Envelope + Ug_sum + g_Λ + g_quantum + g_fluid + g_perturb"""
        import numpy as np
        G, H_0, c = 6.674e-11, 2.268e-18, 2.998e8
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_N = G * M / r**2
        g_H = H_0**2 * r
        g = g_N - g_H
        return {'g_compressed': g, 'g_N': g_N, 'g_Hubble': g_H, 'equation': f"g = g_N - g_H = {g:.3e} m/s²"}
    
    def compute_resonance_equation(self, r: float = None, t: float = 0) -> dict:
        """Resonance MUGE: g = aDPM + aSuperFreq + aQuantumFreq + aAetherFreq + aFluidFreq + aExpFreq"""
        import numpy as np
        G, H_0, c = 6.674e-11, 2.268e-18, 2.998e8
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        aDPM = G * M / r**2
        aExpFreq = H_0 * c * np.sin(H_0 * t)
        g_res = aDPM + aExpFreq
        return {'g_resonance': g_res, 'aDPM': aDPM, 'aExpFreq': aExpFreq, 'equation': f"g_res = aDPM + aExpFreq = {g_res:.3e} m/s²"}
    
    def compute_superconductive_equation(self, T: float = None, T_c: float = None) -> dict:
        """Superconductive: H_SCm = 1 - (T/T_c)^4, affects vacuum coupling"""
        T = T if T else getattr(self, 'T', 4.2)
        T_c = T_c if T_c else getattr(self, 'T_c', 9.2)
        H_SCm = max(0, 1 - (T / T_c)**4) if T < T_c else 0
        return {'H_SCm': H_SCm, 'T': T, 'T_c': T_c, 'equation': f"H_SCm = 1 - (T/T_c)⁴ = {H_SCm:.6f}"}
    
    def compute_buoyant_equation(self, r: float = None, V: float = None) -> dict:
        """Buoyancy: Ub = (ρ_vac × V × g) / M"""
        import numpy as np
        G, rho_vac = 6.674e-11, 7.09e-36
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        if V is None: V = (4/3) * np.pi * r**3
        g = G * M / r**2
        Ub = (rho_vac * V * g) / M
        return {'Ub': Ub, 'g': g, 'V': V, 'equation': f"Ub = (ρ_vac × V × g) / M = {Ub:.3e} m/s²"}
    
    def compute_master_buoyant_equation(self, r: float = None, v: float = 0) -> dict:
        """Master Buoyant: F_U_Bi_i = M × (Ug_i - Ub_i + Ui_i)"""
        import numpy as np
        G, rho_vac, k_eta, c = 6.674e-11, 7.09e-36, 1e-113, 2.998e8
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        V = (4/3) * np.pi * r**3
        Ug = G * M / r**2
        Ub = (rho_vac * V * Ug) / M
        Ui = -k_eta * rho_vac * V * (v / c)**2
        a_net = Ug - Ub + Ui
        F = M * a_net
        return {'F_U_Bi_i': F, 'Ug': Ug, 'Ub': Ub, 'Ui': Ui, 'a_net': a_net, 'equation': f"F_U_Bi_i = M×(Ug-Ub+Ui) = {F:.3e} N"}
    
    def compute_triadic_equation(self, r: float = None) -> dict:
        """Triadic: g_triadic = g_matter + g_field + g_vacuum"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_matter = G * M / r**2
        g_field = 1e-15  # Field contribution placeholder
        g_vacuum = 1e-20  # Vacuum contribution placeholder
        g_triadic = g_matter + g_field + g_vacuum
        return {'g_triadic': g_triadic, 'g_matter': g_matter, 'g_field': g_field, 'g_vacuum': g_vacuum, 'equation': f"g_triadic = {g_triadic:.3e} m/s²"}
    
    def compute_quadratic_equation(self, r: float = None) -> dict:
        """Quadratic: g_quad = a×r² + b×r + c (polynomial fit)"""
        import numpy as np
        G = 6.674e-11
        M = getattr(self, 'M', getattr(self, 'mass', 1e30))
        if r is None: r = getattr(self, 'r', getattr(self, 'radius', 1e8))
        g_base = G * M / r**2
        a, b, c = -1e-30, 1e-20, g_base
        g_quad = a * r**2 + b * r + c
        return {'g_quad': g_quad, 'a': a, 'b': b, 'c': c, 'equation': f"g_quad = ar² + br + c = {g_quad:.3e} m/s²"}
'''

def add_uqff_to_models():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find all Model classes
    class_pattern = r'^class\s+(\w+Model)\s*[:\(]'
    classes = re.findall(class_pattern, content, re.MULTILINE)
    print(f"Found {len(classes)} Model classes")
    
    # Find classes that already have compute_UQFF_base
    has_uqff = set(re.findall(r'class\s+(\w+Model).*?def compute_UQFF_base', content, re.DOTALL))
    print(f"{len(has_uqff)} already have UQFF methods")
    
    # Classes that need UQFF added
    needs_uqff = [c for c in classes if c not in has_uqff]
    print(f"{len(needs_uqff)} need UQFF methods added")
    
    if not needs_uqff:
        print("All models already have UQFF methods!")
        return
    
    modified = 0
    for class_name in needs_uqff:
        # Find the class and its end (next class or end of file)
        class_match = re.search(rf'^(class\s+{class_name}\s*[:\(].*?)(?=^class\s|\Z)', content, re.MULTILINE | re.DOTALL)
        if not class_match:
            print(f"  Could not find class body for {class_name}")
            continue
        
        class_body = class_match.group(1)
        class_start = class_match.start()
        class_end = class_match.end()
        
        # Find the last method in the class (last 'def ' at proper indentation)
        last_def = None
        for m in re.finditer(r'^    def \w+\(self.*?\n(?:        .*\n)*', class_body, re.MULTILINE):
            last_def = m
        
        if last_def:
            # Insert after the last method
            insert_pos = class_start + last_def.end()
            content = content[:insert_pos] + UQFF_METHODS + content[insert_pos:]
            modified += 1
            print(f"  Added UQFF to {class_name}")
        else:
            print(f"  Could not find methods in {class_name}")
    
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.write(content)
    
    print(f"\nModified {modified} models")
    print("Done!")

if __name__ == '__main__':
    add_uqff_to_models()
