#!/usr/bin/env python3
"""
Add 6 missing UQFF Master Equation methods to all 10 May 2025 models.

Missing equations:
1. compute_UQFF_base() - F_U = Σ(Ug_i + Ub_i) + Um
2. compute_superconductive_equation() - [SCm](t) = [SCm]₀ × e^(-κt) × (1 - cos(ωt))
3. compute_buoyant_equation() - F_U_Bi = Σ β_i × (Ug_i - ρ_vac,i)
4. compute_master_buoyant_equation() - F_U_Bi_i = Ω_g × (M_bh/d_g) × Σ(Ug_i + Ub_i)
5. compute_triadic_equation() - g(r,t) = Σ(i=1..26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]
6. compute_quadratic_equation() - Root solutions for g² + bg + c = 0
"""

# Template code for the 6 missing methods
METHODS_TEMPLATE = '''
    def compute_UQFF_base(self, t=None, r=None):
        """
        Compute UQFF Base Unified Field equation.
        
        F_U = Σ(Ug_i + Ub_i) + Um + Ur + Ut + Ui + UA + SCm
        
        This is the fundamental UQFF equation combining all field contributions.
        
        Returns:
            dict with all computed field values
        """
        import math
        if t is None:
            t = self.t_age
        if r is None:
            r = self.r
        
        # Gravitational terms
        g_grav = self.compute_g_grav(r)
        
        # Universal gravity components (26 layers simplified to primary terms)
        k_1, k_2, k_3, k_4 = 1.5, 1.2, 0.8, 0.5
        Ug1 = k_1 * g_grav                          # Magnetic dipole
        Ug2 = k_2 * g_grav * (1 + self.z)           # Charge-reactivity
        Ug3 = k_3 * g_grav * math.sin(2 * math.pi * t / self.t_age) if self.t_age > 0 else 0  # String rotation
        Ug4 = k_4 * g_grav * self.rho_vac_UA / 1e-36 if self.rho_vac_UA else k_4 * g_grav    # Vacuum concentration
        
        # Buoyancy term
        beta = CONSTANTS.get('beta_i', 0.6)
        Ub = beta * (g_grav - self.rho_vac_UA * self.G) if self.rho_vac_UA else beta * g_grav
        
        # Magnetism term
        Um = self.q_e * self.v_wind * self.B / self.m_p if hasattr(self, 'v_wind') else 0
        
        # Additional contributions (simplified)
        Ur = 0.01 * g_grav              # Radiation pressure
        Ut = 0.001 * g_grav             # Thermal contribution
        Ui = 0.0001 * g_grav            # Interstellar medium
        UA = self.rho_vac_UA * self.G * self.M / (r * r) if self.rho_vac_UA else 0
        SCm = self.rho_vac_SCm * self.G * self.M / (r * r) if self.rho_vac_SCm else 0
        
        # Total unified field
        F_U = Ug1 + Ug2 + Ug3 + Ug4 + Ub + Um + Ur + Ut + Ui + UA + SCm
        
        return {
            'F_U': F_U,
            'Ug1': Ug1,
            'Ug2': Ug2,
            'Ug3': Ug3,
            'Ug4': Ug4,
            'Ub': Ub,
            'Um': Um,
            'Ur': Ur,
            'Ut': Ut,
            'Ui': Ui,
            'UA': UA,
            'SCm': SCm,
            'g_grav': g_grav
        }
    
    def compute_superconductive_equation(self, t=None):
        """
        Compute UQFF Superconductive equation.
        
        [SCm](t) = [SCm]₀ × e^(-κt) × (1 - cos(ωt))
        
        Models SCm vacuum modulation over time.
        
        Returns:
            dict with superconductive field values
        """
        import math
        if t is None:
            t = self.t_age
        
        # UQFF calibrated constants
        kappa = CONSTANTS.get('kappa', 0.0005)      # Reactivity decay (day⁻¹)
        SSq = CONSTANTS.get('SSq', 0.57)            # Quantum state factor
        
        # Convert kappa from day⁻¹ to s⁻¹
        kappa_s = kappa / 86400
        
        # SCm base density
        SCm_0 = self.rho_vac_SCm if self.rho_vac_SCm else 7.09e-37
        
        # Angular frequency (based on characteristic timescale)
        omega_SCm = 2 * math.pi / (self.t_age if self.t_age > 0 else 1e13)
        
        # Exponential decay
        decay = math.exp(-kappa_s * t)
        
        # Oscillation factor
        oscillation = 1 - math.cos(omega_SCm * t)
        
        # SCm(t) modulation
        SCm_t = SCm_0 * decay * oscillation
        
        # Effective field contribution
        g_SCm = SCm_t * self.G * self.M / (self.r * self.r) if self.r > 0 else 0
        
        return {
            'SCm_t': SCm_t,
            'SCm_0': SCm_0,
            'decay_factor': decay,
            'oscillation_factor': oscillation,
            'g_SCm': g_SCm,
            'kappa': kappa,
            'omega_SCm': omega_SCm,
            'SSq': SSq
        }
    
    def compute_buoyant_equation(self, t=None, r=None):
        """
        Compute UQFF Buoyant equation (F_U_Bi).
        
        F_U_Bi = Σ β_i × (Ug_i - ρ_vac,i × G)
        
        Inside→Out atomic scale buoyancy (pushing outward).
        
        Returns:
            dict with buoyant field values
        """
        import math
        if t is None:
            t = self.t_age
        if r is None:
            r = self.r
        
        # Beta coupling factors
        beta = CONSTANTS.get('beta_i', 0.6)
        
        # Gravitational base
        g_grav = self.compute_g_grav(r)
        
        # Universal gravity components
        k_1 = 1.5
        Ug = k_1 * g_grav
        
        # Vacuum contribution
        rho_vac = self.rho_vac_UA if self.rho_vac_UA else 7.09e-36
        vac_term = rho_vac * self.G
        
        # Momentum term (electron pressure)
        m_e = CONSTANTS.get('m_e', 9.109e-31)
        c = self.c
        DPM_momentum = CONSTANTS.get('DPM_momentum', 1e-8)
        theta = math.pi / 4  # Default angle
        momentum_term = (m_e * c ** 2 / (r ** 2)) * DPM_momentum * math.cos(theta)
        
        # Gravity term
        DPM_gravity = CONSTANTS.get('DPM_gravity', 1e-10)
        gravity_term = g_grav * DPM_gravity
        
        # Buoyant fraction
        F_U_Bi = beta * (Ug - vac_term) + momentum_term + gravity_term
        
        return {
            'F_U_Bi': F_U_Bi,
            'beta': beta,
            'Ug': Ug,
            'vac_term': vac_term,
            'momentum_term': momentum_term,
            'gravity_term': gravity_term,
            'direction': 'Inside → Out',
            'scale': 'Atomic'
        }
    
    def compute_master_buoyant_equation(self, t=None):
        """
        Compute UQFF Master Buoyant equation (F_U_Bi_i).
        
        F_U_Bi_i = Ω_g × (M_bh/d_g) × Σ(Ug_i + Ub_i)
        
        Outside→In cosmic scale buoyancy (supporting from vacuum).
        
        Returns:
            dict with master buoyant field values
        """
        import math
        if t is None:
            t = self.t_age
        
        # Galactic parameters (defaults for star-forming regions)
        M_bh = CONSTANTS.get('M_bh_MW', 4e6 * self.M_sun)      # MW central BH
        d_g = CONSTANTS.get('d_g_sun', 2.55e20)                 # Distance to galactic center
        Omega_g = CONSTANTS.get('Omega_g_MW', 2.5e-8)           # Galactic rotation rate
        
        # Get base UQFF values
        base = self.compute_UQFF_base(t)
        Ug_sum = base['Ug1'] + base['Ug2'] + base['Ug3'] + base['Ug4']
        Ub = base['Ub']
        
        # Master buoyant integral
        galactic_factor = Omega_g * (M_bh / d_g)
        F_U_Bi_i = galactic_factor * (Ug_sum + Ub)
        
        return {
            'F_U_Bi_i': F_U_Bi_i,
            'galactic_factor': galactic_factor,
            'Omega_g': Omega_g,
            'M_bh': M_bh,
            'd_g': d_g,
            'Ug_sum': Ug_sum,
            'Ub': Ub,
            'direction': 'Outside → In',
            'scale': 'Cosmic'
        }
    
    def compute_triadic_equation(self, t=None, n_layers=26):
        """
        Compute UQFF Triadic equation (26-layer gravitational scaling).
        
        g(r,t) = Σ(i=1..26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i] × Q_i × [UA]_i × [SCm]_i
        
        Full 26-dimensional polynomial gravity summation.
        
        Returns:
            dict with triadic/26-layer field values
        """
        import math
        if t is None:
            t = self.t_age
        
        # Base gravitational acceleration
        g_grav = self.compute_g_grav()
        
        # Universal gravity couplings
        k_1, k_2, k_3, k_4 = 1.5, 1.2, 0.8, 0.5
        
        # 26-layer summation
        g_triadic = 0.0
        layer_contributions = []
        
        for i in range(1, n_layers + 1):
            # Layer-dependent quantum state factor
            Q_i = 1.0 / (1 + 0.01 * i)  # Decreasing with layer
            
            # Layer-dependent UA factor
            UA_i = self.rho_vac_UA * (1 - 0.02 * i) if self.rho_vac_UA else 7.09e-36 * (1 - 0.02 * i)
            
            # Layer-dependent SCm factor
            SCm_i = self.rho_vac_SCm * (1 - 0.01 * i) if self.rho_vac_SCm else 7.09e-37 * (1 - 0.01 * i)
            
            # Phase factor for this layer
            phase = 2 * math.pi * i / n_layers
            
            # Four Ug components for this layer
            Ug1_i = k_1 * g_grav * Q_i * UA_i / 7.09e-36
            Ug2_i = k_2 * g_grav * Q_i * (1 + self.z) * math.cos(phase)
            Ug3_i = k_3 * g_grav * Q_i * math.sin(phase + 2 * math.pi * t / self.t_age) if self.t_age > 0 else 0
            Ug4_i = k_4 * g_grav * Q_i * SCm_i / 7.09e-37
            
            # Total for this layer
            g_i = (Ug1_i + Ug2_i + Ug3_i + Ug4_i) * Q_i * UA_i * SCm_i * 1e72  # Normalize
            g_triadic += g_i
            
            layer_contributions.append({
                'layer': i,
                'Q_i': Q_i,
                'UA_i': UA_i,
                'SCm_i': SCm_i,
                'g_i': g_i
            })
        
        return {
            'g_triadic': g_triadic,
            'n_layers': n_layers,
            'g_grav_base': g_grav,
            'layer_sum': sum(lc['g_i'] for lc in layer_contributions),
            'first_layer': layer_contributions[0] if layer_contributions else None,
            'last_layer': layer_contributions[-1] if layer_contributions else None
        }
    
    def compute_quadratic_equation(self, t=None):
        """
        Compute UQFF Quadratic equation (root solutions).
        
        Solves g² + b×g + c = 0 for gravitational field regimes.
        
        Returns:
            dict with quadratic solution values
        """
        import math
        if t is None:
            t = self.t_age
        
        # Get compressed and resonance values
        compressed = self.compute_compressed_equation(t)
        resonance = self.compute_resonance_equation(t)
        
        g_comp = compressed.get('g_compressed', 0)
        R_amp = resonance.get('R_amplitude', 0)
        
        # Quadratic coefficients based on UQFF field balance
        # g² - (g_comp + R_amp)g + g_comp × R_amp × correlation = 0
        a = 1.0
        b = -(g_comp + R_amp)
        c = g_comp * R_amp * 0.1  # Correlation factor
        
        # Discriminant
        discriminant = b * b - 4 * a * c
        
        # Solutions
        if discriminant >= 0:
            sqrt_disc = math.sqrt(discriminant)
            x1 = (-b + sqrt_disc) / (2 * a)
            x2 = (-b - sqrt_disc) / (2 * a)
            solution_type = 'real'
        else:
            sqrt_disc = math.sqrt(-discriminant)
            x1 = -b / (2 * a)  # Real part
            x2 = sqrt_disc / (2 * a)  # Imaginary part
            solution_type = 'complex'
        
        return {
            'x1': x1,
            'x2': x2,
            'solution_type': solution_type,
            'discriminant': discriminant,
            'a': a,
            'b': b,
            'c': c,
            'g_compressed': g_comp,
            'R_amplitude': R_amp
        }
'''

def add_methods_to_models():
    """Add the 6 missing methods to all 10 models."""
    import re
    
    # Read the file
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find all 10 model classes
    models = [
        'NGC2264Model',
        'UGC10214Model',
        'NGC4676Model',
        'RedSpiderNebulaModel',
        'NGC3372Model',
        'AGCarinaeModel',
        'M42Model',
        'TarantulaNebulaModel',
        'NGC2841Model',
        'MysticMountainModel'
    ]
    
    # For each model, find the validate_model method and insert new methods before it
    for model in models:
        # Pattern to find "def validate_model(self)" within the class
        pattern = rf'(class {model}:.*?)(    def validate_model\(self.*?\):)'
        
        match = re.search(pattern, content, re.DOTALL)
        if match:
            # Insert the methods template before validate_model
            class_content = match.group(1)
            validate_method = match.group(2)
            
            # Check if methods already exist
            if 'compute_UQFF_base' in class_content:
                print(f'{model}: Methods already exist, skipping')
                continue
            
            # Insert methods
            new_content = class_content + METHODS_TEMPLATE + '\n' + validate_method
            content = content[:match.start()] + new_content + content[match.end():]
            print(f'{model}: Added 6 UQFF methods')
        else:
            print(f'{model}: Could not find validate_model, skipping')
    
    # Write back
    with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
        f.write(content)
    
    print('\nDone! Run validation to verify.')

if __name__ == '__main__':
    add_methods_to_models()
