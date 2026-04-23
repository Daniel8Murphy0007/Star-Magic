from dpm_helpers import dpm_ug1_seed, dpm_ug2_shell
import re

# 8 UQFF Master Equation methods template (adapted for any model)
UQFF_METHODS = '''
    # 
    # 8 UQFF MASTER EQUATIONS
    # 
    
    def compute_UQFF_base(self, params=None) -> dict:
        """UQFF Base: F_U = Σ(Ug_i + Ub_i) + Um + Ur + Ut + Ui + UA + SCm"""
        G = CONSTANTS.get('G', 6.6743e-11)
        M = getattr(self, 'M', getattr(params, 'M', 1.989e30)) if params else getattr(self, 'M', 1.989e30)
        r = getattr(self, 'r', getattr(params, 'r', 1e10)) if params else getattr(self, 'r', 1e10)
        g_base = G * M / (r ** 2) if r > 0 else 0
        Ug1 = 1.5 * g_base; Ug2 = 1.2 * g_base; Ug3 = 0.8 * g_base * np.sin(np.pi/4); Ug4 = 0.5 * g_base
        beta = CONSTANTS.get('beta_i', 0.6)
        rho_vac = CONSTANTS.get('rho_vac_UA', 7.09e-36)
        Ub = beta * (g_base - rho_vac * G)
        Um = CONSTANTS.get('mu_0', 4e-7) * 1e6 / (4 * np.pi * r**2) if r > 0 else 0
        F_U = Ug1 + Ug2 + Ug3 + Ug4 + Ub + Um
        return {'F_U': F_U, 'Ug1': Ug1, 'Ug2': Ug2, 'Ug3': Ug3, 'Ug4': Ug4, 'Ub': Ub, 'Um': Um, 'g_base': g_base}
    
    def compute_compressed_equation(self, params=None, t=None) -> dict:
        """UQFF Compressed: g = g_N  (1+Ht)  (1+M_sf)  (1-E_rad)  (1+f_TRZ) + a_EM"""
        G = CONSTANTS.get('G', 6.6743e-11)
        M = getattr(self, 'M', getattr(params, 'M', 1.989e30)) if params else getattr(self, 'M', 1.989e30)
        r = getattr(self, 'r', getattr(params, 'r', 1e10)) if params else getattr(self, 'r', 1e10)
        t = t if t is not None else getattr(self, 't', 3.156e13)
        g_N = G * M / (r ** 2) if r > 0 else 0
        H_0 = CONSTANTS.get('H_0', 2.268e-18)
        f_TRZ = CONSTANTS.get('f_TRZ', 0.1)
        hubble = 1 + H_0 * t; trz = 1 + f_TRZ
        g_compressed = g_N * hubble * trz
        return {'g_compressed': g_compressed, 'g_N': g_N, 'hubble_factor': hubble, 'trz_factor': trz}
    
    def compute_resonance_equation(self, params=None, t=None) -> dict:
        """UQFF Resonance: R = Σ(a_i  cos(ω_i  t))"""
        t = t if t is not None else getattr(self, 't', 1e8)
        omega_grav = 2 * np.pi / 3.156e13
        omega_mag = 2 * np.pi / 3.156e11
        R_grav = CONSTANTS.get('G', 6.6743e-11) * 1e30 / 1e20 * np.cos(omega_grav * t)
        R_mag = 1e-5 * np.cos(omega_mag * t)
        UA_factor = 1 + CONSTANTS.get('rho_vac_UA', 7.09e-36) / CONSTANTS.get('rho_vac_SCm', 7.09e-37)
        R_total = R_grav + R_mag * UA_factor
        return {'R_total': R_total, 'R_grav': R_grav, 'R_mag': R_mag, 'UA_factor': UA_factor}
    
    def compute_superconductive_equation(self, params=None, t=None) -> dict:
        """UQFF Superconductive: [SCm](t) = [SCm]  e^(-κt)  (1 - cos(ωt))"""
        t = t if t is not None else getattr(self, 't', 1e8)
        SCm_0 = CONSTANTS.get('rho_vac_SCm', 7.09e-37)
        kappa = CONSTANTS.get('kappa', 0.0005) / 86400
        omega = 2 * np.pi / 3.156e13
        SCm_t = SCm_0 * np.exp(-kappa * t) * (1 - np.cos(omega * t))
        return {'SCm_t': SCm_t, 'SCm_0': SCm_0, 'decay': np.exp(-kappa * t)}
    
    def compute_F_U_Bi_simple(self, params=None) -> dict:
        """UQFF Buoyant (F_U_Bi) — Form B: stellar scale. F_Bi = β × (Ug - ρ_vac × G)"""
        G = CONSTANTS.get('G', 6.6743e-11)
        M = getattr(self, 'M', getattr(params, 'M', 1.989e30)) if params else getattr(self, 'M', 1.989e30)
        r = getattr(self, 'r', getattr(params, 'r', 1e10)) if params else getattr(self, 'r', 1e10)
        g_base = G * M / (r ** 2) if r > 0 else 0
        beta = CONSTANTS.get('beta_i', 0.6)
        rho_vac = CONSTANTS.get('rho_vac_UA', 7.09e-36)
        F_Bi = beta * (g_base - rho_vac * G)
        return {'F_Bi': F_Bi, 'g_base': g_base, 'beta': beta, 'rho_vac': rho_vac}
    
    def compute_F_U_Bi_resonant(self, params=None) -> dict:
        """UQFF Master Buoyant (F_U_Bi_i) — Form C-2: resonant/TRZ. F_U_Bi_i = F_Bi × (1+f_TRZ) / (1-Omega_g)"""
        buoyant = self.compute_F_U_Bi_simple(params)
        Omega_g = CONSTANTS.get('Omega_g', 7.3e-16)
        f_TRZ = CONSTANTS.get('f_TRZ', 0.1)
        F_U_Bi_i = buoyant['F_Bi'] * (1 + f_TRZ) / (1 - Omega_g)
        return {'F_U_Bi_i': F_U_Bi_i, 'F_Bi': buoyant['F_Bi'], 'Omega_g': Omega_g, 'f_TRZ': f_TRZ}
    
    def compute_triadic_equation(self, params=None, n_layers=26) -> dict:
        """UQFF Triadic: g_26 = Σ(i=1 to 26) [Ug1_i + Ug2_i + Ug3_i + Ug4_i]"""
        G = CONSTANTS.get('G', 6.6743e-11)
        M = getattr(self, 'M', getattr(params, 'M', 1.989e30)) if params else getattr(self, 'M', 1.989e30)
        r = getattr(self, 'r', getattr(params, 'r', 1e10)) if params else getattr(self, 'r', 1e10)
        g_base = G * M / (r ** 2) if r > 0 else 0
        k_values = [1.5, 1.2, 0.8, 0.5]
        g_layers = []; total = 0
        for i in range(1, n_layers + 1):
            scale = 1 + 0.1 * (i - 13) / 13
            Ug_i = sum(k * g_base * scale for k in k_values)
            g_layers.append(Ug_i); total += Ug_i
        return {'g_triadic': total, 'g_layers': g_layers, 'n_layers': n_layers}
    
    def compute_quadratic_equation(self, params=None) -> dict:
        """UQFF Quadratic: Solves g + bg + c = 0 for unified field roots"""
        G = CONSTANTS.get('G', 6.6743e-11)
        M = getattr(self, 'M', getattr(params, 'M', 1.989e30)) if params else getattr(self, 'M', 1.989e30)
        r = getattr(self, 'r', getattr(params, 'r', 1e10)) if params else getattr(self, 'r', 1e10)
        a = 1; b = -G * M / (r ** 2) if r > 0 else 1e-10; c = CONSTANTS.get('Lambda', 1.1e-52) * r**2 / 3
        disc = b**2 - 4*a*c
        if disc >= 0:
            g_plus = (-b + np.sqrt(disc)) / (2*a); g_minus = (-b - np.sqrt(disc)) / (2*a)
        else:
            g_plus = g_minus = -b / (2*a)
        return {'g_plus': g_plus, 'g_minus': g_minus, 'discriminant': disc, 'a': a, 'b': b, 'c': c}
'''

# Read file
with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Models to update with their last method before end of class
models_to_fix = [
    ('NegativeTimeModel', 'def is_retrocausal'),
    ('AetherVacuumEnergyModel', 'def compute_vacuum_repulsion'),
    ('CosmicEggModel', 'def compute_26D_egg_energy'),
    ('SgrAStarGravityModel', 'def compute_smbh_gravity'),
    ('RetrocausalModel', 'def compute_retrocausal_evolution'),
    ('TRZModel', 'def is_TRZ_active'),
    ('VoidOscillationModel', 'def compute_void_expansion_rate'),
    ('TimeVaryingVacuumModel', 'def compute_combined_vacuum_energy'),
]

added = 0
for model_name, last_method in models_to_fix:
    # Find class block
    class_pattern = rf'^class {model_name}:'
    class_match = re.search(class_pattern, content, re.MULTILINE)
    if not class_match:
        print(f"Could not find class {model_name}")
        continue
    
    # Check if already has 8 methods
    next_class = re.search(r'^class \w+:', content[class_match.end():], re.MULTILINE)
    class_end = class_match.end() + next_class.start() if next_class else len(content)
    class_content = content[class_match.start():class_end]
    
    if 'def compute_UQFF_base' in class_content:
        print(f"{model_name}: Already has UQFF methods")
        continue
    
    # Find the last method and its ending
    last_method_match = re.search(rf'{last_method}\(.*?\).*?(?=\n    def |\n\nclass |\Z)', class_content, re.DOTALL)
    if not last_method_match:
        print(f"{model_name}: Could not find {last_method}")
        continue
    
    # Insert UQFF methods after the last method
    insert_pos = class_match.start() + last_method_match.end()
    content = content[:insert_pos] + UQFF_METHODS + content[insert_pos:]
    added += 1
    print(f"{model_name}: Added 8 UQFF methods")

# Write back
with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f"\nTotal models updated: {added}")
