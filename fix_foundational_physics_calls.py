"""
Quick fix script to update all foundational physics method calls in UQFF calculators.

This script performs search-replace operations to fix:
1. Method name mismatches
2. Result extraction from EquationResult objects
"""

# All replacements needed (old_str -> new_str):
replacements = {
    # Floyd Sweet fixes
    "floyd_result = self.floyd_calc.compute_vacuum_density(t)\n            rho_vac_UA_t = floyd_result.result":
        "floyd_result = self.floyd_calc.compute_time_varying_density(rho_vac_UA_base, t)\n            rho_vac_UA_t = floyd_result.result",
    
    "floyd_result = self.floyd_calc.compute_vacuum_density(t)\n            rho_vac_t = floyd_result.result":
        "floyd_result = self.floyd_calc.compute_time_varying_density(rho_vac_base, t)\n            rho_vac_t = floyd_result.result",
    
    "floyd_result = self.floyd_calc.compute_vacuum_density(t)\n            rho_A_t = floyd_result.result":
        "floyd_result = self.floyd_calc.compute_time_varying_density(rho_A_base, t)\n            rho_A_t = floyd_result.result",
    
    # Heisenberg fixes - these are correct, just need result extraction
    
    # Cosmic Egg fixes
    "egg_result = self.egg_calc.compute_dimension_volume(i, t, R)\n                V_i_factor = egg_result":
        "V_0 = (4/3) * np.pi * (R ** 3)\n                egg_result = self.egg_calc.compute_layer_volume(i, V_0, t)\n                V_i_factor = egg_result.result",
    
    "egg_result = self.egg_calc.compute_dimension_volume(1, t, R)\n            volume_factor = egg_result":
        "V_0 = (4/3) * np.pi * (R ** 3)\n            egg_result = self.egg_calc.compute_layer_volume(1, V_0, t)\n            volume_factor = egg_result.result",
    
    "egg_result = self.egg_calc.compute_dimension_volume(1, t, R)\n            V_sys = (4/3) * np.pi * (R ** 3) * egg_result":
        "V_0 = (4/3) * np.pi * (R ** 3)\n            egg_result = self.egg_calc.compute_layer_volume(1, V_0, t)\n            V_sys = V_0 * egg_result.result",
    
    "egg_result = self.egg_calc.compute_dimension_volume(1, t, R)\n            V_t = (4/3) * np.pi * (R ** 3) * egg_result  # Volume with breathing\n            volume_factor = egg_result":
        "V_0 = (4/3) * np.pi * (R ** 3)\n            egg_result = self.egg_calc.compute_layer_volume(1, V_0, t)\n            V_t = V_0 * egg_result.result  # Volume with breathing\n            volume_factor = egg_result.result",
    
    # Negative Time fixes
    "trz_result = self.neg_time_calc.compute_TRZ_operator(t_n_i)\n                TRZ_i = trz_result.result":
        "trz_result = self.neg_time_calc.compute_negative_time_operator(t_n_i)\n                TRZ_i = trz_result.result",
    
    "trz_result = self.neg_time_calc.compute_TRZ_operator(t_n)\n            TRZ_amp = trz_result.result":
        "trz_result = self.neg_time_calc.compute_negative_time_operator(t_n)\n            TRZ_amp = trz_result.result",
    
    "trz_result = self.neg_time_calc.compute_TRZ_operator(t_n)\n            TRZ_factor = trz_result.result":
        "trz_result = self.neg_time_calc.compute_negative_time_operator(t_n)\n            TRZ_factor = trz_result.result",
}

print("Foundational Physics Method Fix Script")
print("=" * 80)
print(f"\nTotal replacements needed: {len(replacements)}")
print("\nReplacements:")
for i, (old, new) in enumerate(replacements.items(), 1):
    print(f"\n{i}. {old[:60]}...")
    print(f"   → {new[:60]}...")

print("\n" + "=" * 80)
print("Execute these replacements in QCalc.py using multi_replace_string_in_file")
