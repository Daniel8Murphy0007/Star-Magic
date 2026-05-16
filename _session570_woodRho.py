"""S570: Wood (pine) density 500 kg/m^3 = SO5^2 * D_phys + SO5^2 (EXACT)"""
D_p=4; SO=10
v=SO*SO*D_p+SO*SO
target=500
print(f"Wood rho: {v} kg/m^3 vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
