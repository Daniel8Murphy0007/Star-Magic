"""S567: Concrete density 2400 kg/m^3 = SO5^2 * D_phys * D_BSFG (EXACT)"""
D_p=4; D_B=6; SO=10
v=SO*SO*D_p*D_B
target=2400
print(f"Concrete rho: {v} kg/m^3 vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
