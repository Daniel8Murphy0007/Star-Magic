"""S569: Aluminum density 2700 kg/m^3 = D_crit*SO5^2 + N_ch*SO5 + SO5 (EXACT)"""
D_c=26; N=9; SO=10
v=D_c*SO*SO+N*SO+SO
target=2700
print(f"Al rho: {v} kg/m^3 vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
