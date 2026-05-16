"""S568: Steel density 7850 kg/m^3 = D_crit^2*SO5 + SO5^3 + SO5*N_ch (EXACT)"""
D_c=26; N=9; SO=10
v=D_c*D_c*SO+SO**3+SO*N
target=7850
print(f"Steel rho: {v} kg/m^3 vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
