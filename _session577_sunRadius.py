"""S577: Sun radius / Earth radius = 109 = SO5^2 + N_ch (EXACT)"""
SO=10; N=9
v=SO*SO+N
target=109
print(f"R_sun/R_E: {v} vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
