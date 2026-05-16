"""S565: Concrete f'c = 30 MPa = D_crit + D_phys (EXACT)"""
D_c=26; D_p=4
v=D_c+D_p
target=30
print(f"Concrete fc: {v} MPa vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
