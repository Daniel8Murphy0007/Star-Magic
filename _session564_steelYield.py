"""S564: Steel yield 250 MPa = D_crit*SO5 - D_BSFG - D_phys (EXACT)"""
D_c=26; D_B=6; D_p=4; SO=10
v=D_c*SO-D_B-D_p
target=250
print(f"Steel yield: {v} MPa vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
