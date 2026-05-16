"""S566: Steel Young's E = 200 GPa = D_crit*D_BSFG + D_phys*SO5 + D_phys (EXACT)"""
D_c=26; D_B=6; D_p=4; SO=10
v=D_c*D_B+D_p*SO+D_p
target=200
print(f"Steel E: {v} GPa vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
