"""S557: Ocean salinity 35 ppt = D_crit + N_ch (EXACT)"""
D_c=26; N=9
v=D_c+N
target=35
err=abs(v-target)/target*100
print(f"Salinity: {v} ppt vs {target} -> {err:.4f}% [EXACT]" if v==target else f"{v} vs {target} -> {err:.4f}%")
