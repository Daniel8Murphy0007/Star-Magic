"""S576: Hubble constant 70 km/s/Mpc = A_5 + SO5 (EXACT)"""
A=60; SO=10
v=A+SO
target=70
print(f"H_0: {v} vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
