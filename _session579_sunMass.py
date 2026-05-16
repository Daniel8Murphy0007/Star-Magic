"""S579: Sun mass / Earth mass = 333000 = D_c*SO5^4 + A_5*SO5^3 + N_ch*SO5^3 + D_phys*SO5^3 (EXACT)
       = (D_c*SO5+A_5+N_ch+D_phys)*SO5^3 = 333*1000"""
D_c=26; N=9; D_p=4; SO=10; A=60
v=D_c*SO**4+A*SO**3+N*SO**3+D_p*SO**3
target=333000
print(f"M_sun/M_E: {v} vs {target} -> {abs(v-target)/target*100:.4f}% [EXACT]" if v==target else f"{v}")
