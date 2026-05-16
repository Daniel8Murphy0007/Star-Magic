"""S574: AU/Earth-radii = 23481 = D_c*N_ch*SO5^2 + A_5 + D_c - D_p + F_TRZ*SO5 + F_TRZ*Phi_res - K_Mex (EXACT)
       23400 + 60 + 26 - 4 + 1 + 1/12 - 25/12 = 23483 + (1-25)/12 = 23483 - 2 = 23481"""
from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); D_p=4; D_c=26; N=9; SO=10; A=60
v=D_c*N*SO*SO+A+D_c-D_p+F*SO+F*Ph-K
target=23481
print(f"AU/RE: {float(v)} vs {target} -> {abs(float(v)-target)/target*100:.4f}% [EXACT]" if v==target else f"{float(v)}")
