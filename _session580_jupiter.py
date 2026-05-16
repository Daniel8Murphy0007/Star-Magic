"""S580: Jupiter mass / Earth mass = 317.8 ; D_c*SO5 + SSq*SO5 + SO5*D_p + SO5 + K_Mex"""
from fractions import Fraction
K=Fraction(25,12); SSq=Fraction(57,100); D_c=26; D_p=4; SO=10
v=D_c*SO+SSq*SO+SO*D_p+SO+K
target=317.8
print(f"M_J/M_E: {float(v):.6f} vs {target} -> {abs(float(v)-target)/target*100:.4f}%")
