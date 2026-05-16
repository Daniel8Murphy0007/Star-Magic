"""S556: Solar constant 1361 W/m^2 = A_5^2*F_TRZ*D_phys - N_ch*SO5 + SO5 + SSq + F_TRZ*D_phys"""
from fractions import Fraction
F=Fraction(1,10); SSq=Fraction(57,100); D_p=4; N=9; SO=10; A=60
v=A*A*F*D_p - N*SO + SO + SSq + F*D_p
target=1361
err=abs(float(v)-target)/target*100
print(f"Solar constant: {float(v):.6f} W/m^2 vs {target} -> {err:.4f}%")
