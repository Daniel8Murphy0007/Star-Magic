from fractions import Fraction as F
FT,SQ,K,Dp=F(1,10),F(57,100),F(25,12),4
v=K + FT*Dp - FT*SQ
t=2.426
print(f"Compton wavelength lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
