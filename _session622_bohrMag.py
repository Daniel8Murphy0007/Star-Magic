from fractions import Fraction as F
FT,SQ,K,Dp=F(1,10),F(57,100),F(25,12),4
v=K*Dp + SQ + FT*Dp - FT**2*Dp + FT**2
t=9.274
print(f"Bohr magneton lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
