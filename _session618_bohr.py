from fractions import Fraction as F
FT,Ph,SQ,Dp=F(1,10),F(5,6),F(57,100),4
v=Dp + Ph + FT*Dp + FT*SQ
t=5.292
print(f"Bohr radius lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
