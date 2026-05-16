from fractions import Fraction as F
FT,Ph,SQ,K,Dp=F(1,10),F(5,6),F(57,100),F(25,12),4
v=Ph + SQ - FT**2*SQ*Dp
t=1.381
print(f"k_B lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
