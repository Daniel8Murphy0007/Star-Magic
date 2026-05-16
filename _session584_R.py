from fractions import Fraction as F
FT,Ph,SQ,K,Dp,DB=F(1,10),F(5,6),F(57,100),F(25,12),4,6
v=K*(Dp - FT**2)
t=8.314
print(f"R gas: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
