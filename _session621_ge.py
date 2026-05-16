from fractions import Fraction as F
FT,SQ,K,Dp=F(1,10),F(57,100),F(25,12),4
v=K - FT*SQ - FT**2*Dp + FT**2*SQ + FT**2*K - FT**2
t=2.0023
print(f"g_e factor: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
