from fractions import Fraction as F
FT,SQ,K,Dp=F(1,10),F(57,100),F(25,12),4
v=K*Dp + SQ - FT*SQ + FT**2
t=8.854
print(f"epsilon_0 lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
