from fractions import Fraction as F
FT,Ph,SQ,K,Dp,DB,SO=F(1,10),F(5,6),F(57,100),F(25,12),4,6,10
v=FT*SO + FT*SQ*Ph/DB
t=1.008
print(f"H mass: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
