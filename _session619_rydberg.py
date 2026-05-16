from fractions import Fraction as F
FT,SQ,Dp,SO=F(1,10),F(57,100),4,10
v=FT*SO + FT*SQ + FT**2*Dp
t=1.0974
print(f"Rydberg lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
