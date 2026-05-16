from fractions import Fraction as F
FT,SQ,DB=F(1,10),F(57,100),6
v=DB + SQ + FT*SQ
t=6.626
print(f"h Planck lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
