from fractions import Fraction as F
FT,Ph,SQ,DB,SO,A=F(1,10),F(5,6),F(57,100),6,10,60
v=A*DB + SO + DB + Ph - FT*Ph - FT**2*SQ
t=376.73
print(f"Z_0 Ohm: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
