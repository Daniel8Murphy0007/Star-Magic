from fractions import Fraction as F
FT,Ph,SQ,Dp,N,SO=F(1,10),F(5,6),F(57,100),4,9,10
v=N + SO + SO + Ph - FT**2*Dp - FT**2*SQ
t=29.78
print(f"Earth orbital v km/s: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
