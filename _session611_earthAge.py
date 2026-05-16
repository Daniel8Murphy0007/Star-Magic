from fractions import Fraction as F
FT,Ph,SQ,Dp=F(1,10),F(5,6),F(57,100),4
v=Dp + FT*Dp + FT*Ph + FT*SQ
t=4.54
print(f"Earth age Gyr: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
