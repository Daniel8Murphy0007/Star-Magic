from fractions import Fraction as F
FT,SQ,Dp,N=F(1,10),F(57,100),4,9
v=N - FT*SQ + FT**2*Dp
t=8.988
print(f"k_e lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
