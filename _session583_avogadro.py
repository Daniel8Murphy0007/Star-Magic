from fractions import Fraction as F
FT,Ph,SQ,K,Dp,DB,Dc,N,SO,A=F(1,10),F(5,6),F(57,100),F(25,12),4,6,26,9,10,60
v=DB + FT**2*SQ*Dp
t=6.022
print(f"Avogadro lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
