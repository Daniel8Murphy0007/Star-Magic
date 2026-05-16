from fractions import Fraction as F
FT,K,Dp,N,SO,A=F(1,10),F(25,12),4,9,10,60
v=A*K + N + Dp - FT*SO + FT**2*Dp
t=137.036
print(f"alpha^-1: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
