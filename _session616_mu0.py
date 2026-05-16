from fractions import Fraction as F
FT,Ph,SQ,K=F(1,10),F(5,6),F(57,100),F(25,12)
v=K - Ph + FT**2*SQ
t=1.257
print(f"mu_0 lead: {float(v):.6f} vs {t} -> {abs(float(v)-t)/t*100:.4f}%")
