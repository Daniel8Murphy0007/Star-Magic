"""S682 Earth bulk density closure (Tier JJ).
Expr: beta^2 + beta^5 + F*beta + 5
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = beta**2 + beta**5 + F*beta + 5
obs  = Fraction(5514,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth density (g/cm^3): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
