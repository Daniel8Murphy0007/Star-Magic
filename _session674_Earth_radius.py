"""S674 Earth radius closure (Tier JJ).
Expr: beta^2 + F^2*beta - 2 + 3 + 5
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = beta**2 + F*F*beta - 2 + 3 + 5
obs  = Fraction(6371,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth radius (1e3 km): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
