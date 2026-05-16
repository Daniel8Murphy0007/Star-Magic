"""S679 Earth sidereal year closure (Tier JJ).
Scale: 1e2 days. Expr: beta + F*beta^2 + F^2*beta + 3
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = beta + F*beta**2 + F*F*beta + 3
obs  = Fraction(36525,10000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth year (1e2 days): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
