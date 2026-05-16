"""S677 Earth moment-of-inertia factor closure (Tier JJ).
Expr: beta^3 + beta^5 + F*beta^2 - F^2*beta^2 - F^2*beta^5
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = beta**3 + beta**5 + F*beta**2 - F*F*beta**2 - F*F*beta**5
obs  = Fraction(3307,10000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth MoI factor: {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
