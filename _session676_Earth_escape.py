"""S676 Earth escape velocity closure (Tier JJ).
Expr: F*K^5 + F^3*K^5 + beta^3 + 2 + 5
"""
from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
pred = F*K**5 + F**3*K**5 + beta**3 + 2 + 5
obs  = Fraction(11186,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth escape (km/s): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
