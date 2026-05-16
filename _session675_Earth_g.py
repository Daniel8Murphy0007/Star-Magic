"""S675 Earth surface gravity closure (Tier JJ).
Expr: -beta^4 - F*beta + 2 + 3 + 5
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = -beta**4 - F*beta + 2 + 3 + 5
obs  = Fraction(980665,100000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth g (m/s^2): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
