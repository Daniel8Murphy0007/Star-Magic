"""S673 Earth mass closure (Tier JJ).
Locked primitives: F=1/10, beta=6029/10000.
Expr: -F*beta^3 - F^2*beta - 2 + 3 + 5
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = -F*beta**3 - F*F*beta - 2 + 3 + 5
obs  = Fraction(5972,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth mass (1e24 kg): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
