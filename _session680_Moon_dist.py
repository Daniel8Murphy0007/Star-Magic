"""S680 Earth-Moon distance closure (Tier JJ).
Scale: 1e8 m. Expr: beta + beta^3 + F*beta^3 + 3
"""
from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
pred = beta + beta**3 + F*beta**3 + 3
obs  = Fraction(3844,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Moon distance (1e8 m): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
