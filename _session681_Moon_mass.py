"""S681 Moon mass closure (Tier JJ).
Scale: 1e22 kg. Expr: beta^3 + beta^4 + 2 + 5
"""
from fractions import Fraction
beta=Fraction(6029,10000)
pred = beta**3 + beta**4 + 2 + 5
obs  = Fraction(7342,1000)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Moon mass (1e22 kg): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
