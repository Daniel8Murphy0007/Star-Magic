"""S678 Earth orbital velocity closure (Tier JJ).
Expr: K^5 - F*K^5 - beta + F*beta - 5
"""
from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
pred = K**5 - F*K**5 - beta + F*beta - 5
obs  = Fraction(2978,100)
err  = abs(float(pred-obs)/float(obs))*100
print(f"Earth orbital v (km/s): {float(pred):.6f} vs {float(obs)} -> {err:.4f}%")
