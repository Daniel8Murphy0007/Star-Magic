"""S693 (Tier KK) — Mercury year / 10 day.

Locked-primitive closure derived by _tier_KK_search.py.
expr = +beta +beta3 +3 +5
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = +(beta)+(beta**3)+Fraction(3,1)+Fraction(5,1)
pred = float(val)
obs  = 8.797
err  = abs(pred-obs)/obs*100
print(f"Mercury_year_per_10_day: {pred} vs {obs} -> {err:.4f}%")
