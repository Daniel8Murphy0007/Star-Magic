"""S692 (Tier KK) — Mars semi-major axis (AU).

Locked-primitive closure derived by _tier_KK_search.py.
expr = -beta2 -beta5 -F.beta2 +2
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = -(beta**2)-(beta**5)-(F*beta**2)+Fraction(2,1)
pred = float(val)
obs  = 1.524
err  = abs(pred-obs)/obs*100
print(f"Mars_orbit_AU: {pred} vs {obs} -> {err:.4f}%")
