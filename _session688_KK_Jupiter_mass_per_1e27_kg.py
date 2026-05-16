"""S688 (Tier KK) — M_J / 10^27 kg.

Locked-primitive closure derived by _tier_KK_search.py.
expr = -F.beta -F.beta3 -F.beta4 +2
"""
from fractions import Fraction
F=Fraction(1,10); Phi=Fraction(5,6); SSq=Fraction(57,100); K=Fraction(25,12); beta=Fraction(6029,10000)
val = -(F*beta)-(F*beta**3)-(F*beta**4)+Fraction(2,1)
pred = float(val)
obs  = 1.898
err  = abs(pred-obs)/obs*100
print(f"Jupiter_mass_per_1e27_kg: {pred} vs {obs} -> {err:.4f}%")
