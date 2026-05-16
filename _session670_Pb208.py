from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 + beta + beta**2 - F*beta**3 + 3
print(f"Pb208 BE/A: {float(v):.6f} vs 7.867500 -> {abs(float(v)-7.8675)/7.8675*100:.4f}%")
