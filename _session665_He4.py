from fractions import Fraction
F=Fraction(1,10); K=Fraction(25,12); beta=Fraction(6029,10000)
v = F*K**5 + beta**5 + F*beta + F**2*beta + 3
print(f"He4 BE/A: {float(v):.6f} vs 7.073900 -> {abs(float(v)-7.0739)/7.0739*100:.4f}%")
