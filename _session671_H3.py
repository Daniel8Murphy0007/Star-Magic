from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
v = -beta**5 - F*beta - F*beta**2 + F**2*beta**3 + 3
print(f"H3 BE/A: {float(v):.6f} vs 2.827000 -> {abs(float(v)-2.827)/2.827*100:.4f}%")
