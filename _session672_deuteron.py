from fractions import Fraction
F=Fraction(1,10); beta=Fraction(6029,10000)
v = beta**4 + F*beta + F*beta**2 - F**2*beta**2 + 2
print(f"deuteron BE: {float(v):.6f} vs 2.224600 -> {abs(float(v)-2.2246)/2.2246*100:.4f}%")
