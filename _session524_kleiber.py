from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6)
v = Ph - F*Ph; target = Fraction(3,4)
print(f"S524 COMPLETE. Kleiber allometric exponent = {float(v):.4f} = Phi_res-F_TRZ*Phi_res = 3/4; target {float(target)}; match {abs(float(v)-float(target))/float(target)*100:.4f}%.")
