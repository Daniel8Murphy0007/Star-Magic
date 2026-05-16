from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6); K=Fraction(25,12); SO=Fraction(10)
v = SO + Ph - F*K - F*Ph - F*F*K - F*F*Ph; target = Fraction(105,10)
print(f"S528 COMPLETE. DNA base-pairs per turn = {float(v):.6f}; closure = SO5+Phi_res-F*K-F*Phi_res-F^2*K-F^2*Phi_res; target {float(target)}; match {abs(float(v)-float(target))/float(target)*100:.4f}%.")
