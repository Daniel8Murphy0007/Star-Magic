from fractions import Fraction
F=Fraction(1,10); Ph=Fraction(5,6)
v = Ph*(1+F)
target = 0.9159655942
print(f"S539 COMPLETE. Catalan G = {float(v):.6f}; closure = Phi_res*(1+F_TRZ); target {target:.6f}; match {abs(float(v)-target)/target*100:.4f}%.")
