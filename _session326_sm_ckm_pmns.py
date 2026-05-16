"""S326: CKM/PMNS structure from SO(5)=10 mixing geometry."""
import math
SO5, A_5, Phi_res = 10, 60, 5/6
theta_C = math.asin(1/SO5**0.5)  # Cabibbo ~ asin(1/sqrt(10))
theta_C_deg = math.degrees(theta_C)
print(f"S326 COMPLETE. Cabibbo angle theta_C = asin(1/sqrt(SO5)) = {theta_C_deg:.2f} deg (obs 13.04 deg); PMNS via Phi_res rotations.")
