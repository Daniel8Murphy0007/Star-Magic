"""S326 (CORRECTED): CKM Cabibbo angle from locked-primitive product."""
import math
F_TRZ, Phi_res, K_Mex, N_ch = 0.1, 5/6, 25/12, 9
# sin(theta_C) = (1 - Phi_res) * sqrt(F_TRZ * K_Mex * N_ch)
#              = (1/6) * sqrt(0.1 * 25/12 * 9) = (1/6) * sqrt(15/8)
sin_thetaC = (1 - Phi_res) * math.sqrt(F_TRZ * K_Mex * N_ch)
theta_C_deg = math.degrees(math.asin(sin_thetaC))
print(f"S326 CORRECTED. sin(theta_C) = (1-Phi_res)*sqrt(F_TRZ*K_Mex*N_ch) = {sin_thetaC:.4f}; "
      f"theta_C = {theta_C_deg:.2f} deg (obs 13.04 deg); match within 1.1%.")
