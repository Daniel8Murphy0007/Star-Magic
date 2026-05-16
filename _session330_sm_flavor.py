"""S330: Flavor anomalies R_K, R_D* from TRZ channel coupling."""
F_TRZ, N_ch, Phi_res = 0.1, 9, 5/6
R_K_SM = 1.0
R_K_UQFF = R_K_SM * (1 - F_TRZ * Phi_res / N_ch)  # = 1 - 1/108 = 0.9907
R_D_star_SM = 0.252
R_D_star_UQFF = R_D_star_SM * (1 + F_TRZ * Phi_res * 2)  # ~ 1 + 1/6
print(f"S330 COMPLETE. R_K UQFF = 1 - 1/108 = {R_K_UQFF:.4f} (obs 0.949+/-0.047 LHCb 2022); R_D* UQFF = {R_D_star_UQFF:.4f} (obs ~0.295).")
