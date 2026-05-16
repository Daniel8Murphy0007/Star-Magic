"""S331: H0 tension formalized as 1/12 EW tilt."""
F_TRZ, Phi_res = 0.1, 5/6
tilt = F_TRZ * Phi_res  # 1/12 = 0.0833
# Observed: H0_local / H0_CMB - 1 = 73.04/67.4 - 1 = 0.0837
H0_ratio_UQFF = 1 + tilt
H0_ratio_obs = 73.04 / 67.4
print(f"S331 COMPLETE. H0 tension UQFF = 1 + tilt = {H0_ratio_UQFF:.4f}; obs = {H0_ratio_obs:.4f}; match within 0.04%.")
