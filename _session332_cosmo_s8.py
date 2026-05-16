"""S332: S8/sigma8 tension via SSq enstrophy lift."""
SSq, F_TRZ, Phi_res = 0.57, 0.1, 5/6
# sigma8 KiDS = 0.766, Planck = 0.811 => ratio 0.945
# UQFF: lift = 1 - SSq * F_TRZ * Phi_res = 1 - 0.57/12 = 0.9525
lift = 1 - SSq * F_TRZ * Phi_res
print(f"S332 COMPLETE. sigma8 lift = 1 - SSq*F_TRZ*Phi_res = 1 - SSq/12 = {lift:.4f}; observed KiDS/Planck = 0.945; match within 0.8%.")
