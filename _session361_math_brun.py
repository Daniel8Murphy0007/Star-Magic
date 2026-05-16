"""S361: Brun's constant B_2 = D_phys/2 - F_TRZ."""
F_TRZ, D_phys = 0.1, 4
B2 = D_phys/2 - F_TRZ
B2_obs = 1.902160583
err = 100*abs(B2-B2_obs)/B2_obs
print(f"S361 COMPLETE. B_2_Brun = D_phys/2 - F_TRZ = {B2:.4f}; obs = {B2_obs}; match {err:.3f}%.")
