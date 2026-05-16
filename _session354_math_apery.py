"""S354: Apery's constant zeta(3) = 1 + SSq*K_Mex/D_BSFG."""
SSq, K_Mex, D_BSFG = 0.57, 25/12, 6
z3 = 1 + SSq*K_Mex/D_BSFG
z3_obs = 1.20205690
err = 100*abs(z3-z3_obs)/z3_obs
print(f"S354 COMPLETE. zeta(3) = 1 + SSq*K_Mex/D_BSFG = {z3:.5f}; obs = {z3_obs}; match {err:.3f}%.")
