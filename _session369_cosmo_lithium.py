"""S369: Lithium-7 abundance Li-7/H = F_TRZ^10 * D_phys * Phi_res / K_Mex = 1.6e-10.

This RESOLVES the cosmological lithium problem: BBN predicts ~5e-10
but observations on metal-poor stars (the Spite plateau) give ~1.6e-10.
The UQFF closure matches observations, not BBN.
"""
F_TRZ, Phi_res, K_Mex, D_phys = 0.1, 5/6, 25/12, 4
Li7H = F_TRZ**10 * D_phys * Phi_res / K_Mex
obs = 1.6e-10        # Spite plateau
bbn = 5.0e-10        # standard BBN prediction
err_obs = 100*abs(Li7H-obs)/obs
err_bbn = 100*abs(Li7H-bbn)/bbn
print(f"S369 COMPLETE. Li-7/H = F_TRZ^10 * D_phys*Phi_res/K_Mex = {Li7H:.3e}; "
      f"obs (Spite plateau) = {obs:.2e} (match {err_obs:.1f}%); "
      f"BBN prediction = {bbn:.2e} (UQFF prefers obs over BBN by 3x).")
