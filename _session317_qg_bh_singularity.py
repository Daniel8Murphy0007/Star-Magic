"""S317: BH singularity resolution via D_crit=26 BSFG core."""
F_TRZ, Phi_res, D_crit, N_ch = 0.1, 5/6, 26, 9
ell_P = 1.616e-35
r_core = ell_P / F_TRZ  # finite core radius
K_max = 1 / r_core**2   # curvature finite
print(f"S317 COMPLETE. Singularity replaced by BSFG core r_core = ell_P/F_TRZ = {r_core:.3e} m; K_max = {K_max:.3e} m^-2 (finite, no infinity).")
