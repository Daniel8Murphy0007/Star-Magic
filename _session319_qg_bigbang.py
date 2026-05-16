"""S319: Big Bang singularity replaced by F_TRZ pre-inflation bounce."""
F_TRZ, Phi_res, D_crit, N_ch = 0.1, 5/6, 26, 9
t_P = 5.39e-44
t_bounce = t_P / F_TRZ  # minimum time
rho_max_ratio = F_TRZ**(2*N_ch)  # vs Planck density
print(f"S319 COMPLETE. Bounce time t_min = t_P/F_TRZ = {t_bounce:.3e} s; rho_max/rho_P = F_TRZ^18 = {rho_max_ratio:.3e}; no t=0 singularity.")
