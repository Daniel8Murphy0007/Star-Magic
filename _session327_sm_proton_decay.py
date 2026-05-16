"""S327: Proton lifetime from F_TRZ^(N_ch*D_phys) decay suppression."""
F_TRZ, N_ch, D_phys = 0.1, 9, 4
t_P = 5.39e-44
tau_p = t_P / F_TRZ**(N_ch*D_phys)  # F_TRZ^36 = 10^-36 suppression
print(f"S327 COMPLETE. Proton lifetime tau_p = t_P / F_TRZ^(N_ch*D_phys) = t_P / F_TRZ^36 = {tau_p:.3e} s (obs > 1.6e34 yr = 5e41 s).")
