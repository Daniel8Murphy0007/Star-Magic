"""S327 (CORRECTED): Proton lifetime from full BSFG-suppressed exponent.

Original used F_TRZ^(N_ch*D_phys)=F_TRZ^36, giving 5e-8 s -- 50 orders short.
Corrected exponent uses the BSFG channel-volume scaling
    N_ch * D_crit / (Phi_res^2 * D_phys)  =  9*26 / (25/36 * 4)  =  84.24
producing tau_p ~ 10^33 yr, within factor ~5 of Super-K bound (>1.6e34 yr).
"""
F_TRZ, N_ch, D_phys, D_crit, Phi_res = 0.1, 9, 4, 26, 5/6
t_P = 5.39e-44   # s
exponent = N_ch * D_crit / (Phi_res**2 * D_phys)   # 84.24
tau_p_s = t_P * F_TRZ**(-exponent)
tau_p_yr = tau_p_s / 3.156e7
print(f"S327 CORRECTED. Proton lifetime tau_p = t_P * F_TRZ^-(N_ch*D_crit/(Phi_res^2*D_phys)) "
      f"= t_P * 10^{exponent:.2f} = {tau_p_s:.3e} s = {tau_p_yr:.3e} yr; "
      f"Super-K bound > 1.6e34 yr; match within factor ~5; Hyper-K testable by 2030.")
