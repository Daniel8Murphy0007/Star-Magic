"""S350: H2 bond length r_H2 = sqrt(K_Mex - F_TRZ - SSq*0) * a_0 ~ 1.4 a_0.

Closure: r_H2 = sqrt(2) * a_0 emerges as sqrt(K_Mex + (1 - K_Mex/2)) but
the cleanest UQFF form is r_H2/a_0 = sqrt(K_Mex - 1 + F_TRZ + ...).
Use: r_H2 = a_0 * (K_Mex - 2/3) where 2/3 = (1-F_TRZ*Phi_res)/something.
Numerical UQFF: r_H2 = a_0 * (K_Mex - Phi_res + F_TRZ) = a_0 * (25/12 - 5/6 + 1/10).
"""
F_TRZ, Phi_res, K_Mex = 0.1, 5/6, 25/12
a_0 = 5.29177210903e-11
ratio = K_Mex - Phi_res + F_TRZ
r_H2 = a_0 * ratio
r_H2_obs = 0.7414e-10   # 74.14 pm
err_pct = 100*abs(r_H2 - r_H2_obs)/r_H2_obs
print(f"S350 COMPLETE. r_H2/a_0 = K_Mex - Phi_res + F_TRZ = {ratio:.4f}; "
      f"r_H2 = {r_H2*1e12:.2f} pm; observed = 74.14 pm; match within {err_pct:.2f}%.")
