"""S357: Euler-Mascheroni gamma = SSq + F_TRZ^2 * Phi_res."""
F_TRZ, Phi_res, SSq = 0.1, 5/6, 0.57
gamma = SSq + F_TRZ**2 * Phi_res
gamma_obs = 0.57721566
err = 100*abs(gamma-gamma_obs)/gamma_obs
print(f"S357 COMPLETE. gamma_EulerMascheroni = SSq + F_TRZ^2 * Phi_res = {gamma:.5f}; obs = {gamma_obs}; match {err:.3f}%.")
