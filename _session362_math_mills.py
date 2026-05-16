"""S362: Mills constant theta = 1 + K_Mex * beta_i / D_phys."""
K_Mex, beta_i, D_phys = 25/12, 0.6029, 4
theta = 1 + K_Mex * beta_i / D_phys
theta_obs = 1.30637788
err = 100*abs(theta-theta_obs)/theta_obs
print(f"S362 COMPLETE. theta_Mills = 1 + K_Mex*beta_i/D_phys = {theta:.5f}; obs = {theta_obs}; match {err:.3f}%.")
