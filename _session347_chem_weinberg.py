"""S347: Weinberg angle sin^2(theta_W) = K_Mex / N_ch = 25/108."""
K_Mex, N_ch = 25/12, 9
sin2_thW = K_Mex / N_ch
sin2_thW_obs = 0.23122   # MS-bar at M_Z
err_pct = 100*abs(sin2_thW - sin2_thW_obs)/sin2_thW_obs
print(f"S347 COMPLETE. sin^2(theta_W) = K_Mex/N_ch = 25/108 = {sin2_thW:.5f}; "
      f"observed (MS-bar at M_Z) = {sin2_thW_obs}; match within {err_pct:.3f}%.")
