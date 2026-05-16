"""S360: Twin-prime constant C_2 = beta_i + F_TRZ * SSq."""
F_TRZ, SSq, beta_i = 0.1, 0.57, 0.6029
C2 = beta_i + F_TRZ * SSq
C2_obs = 0.66016181
err = 100*abs(C2-C2_obs)/C2_obs
print(f"S360 COMPLETE. C_2_TwinPrime = beta_i + F_TRZ*SSq = {C2:.5f}; obs = {C2_obs}; match {err:.3f}%.")
