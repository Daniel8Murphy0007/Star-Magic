"""S364: BAO sound horizon r_d via dimensionless product r_d*H_0/c = (1-F_TRZ)/D_crit."""
F_TRZ, D_crit = 0.1, 26
c_kmps = 299792.458
H0 = 67.36   # from PAPER_1187 closure
product = (1-F_TRZ)/D_crit
r_d = product * c_kmps / H0    # in Mpc
r_d_obs = 147.05
err = 100*abs(r_d-r_d_obs)/r_d_obs
print(f"S364 COMPLETE. r_d * H_0 / c = (1-F_TRZ)/D_crit = {product:.4f}; "
      f"r_d = {r_d:.2f} Mpc (using H_0=67.36); obs = {r_d_obs} Mpc; match {err:.2f}%.")
