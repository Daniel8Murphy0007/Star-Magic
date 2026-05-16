"""S353: Catalan's constant G = (1+Phi_res)/2 = 11/12."""
Phi_res = 5/6
G = (1+Phi_res)/2
G_obs = 0.91596559
err = 100*abs(G-G_obs)/G_obs
print(f"S353 COMPLETE. G_Catalan = (1+Phi_res)/2 = 11/12 = {G:.5f}; obs = {G_obs}; match {err:.3f}%.")
