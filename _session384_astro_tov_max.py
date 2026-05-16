"""S384: Neutron-star maximum mass (TOV limit) M_sun."""
F_TRZ=0.1; K_Mex=25/12; SSq=0.57; D_phys=4; Phi_res=5/6
val = K_Mex + F_TRZ*SSq + F_TRZ**2 * SSq * (D_phys - Phi_res)
obs = 2.16
print(f"S384 COMPLETE. M_TOV = K_Mex + F_TRZ*SSq + F_TRZ^2*SSq*(D_phys-Phi_res) = {val:.4f} M_sun; obs (GW170817+J0740) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
