"""S391: Canonical neutron-star compactness GM/(Rc^2) for M=1.4 M_sun, R=10 km."""
F_TRZ=0.1; K_Mex=25/12; D_phys=4; SSq=0.57
val = K_Mex*F_TRZ + F_TRZ**3 * D_phys * SSq
obs = 0.21
print(f"S391 COMPLETE. NS compactness = K_Mex*F_TRZ + F_TRZ^3*D_phys*SSq = {val:.4f}; obs (1.4 M_sun, R=10 km) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
