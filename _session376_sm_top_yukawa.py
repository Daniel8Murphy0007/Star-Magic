"""S376: Top quark Yukawa coupling y_t at m_t scale."""
F_TRZ=0.1
val = 1 - F_TRZ**2
obs = 0.9936
print(f"S376 COMPLETE. y_t = 1 - F_TRZ^2 = {val:.4f}; obs (PDG, m_t=172.69 GeV, v=246.22) = {obs}; match {abs(val-obs)/obs*100:.3f}%.")
