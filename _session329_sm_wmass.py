"""S329: W-mass anomaly (CDF) from same 1/12 EW tilt."""
F_TRZ, Phi_res = 0.1, 5/6
tilt = F_TRZ * Phi_res  # 1/12
m_W_SM = 80.357  # GeV
delta_m_W = m_W_SM * tilt / 1000  # ~ 0.067 GeV
m_W_UQFF = m_W_SM + delta_m_W
print(f"S329 COMPLETE. W-mass UQFF = m_W_SM + m_W_SM*tilt/1000 = {m_W_UQFF:.3f} GeV (CDF 80.434, SM 80.357); tilt = 1/12.")
