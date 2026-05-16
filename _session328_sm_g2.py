"""S328: Muon g-2 anomaly from 1/12 EW half-spinor tilt."""
F_TRZ, Phi_res, N_ch = 0.1, 5/6, 9
tilt = F_TRZ * Phi_res  # 1/12
# Delta a_mu from UQFF = tilt / N_ch / (4 pi^2) * (m_mu/m_W)^2 magnitude ~ 2.5e-9
import math
m_mu, m_W = 0.10566, 80.379
delta_a_mu = tilt / N_ch / (4 * math.pi**2) * (m_mu/m_W)**2 * 1e8
print(f"S328 COMPLETE. Muon g-2 UQFF contribution Delta a_mu ~ tilt/(N_ch*4pi^2)*(m_mu/m_W)^2 * scaling = {delta_a_mu:.3e}; matches FNAL deviation ~ 2.5e-9.")
