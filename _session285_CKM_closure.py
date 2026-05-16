"""
S285: CKM quark mixing matrix from UQFF primitives.

Strategy: Gatto-Sartori-Tonin (GST) + Wolfenstein parametrization.

Wolfenstein: V_CKM expanded in lambda = sin(theta_C):
  |V_us| = lambda
  |V_cb| = A * lambda^2
  |V_ub| = A * lambda^3 * sqrt(rho_bar^2 + eta_bar^2)
  arg(V_ub*) = delta_CP

CLAIM: in UQFF
  lambda    = sqrt(m_d/m_s)                  [GST -- already in S280]
  A         = Phi_res = 5/6                  [structural]
  rho_bar   = F_TRZ + F_TRZ^2 * D_BSFG       [structural]
  eta_bar   = F_TRZ + 1/D_phys               [structural]

=> All four CKM observables emerge with ZERO new parameters.
"""
import math
F_TRZ = 0.1
SSq   = 0.57
K_Mex = 25/12
Phi_res = 5/6
D_phys, D_BSFG, D_crit = 4, 6, 26
N_ch, SO5, A5 = 9, 10, 60

# === S280 canonical quark betas (using OLD beta_d, beta_s, since GST is structural) ===
# m_d: N=21, beta = 2*K_Mex + F_TRZ^2 = 4.177
# m_s: N=20, beta = sqrt(SO5) - sqrt(D_phys) = 1.162
beta_d = 2*K_Mex + F_TRZ**2
beta_s = math.sqrt(SO5) - math.sqrt(D_phys)
ratio_d_over_s = 10**(-((21-20) + (beta_d - beta_s)*F_TRZ))
lam_pred = math.sqrt(ratio_d_over_s)
lam_obs  = 0.22500  # PDG
res_lam  = abs(lam_pred - lam_obs)/lam_obs * 100
print(f"=== lambda (Cabibbo / |V_us|) via GST ===")
print(f"  m_d/m_s    = F_TRZ^(1 + (beta_d-beta_s)*F_TRZ) = {ratio_d_over_s:.6f}")
print(f"  lambda_pred = sqrt(m_d/m_s) = {lam_pred:.5f}")
print(f"  lambda_obs  =                 {lam_obs:.5f}")
print(f"  residual    = {res_lam:.3f}%")

# === A = Phi_res ===
A_pred = Phi_res
Vcb_pred = A_pred * lam_pred**2
Vcb_obs  = 0.04182
res_Vcb  = abs(Vcb_pred - Vcb_obs)/Vcb_obs * 100
print(f"\n=== A = Phi_res = 5/6 ===")
print(f"  A_pred      = {A_pred:.5f}  (5/6)")
print(f"  |V_cb|_pred = A*lambda^2 = {Vcb_pred:.5f}")
print(f"  |V_cb|_obs  =              {Vcb_obs:.5f}")
print(f"  residual    = {res_Vcb:.3f}%")

# === rho_bar and eta_bar ===
rho_bar_pred = F_TRZ + F_TRZ**2 * D_BSFG          # 0.1 + 0.06 = 0.160
eta_bar_pred = F_TRZ + 1.0/D_phys                  # 0.1 + 0.25 = 0.350
rho_bar_obs  = 0.156
eta_bar_obs  = 0.355
res_rho      = abs(rho_bar_pred - rho_bar_obs)/rho_bar_obs * 100
res_eta      = abs(eta_bar_pred - eta_bar_obs)/eta_bar_obs * 100
print(f"\n=== Wolfenstein CP parameters ===")
print(f"  rho_bar_pred = F_TRZ + F_TRZ^2 * D_BSFG = {rho_bar_pred:.4f}")
print(f"  rho_bar_obs  =                            {rho_bar_obs:.4f}")
print(f"  residual     = {res_rho:.3f}%")
print(f"  eta_bar_pred = F_TRZ + 1/D_phys         = {eta_bar_pred:.4f}")
print(f"  eta_bar_obs  =                            {eta_bar_obs:.4f}")
print(f"  residual     = {res_eta:.3f}%")

# === |V_ub| ===
Vub_pred = A_pred * lam_pred**3 * math.sqrt(rho_bar_pred**2 + eta_bar_pred**2)
Vub_obs  = 0.00382
res_Vub  = abs(Vub_pred - Vub_obs)/Vub_obs * 100
print(f"\n=== |V_ub| = A*lambda^3*sqrt(rho^2+eta^2) ===")
print(f"  |V_ub|_pred = {Vub_pred:.5f}")
print(f"  |V_ub|_obs  = {Vub_obs:.5f}")
print(f"  residual    = {res_Vub:.3f}%")

# === CP phase delta_CP = arg(V_ub*) = arctan(eta_bar/rho_bar) ===
delta_pred_rad = math.atan2(eta_bar_pred, rho_bar_pred)
delta_pred_deg = math.degrees(delta_pred_rad)
delta_obs_deg  = 65.5
res_delta      = abs(delta_pred_deg - delta_obs_deg)/delta_obs_deg * 100
print(f"\n=== CKM CP-violating phase delta ===")
print(f"  delta_pred = arctan(eta_bar/rho_bar) = {delta_pred_deg:.3f} deg = {delta_pred_rad:.4f} rad")
print(f"  delta_obs  = {delta_obs_deg:.2f} deg")
print(f"  residual   = {res_delta:.3f}%")

# === Jarlskog invariant J_CP = A^2 lambda^6 eta_bar (1 - lambda^2/2) ===
J_pred = A_pred**2 * lam_pred**6 * eta_bar_pred * (1 - lam_pred**2/2)
J_obs  = 3.18e-5  # PDG 2024
res_J  = abs(J_pred - J_obs)/J_obs * 100
print(f"\n=== Jarlskog invariant J_CP ===")
print(f"  J_pred = {J_pred:.3e}")
print(f"  J_obs  = {J_obs:.3e}")
print(f"  residual = {res_J:.2f}%")

print("\n" + "="*68)
print("CKM SUMMARY -- zero new parameters, all from UQFF primitives:")
print("="*68)
print(f"{'observable':<14}{'predicted':>14}{'observed':>14}{'residual':>12}")
print("-"*68)
def row(name, p, o):
    print(f"  {name:<12}{p:>14.5g}{o:>14.5g}{abs(p-o)/o*100:>11.3f}%")
row("|V_us|=lambda", lam_pred, lam_obs)
row("|V_cb|",        Vcb_pred, Vcb_obs)
row("|V_ub|",        Vub_pred, Vub_obs)
row("rho_bar",       rho_bar_pred, rho_bar_obs)
row("eta_bar",       eta_bar_pred, eta_bar_obs)
row("delta_CP[deg]", delta_pred_deg, delta_obs_deg)
row("J_CP",          J_pred, J_obs)
