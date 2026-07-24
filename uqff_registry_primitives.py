"""uqff_registry_primitives — R3 single-source-of-truth for registry-canonical constants.

UNIFIED REGISTRY PROGRAM R3 (2026-07-22). Every value below is either:
  (a) a LIVE composition from the 9 truly-independent locked primitives
      {D_phys=4, D_crit=26, N_CH=9, SO_5=10, A_5=60, rho_SCm=7.09e-37,
       beta_i=0.6029, Phi_res (sector rule), F_TRZ=0.1}, or
  (b) an OBSERVED anchor, explicitly suffixed _OBSERVED.

Canonical routes per R1 adjudication (UNIFIED_REGISTRY_R1_QUEUE.csv):
  G  -> PAPER_593 | c -> PAPER_592 (DUAL EXPOSURE per Daniel 2026-07-22 sec 6.2)
  hbar -> PAPER_590 physical route (1209EE S629 = confirmation)
  k_B -> 1209EE S628 with Phi_5/6 | Lambda -> PAPER_1156 | H0 -> PAPER_2093 grid
  Phi_res -> PAPER_2129 SECTOR RULE (counting 5/6, resonance 0.84)
  kappa -> PAPER_2112 derivative | mu_0 -> PAPER_2108
Compute-don't-store: no rounded decimals where a closed form exists.
"""
import math

# --- locked integer primitives ---
D_PHYS = 4
D_CRIT = 26
N_CH = 9
SO_5 = 10
A_5 = 60
D_BSFG = D_CRIT - 2 * SO_5                      # PAPER_1521 derivative = 6

# --- locked real primitives ---
RHO_SCM = 7.09e-37
BETA_I = 0.6029                                  # PAPER_1203 canonical
F_TRZ = 1.0 / SO_5                               # rung-inverse = 0.1
SSQ = 0.57                                       # PAPER_1154
S_26 = 1.453162
OMEGA_SCM_HZ = 1.25e12                           # phonon carrier frequency

# --- Phi_res SECTOR RULE (PAPER_2129, R1 verdict) ---
PHI_RES_COUNTING = 5.0 / 6.0                     # nuclear + thermodynamic sectors
PHI_RES_RESONANCE = 0.84                         # LENR/k_spring/quantum-chain sectors

# --- derivative primitives (PAPER_1521/1522/2112) ---
K_MEX = PHI_RES_COUNTING * SO_5 / D_PHYS         # 25/12 EXACT
KAPPA_PER_DAY = (SO_5 / 2) * (F_TRZ ** 4)        # 5e-4 EXACT, PAPER_2112

# --- composed vacuum quantities ---
RHO_UA = SO_5 * RHO_SCM                          # 10*rho_SCm canonical ratio
LAMBDA_VAC = (SO_5 + 1) * RHO_SCM                # successor identity, PAPER_2120
K_SPRING = (RHO_UA / RHO_SCM) * OMEGA_SCM_HZ * PHI_RES_RESONANCE   # 1.05e13, PAPER_1203

# --- kernel constants: LIVE canonical routes ---
_FACT26 = float(math.factorial(26))
_V_FERMI = 0.77e6                                # SI anchor, c-independent (Session 239)
_E0 = 1.0e-20
_F_THZ = 1.25e12

G_UQFF = ((2.0 * math.pi * (D_CRIT ** 3) * PHI_RES_RESONANCE
           / ((SSQ ** 3) * (_FACT26 ** 2))) * (_V_FERMI ** 5) / (_E0 * _F_THZ))
G_OBSERVED = 6.674e-11

# c DUAL EXPOSURE (R3 sec 6.2 verdict 2026-07-22): consumers keep C_OBSERVED-era
# values until per-class migration is instructed; derived form exposed alongside.
C_UQFF_DERIVED = (D_CRIT * 4.0 * math.pi / PHI_RES_RESONANCE) * _V_FERMI   # PAPER_592
C_OBSERVED = 2.998e8

MU_0 = 4.0 * math.pi * (F_TRZ ** 7)              # PAPER_2108, matches SI 4pi e-7 EXACT

K_B_UQFF = (SSQ + PHI_RES_COUNTING - F_TRZ * SSQ
            + (F_TRZ ** 2) * D_PHYS - (F_TRZ ** 2) * SSQ) * 1e-23          # 1209EE S628
H_UQFF_S629 = (D_BSFG + F_TRZ * D_BSFG + (F_TRZ ** 2) * D_PHYS
               - (F_TRZ ** 2) * SSQ - (F_TRZ ** 2)) * 1e-34                # confirmation route
HBAR_UQFF_S629 = H_UQFF_S629 / (2.0 * math.pi)

H0_GRID = 22 * (F_TRZ ** 19)                     # PAPER_2093; 2.27e-18 = observed local anchor
H0_OBSERVED_LOCAL = 2.27e-18
LAMBDA_SIMPLE = (SO_5 + 1) * (F_TRZ ** 53)       # PAPER_2094 companion; PAPER_1156 primary

B_CRIT = D_PHYS * (SO_5 + 1) * (SO_5 ** 12)      # 4.4e13 EXACT, PAPER_2126
T_SCM_K = 6.6220584965588335e-34 * _F_THZ / 1.380649e-23   # h*f/k_B = 59.95 K, R1 canonical

# --- observed anchors (observations, not SM) ---
M_SUN_OBSERVED = 1.989e30
R_SUN_OBSERVED = 6.96e8
