"""
WSTP Symbolic Exporter — Unified 9-Sector UQFF Lagrangian → Wolfram .wl Export

Session 204 | PAPER_EXP_WSTP | Daniel Murphy
PURPOSE: Export the complete 9-sector UQFF Lagrangian from uqff_lagrangian_derivation.py
         into a single importable Wolfram Language (.wl) package.

ARCHITECTURE:
  uqff_lagrangian_derivation.py (9-sector Python)
        ↓ this module
  uqff_lagrangian_unified.wl (Wolfram Language package)
        ↓
  Wolfram Mathematica / WSTP kernel  (<< uqff_lagrangian_unified.wl)

SECTORS:
  1. Einstein-Hilbert (L_EH)      2. Yang-Mills (L_YM)
  3. Dirac (L_Dirac)              4. Scalar-Higgs-Vacuum (L_phi)
  5. Magnetic-Dipole (L_mag)      6. Buoyancy-Archimedes (L_buoy)
  7. Aether-Tensor (L_aether)     8. LENR-Resonance (L_LENR)
  9. Kaluza-Klein-26D (L_KK)
"""

import math
import json
import os
import sys
from typing import Dict, List, Optional
from datetime import datetime

# ── §1  UQFF CALIBRATED CONSTANTS ──────────────────────────────────────────

G       = 6.67430e-11       # m³ kg⁻¹ s⁻²
c       = 2.99792e8         # m/s
hbar    = 1.05457e-34       # J·s
PI      = math.pi
M_sun   = 1.98892e30        # kg

KAPPA   = 5.787e-9          # s⁻¹ (0.0005/day)
SSQ     = 0.57              # [SSq] dimensionless
H_SCM   = 0.99              # superconductive manifold metric
BETA_I  = 0.603             # buoyancy reversal coefficient
U_UA    = 1e-4              # aether velocity fraction
RHO_UA  = 7.09e-36          # kg/m³ aether density
RHO_SCM = 7.09e-37          # kg/m³ SCm density
ETA_AETHER = 1e-22          # aether tensor coupling
K_ETA   = 1e-113            # J·s/m³ quantum coupling
E_REACT = 1e46              # J reactor energy scale
F_TRZ   = 0.1               # time-reversal zone factor

# GW damping constants (PAPER_001–009)
D_AETHER = 1.000
D_SCM    = 1.000
D_TRZ    = 0.900
D_STRING = 0.370
F_COMBINED = D_AETHER * D_SCM * D_TRZ * D_STRING  # 0.333


# ── §2  WOLFRAM LANGUAGE SECTOR DEFINITIONS ─────────────────────────────────

def _wl_header() -> str:
    """Package header with metadata."""
    ts = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    return f"""\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* UQFF Unified 9-Sector Lagrangian — Wolfram Language Package           *)
(* Generated: {ts}                                                       *)
(* Source: wstp_symbolic_exporter.py (Star-Magic Session 204)            *)
(* Author: Daniel Murphy                                                  *)
(* Usage: << "uqff_lagrangian_unified.wl"                                *)
(* ═══════════════════════════════════════════════════════════════════════ *)

BeginPackage["UQFFLagrangian`"];

(* Public symbols *)
LUQFF::usage = "LUQFF[params] returns the unified 9-sector UQFF Lagrangian density.";
LEH::usage   = "LEH[R] returns the Einstein-Hilbert Lagrangian density.";
LYM::usage   = "LYM[FmunuSq] returns the Yang-Mills Lagrangian density.";
LDirac::usage = "LDirac[m, yij] returns the Dirac sector Lagrangian density.";
LPhi::usage  = "LPhi[phiH, phi4, v, lambda, kappa, SSq] returns the scalar-Higgs-vacuum density.";
LMag::usage  = "LMag[BField, rhoSCm, vSCm, r, Rb] returns the magnetic-dipole density.";
LBuoy::usage = "LBuoy[Ug1, Ug2, Ug3, Ug4, betai, OmegaG, M, dg, UA, tn] returns the buoyancy density.";
LAether::usage = "LAether[rhoA, vUA, eta, tn] returns the aether-tensor density.";
LLENR::usage = "LLENR[chi, omegaLENR, omegaAct, kLENR, lambdaAct, sigmaN, t] returns the LENR-resonance density.";
LKK::usage   = "LKK[R22, kappa22, a, ma] returns the Kaluza-Klein 26D density.";

(* Force functions *)
FUBii::usage = "FUBii[params] returns the total F_U_Bi_i force from E-L derivation.";
Ug1Func::usage  = "Ug1Func[G, muS, r] returns Ug1 magnetic dipole gravity.";
Ug2Func::usage  = "Ug2Func[G, muS, Ereact, r] returns Ug2 charge-reactivity gravity.";
Ug3Func::usage  = "Ug3Func[G, M, r, omegaS, t] returns Ug3 string rotation gravity.";
Ug4Func::usage  = "Ug4Func[G, M, r, SSq, rhoSCm, rhoUA] returns Ug4 vacuum concentration.";
UbiFunc::usage  = "UbiFunc[betai, Ugi, OmegaG, M, dg, UA, tn] returns buoyancy force Ubi_i.";
UmFunc::usage   = "UmFunc[G, M, r, SCm] returns Um magnetic torque.";

(* GW damping *)
DTotal::usage   = "DTotal[DAether, DSCm, DTRZ, DString] returns total GW damping factor.";
HUQFF::usage    = "HUQFF[hGR, Dtotal] returns UQFF-damped GW strain amplitude.";

(* Constants *)
UQFFConstants::usage = "UQFFConstants[] returns Association of all calibrated UQFF constants.";

Begin["`Private`"];
"""


def _wl_constants() -> str:
    """Calibrated constants block."""
    return f"""\
(* ── §2.1 Calibrated Constants ── *)
cLight     = {c};             (* m/s *)
GNewton    = {G};             (* m^3 kg^-1 s^-2 *)
hbarConst  = {hbar};          (* J s *)
piConst    = Pi;
MSun       = {M_sun};         (* kg *)
kappaUQFF  = {KAPPA};         (* s^-1,  0.0005/day *)
SSqConst   = {SSQ};           (* [SSq] string squeezing *)
HScmConst  = {H_SCM};         (* SCm metric *)
betaI      = {BETA_I};        (* buoyancy reversal *)
uUA        = {U_UA};          (* aether velocity fraction *)
rhoUA      = {RHO_UA};        (* kg/m^3 aether density *)
rhoSCm     = {RHO_SCM};       (* kg/m^3 SCm density *)
etaAether  = {ETA_AETHER};    (* aether tensor coupling *)
kEta       = {K_ETA};         (* J s/m^3 quantum coupling *)
Ereact     = {E_REACT};       (* J reactor energy scale *)
fTRZ       = {F_TRZ};         (* time-reversal zone factor *)

(* GW damping (PAPER_001-009) *)
dAether     = {D_AETHER};
dSCm        = {D_SCM};
dTRZ        = {D_TRZ};
dString     = {D_STRING};
fCombined   = {F_COMBINED};   (* 0.333 = 66.7% strain reduction *)

UQFFConstants[] := <|
  "G" -> GNewton, "c" -> cLight, "hbar" -> hbarConst,
  "kappa" -> kappaUQFF, "SSq" -> SSqConst, "H_SCm" -> HScmConst,
  "beta_i" -> betaI, "U_UA" -> uUA, "rho_UA" -> rhoUA,
  "rho_SCm" -> rhoSCm, "eta" -> etaAether, "k_eta" -> kEta,
  "D_TRZ" -> dTRZ, "D_String" -> dString, "F_combined" -> fCombined
|>;
"""


def _wl_sector_1() -> str:
    """Sector 1: Einstein-Hilbert."""
    return """\
(* ── §3.1 Sector 1: Einstein-Hilbert ── *)
(* L_EH = c^4/(16 pi G) R *)
LEH[R_] := cLight^4 / (16 piConst GNewton) * R;

(* EOM: G_munu = 8 pi G / c^4  T_munu *)
EinsteinEOM[Tmunu_] := 8 piConst GNewton / cLight^4 * Tmunu;

(* Baseline gravity *)
FGravityBaseline[M_, r_] := GNewton * M / r^2;
"""


def _wl_sector_2() -> str:
    """Sector 2: Yang-Mills."""
    return """\
(* ── §3.2 Sector 2: Yang-Mills ── *)
(* L_YM = -1/4 F^a_munu F_a^munu *)
LYM[FmunuSq_] := -1/4 * FmunuSq;

(* EOM: D_nu F^{a mu nu} = J^{a mu} *)
(* Yields Ug3 — string rotation gravity *)
Ug3Func[G_, M_, r_, omegaS_, t_] := G * M / r^2 * Sin[omegaS * t];

(* F_quark from QCD coupling *)
FQuark[alphaS_, r_] := alphaS / r^2 * (1 + alphaS / piConst);
"""


def _wl_sector_3() -> str:
    """Sector 3: Dirac."""
    return """\
(* ── §3.3 Sector 3: Dirac ── *)
(* L_Dirac = psibar (i gamma^mu D_mu - m) psi + y_ij Lbar_i Htilde N_Rj *)
LDirac[m_, yij_] := Module[{mNu},
  mNu = yij^2 * (246)^2 / (2 * 10^15);  (* seesaw: m_nu ~ y^2 v^2 / M_R *)
  <|"m_fermion" -> m, "m_neutrino_eV" -> mNu, "y_ij" -> yij|>
];

(* Neutrino force *)
FNeutrino[mNu_, r_] := GNewton * mNu / r^2;

(* Neutron interaction *)
FNeutron[sigmaN_, nDensity_] := sigmaN * nDensity * cLight;
"""


def _wl_sector_4() -> str:
    """Sector 4: Scalar-Higgs-Vacuum."""
    return """\
(* ── §3.4 Sector 4: Scalar-Higgs-Vacuum ── *)
(* L_phi = |grad phi_H|^2 - lambda(phi_H^2 - v^2/2)^2 + |d phi_4|^2 - V(phi_4) + kappa [SSq] phi_4^2 *)
LPhi[phiH_, phi4_, v_, lambda_, kappa_, SSq_] :=
  phiH^2 - lambda * (phiH^2 - v^2/2)^2 + phi4^2 - lambda/4 * phi4^4 + kappa * SSq * phi4^2;

(* EOM: Box phi_4 + V'(phi_4) - kappa [SSq] phi_4 = 0 *)
(* Yields Ug4 — vacuum concentration gravity *)
Ug4Func[G_, M_, r_, SSq_, rhoSCm_, rhoUA_] :=
  G * M / r^2 * SSq * (rhoSCm / rhoUA);

(* Dark matter force *)
FDark[M_, r_, SSq_] := GNewton * M / r^2 * SSq;
"""


def _wl_sector_5() -> str:
    """Sector 5: Magnetic-Dipole."""
    return """\
(* ── §3.5 Sector 5: Magnetic-Dipole ── *)
(* L_mag = mu0/(8 pi) |curl A_SCm|^2 - 1/2 rho_SCm |v_SCm|^2 Theta(r - R_b) *)
LMag[BField_, rhoSCm_, vSCm_, r_, Rb_] := Module[{mu0 = 4 piConst * 10^(-7)},
  mu0 / (8 piConst) * BField^2 - 1/2 * rhoSCm * vSCm^2 * If[r > Rb, 1, 0]
];

(* Yields Ug1 — magnetic dipole gravity *)
Ug1Func[G_, muS_, r_] := G * muS^2 / r^4;

(* Yields Ug2 — charge-reactivity gravity *)
Ug2Func[G_, muS_, Ereact_, r_] := G * muS * Ereact / (r^3 * cLight^2);

(* Magnetic torque Um *)
UmFunc[G_, M_, r_, SCm_] := G * M / r^2 * SCm;

(* Dark energy force *)
FDE[Lambda_, r_] := Lambda * cLight^2 / 3 * r;

(* Torque force *)
FTorque[muS_, BField_, r_] := muS * BField / r^3;
"""


def _wl_sector_6() -> str:
    """Sector 6: Buoyancy-Archimedes."""
    return """\
(* ── §3.6 Sector 6: Buoyancy-Archimedes ── *)
(* L_buoy = -beta_i Sum_i Ug_i * Omega_g * M/d_g * (1 + eps rho_sw) [UA] Cos[pi t_n] + ... *)
LBuoy[Ug1_, Ug2_, Ug3_, Ug4_, betai_, OmegaG_, M_, dg_, UA_, tn_] := Module[
  {SumUg, buoyancy},
  SumUg = Ug1 + Ug2 + Ug3 + Ug4;
  buoyancy = -betai * SumUg * OmegaG * M / dg * UA * Cos[piConst * tn];
  buoyancy
];

(* Buoyancy force Ubi_i for each channel *)
UbiFunc[betai_, Ugi_, OmegaG_, M_, dg_, UA_, tn_] :=
  -betai * Ugi * OmegaG * M / dg * UA * Cos[piConst * tn];

(* Total buoyancy: Sum Ubi_1..4 *)
UbiTotal[Ug1_, Ug2_, Ug3_, Ug4_, betai_, OmegaG_, M_, dg_, UA_, tn_] :=
  Sum[UbiFunc[betai, {Ug1, Ug2, Ug3, Ug4}[[i]], OmegaG, M, dg, UA, tn], {i, 1, 4}];
"""


def _wl_sector_7() -> str:
    """Sector 7: Aether-Tensor."""
    return """\
(* ── §3.7 Sector 7: Aether-Tensor ── *)
(* L_aether = 1/2 eta rho_A v_UA^2 Cos[pi t_n] g^munu g_munu *)
LAether[rhoA_, vUA_, eta_, tn_] :=
  1/2 * eta * rhoA * vUA^2 * Cos[piConst * tn] * 4;  (* Tr(g^mu_nu) = 4 in 4D *)

(* Aether trace correction to metric *)
AetherTrace[eta_, Ts00_, tn_] := 1 + eta * Ts00 * Cos[piConst * tn];

(* Force from aether-tensor *)
FAetherTrace[eta_, Ts00_, M_, r_, tn_] :=
  GNewton * M / r^2 * eta * Ts00 * Cos[piConst * tn];
"""


def _wl_sector_8() -> str:
    """Sector 8: LENR-Resonance."""
    return """\
(* ── §3.8 Sector 8: LENR-Resonance ── *)
(* L_LENR = 1/2 k chi'^2 - 1/2 omega^2 chi^2 + lambda chi Cos[omega_act t]
            + 1/2 sigma_n(omega) chi^2 Exp[-(omega - omega_LENR)^2 / (2 Domega^2)] *)
LLENR[chi_, omegaLENR_, omegaAct_, kLENR_, lambdaAct_, sigmaN_, t_] :=
  1/2 * kLENR * chi^2 - 1/2 * omegaLENR^2 * chi^2 +
  lambdaAct * chi * Cos[omegaAct * t] +
  1/2 * sigmaN * chi^2;

(* Resonance forces *)
FLENR[kLENR_, chi_, omegaLENR_] := kLENR * omegaLENR^2 * chi;
FAct[lambdaAct_, omegaAct_, t_] := lambdaAct * Cos[omegaAct * t];
FRes[sigmaN_, chi_, omegaLENR_, omega_, Domega_] :=
  sigmaN * chi * Exp[-(omega - omegaLENR)^2 / (2 Domega^2)];
"""


def _wl_sector_9() -> str:
    """Sector 9: Kaluza-Klein 26D."""
    return """\
(* ── §3.9 Sector 9: Kaluza-Klein 26D ── *)
(* L_KK = 1/V_22 * Integral[Sqrt[-g22] (R22/(2 kappa22^2) + |da|^2 - ma^2 a^2), {y, S22}] *)
LKK[R22_, kappa22_, a_, ma_] :=
  R22 / (2 kappa22^2) + a^2 - ma^2 * a^2;

(* ALP-photon coupling *)
FALP[gaGammaGamma_, EdotB_] := gaGammaGamma * EdotB;

(* Large extra dimension force *)
FLED[M_, r_, n_: 22] := GNewton * M / r^(2 + n) * r^n;  (* simplification *)
"""


def _wl_unified_lagrangian() -> str:
    """Unified L_UQFF = sum of all 9 sectors."""
    return """\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §4  UNIFIED LAGRANGIAN ASSEMBLY                                        *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* L_UQFF = Sqrt[-g] [L_EH + L_YM + L_Dirac + L_phi + L_mag + L_buoy + L_aether + L_LENR + L_KK] *)
LUQFF[params_Association] := Module[
  {R, FmunuSq, m, yij, phiH, phi4, v, lambda, kappa, SSq,
   BField, rhoSCmLocal, vSCm, r, Rb,
   Ug1, Ug2, Ug3, Ug4, OmegaG, M, dg, UA, tn,
   rhoA, vUA, eta,
   chi, omegaLENR, omegaAct, kLENR, lambdaAct, sigmaN, t,
   R22, kappa22, a, ma,
   sectors},

  (* Extract parameters with defaults *)
  R         = Lookup[params, "R", 0];
  FmunuSq   = Lookup[params, "FmunuSq", 0];
  m         = Lookup[params, "m", 0.511 * 10^6 * 1.602*^-19 / cLight^2];
  yij       = Lookup[params, "yij", 10^(-6)];
  phiH      = Lookup[params, "phiH", 246];
  phi4      = Lookup[params, "phi4", 10^(-10)];
  v         = Lookup[params, "v", 246];
  lambda    = Lookup[params, "lambda", 0.13];
  kappa     = Lookup[params, "kappa", kappaUQFF];
  SSq       = Lookup[params, "SSq", SSqConst];
  BField    = Lookup[params, "BField", 10^4];
  rhoSCmLocal = Lookup[params, "rhoSCm", rhoSCm];
  vSCm      = Lookup[params, "vSCm", 10^3];
  r         = Lookup[params, "r", 2.44*^20];
  Rb        = Lookup[params, "Rb", 10^4];
  M         = Lookup[params, "M", 4.3*^6 * MSun];
  dg        = Lookup[params, "dg", 2.44*^20];
  OmegaG    = Lookup[params, "OmegaG", 7.27*^-12];
  UA        = Lookup[params, "UA", uUA];
  tn        = Lookup[params, "tn", 0.5];
  rhoA      = Lookup[params, "rhoA", rhoUA];
  vUA       = Lookup[params, "vUA", uUA * cLight];
  eta       = Lookup[params, "eta", etaAether];
  chi       = Lookup[params, "chi", 10^(-15)];
  omegaLENR = Lookup[params, "omegaLENR", 1.25*^12];
  omegaAct  = Lookup[params, "omegaAct", 300 * 2 piConst];
  kLENR     = Lookup[params, "kLENR", 1];
  lambdaAct = Lookup[params, "lambdaAct", 10^(-20)];
  sigmaN    = Lookup[params, "sigmaN", 10^(-28)];
  t         = Lookup[params, "t", 0];
  R22       = Lookup[params, "R22", 0];
  kappa22   = Lookup[params, "kappa22", 10^10];
  a         = Lookup[params, "a", 10^(-30)];
  ma        = Lookup[params, "ma", 10^(-22)];

  (* Compute individual Ug channels for buoyancy feedback *)
  Ug1 = Ug1Func[GNewton, Lookup[params, "muS", 10^25], r];
  Ug2 = Ug2Func[GNewton, Lookup[params, "muS", 10^25], Ereact, r];
  Ug3 = Ug3Func[GNewton, M, r, Lookup[params, "omegaS", 10^(-7)], t];
  Ug4 = Ug4Func[GNewton, M, r, SSq, rhoSCmLocal, rhoA];

  (* Assemble sectors *)
  sectors = <|
    "L_EH"     -> LEH[R],
    "L_YM"     -> LYM[FmunuSq],
    "L_Dirac"  -> 0,  (* symbolic — needs spinor fields *)
    "L_phi"    -> LPhi[phiH, phi4, v, lambda, kappa, SSq],
    "L_mag"    -> LMag[BField, rhoSCmLocal, vSCm, r, Rb],
    "L_buoy"   -> LBuoy[Ug1, Ug2, Ug3, Ug4, betaI, OmegaG, M, dg, UA, tn],
    "L_aether" -> LAether[rhoA, vUA, eta, tn],
    "L_LENR"   -> LLENR[chi, omegaLENR, omegaAct, kLENR, lambdaAct, sigmaN, t],
    "L_KK"     -> LKK[R22, kappa22, a, ma]
  |>;

  <|"sectors" -> sectors,
    "L_UQFF_total" -> Total[Values[sectors]],
    "Ug_channels" -> <|"Ug1" -> Ug1, "Ug2" -> Ug2, "Ug3" -> Ug3, "Ug4" -> Ug4|>,
    "forces" -> <|
      "Ubi_total" -> UbiTotal[Ug1, Ug2, Ug3, Ug4, betaI, OmegaG, M, dg, UA, tn],
      "Um" -> UmFunc[GNewton, M, r, HScmConst],
      "A_trace" -> AetherTrace[eta, 1.27*^3, tn]
    |>
  |>
];
"""


def _wl_force_assembly() -> str:
    """F_U_Bi_i total force from E-L derivation."""
    return """\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §5  EULER-LAGRANGE FORCE ASSEMBLY                                      *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* F_U_Bi_i = Sum Ug_i + Sum Ubi_i + Um + Tr(A_munu) + F_LENR + F_act
              + F_res + F_quark + F_neutrino + F_ALP + F_dark + F_LED
              + F_neutron + F_torque + F_DE *)

FUBii[params_Association] := Module[
  {r, M, muS, SSq, OmegaG, dg, UA, tn, SCm, eta, Ts00,
   chi, kLENR, omegaLENR, omegaAct, lambdaAct, sigmaN, omega, Domega, t,
   alphaS, mNu, nDensity, gaGG, EdotB, BField, Lambda,
   Ug1, Ug2, Ug3, Ug4, Ubi1, Ubi2, Ubi3, Ubi4, Um, Atrace,
   fLENR, fAct, fRes, fQuark, fNeutrino, fALP, fDark, fLED, fNeutron, fTorque, fDE,
   total},

  (* Parameters *)
  r     = Lookup[params, "r", 2.44*^20];
  M     = Lookup[params, "M", 4.3*^6 * MSun];
  muS   = Lookup[params, "muS", 10^25];
  SSq   = Lookup[params, "SSq", SSqConst];
  OmegaG = Lookup[params, "OmegaG", 7.27*^-12];
  dg    = Lookup[params, "dg", 2.44*^20];
  UA    = Lookup[params, "UA", uUA];
  tn    = Lookup[params, "tn", 0.5];
  SCm   = Lookup[params, "SCm", HScmConst];
  eta   = Lookup[params, "eta", etaAether];
  Ts00  = Lookup[params, "Ts00", 1.27*^3];
  chi   = Lookup[params, "chi", 10^(-15)];
  kLENR = Lookup[params, "kLENR", 1];
  omegaLENR = Lookup[params, "omegaLENR", 1.25*^12];
  omegaAct = Lookup[params, "omegaAct", 300 * 2 piConst];
  lambdaAct = Lookup[params, "lambdaAct", 10^(-20)];
  sigmaN = Lookup[params, "sigmaN", 10^(-28)];
  omega  = Lookup[params, "omega", 1.25*^12];
  Domega = Lookup[params, "Domega", 10^9];
  t      = Lookup[params, "t", 0];
  alphaS = Lookup[params, "alphaS", 0.118];
  mNu    = Lookup[params, "mNu", 0.1 * 1.602*^-19 / cLight^2];
  nDensity = Lookup[params, "nDensity", 10^18];
  gaGG   = Lookup[params, "gaGammaGamma", 10^(-12)];
  EdotB  = Lookup[params, "EdotB", 0];
  BField = Lookup[params, "BField", 10^4];
  Lambda = Lookup[params, "Lambda", 1.1056*^-52];

  (* 4 gravity channels *)
  Ug1 = Ug1Func[GNewton, muS, r];
  Ug2 = Ug2Func[GNewton, muS, Ereact, r];
  Ug3 = Ug3Func[GNewton, M, r, Lookup[params, "omegaS", 10^(-7)], t];
  Ug4 = Ug4Func[GNewton, M, r, SSq, rhoSCm, rhoUA];

  (* 4 buoyancy channels *)
  Ubi1 = UbiFunc[betaI, Ug1, OmegaG, M, dg, UA, tn];
  Ubi2 = UbiFunc[betaI, Ug2, OmegaG, M, dg, UA, tn];
  Ubi3 = UbiFunc[betaI, Ug3, OmegaG, M, dg, UA, tn];
  Ubi4 = UbiFunc[betaI, Ug4, OmegaG, M, dg, UA, tn];

  (* Magnetic torque *)
  Um = UmFunc[GNewton, M, r, SCm];

  (* Aether trace *)
  Atrace = FAetherTrace[eta, Ts00, M, r, tn];

  (* LENR forces *)
  fLENR = FLENR[kLENR, chi, omegaLENR];
  fAct  = FAct[lambdaAct, omegaAct, t];
  fRes  = FRes[sigmaN, chi, omegaLENR, omega, Domega];

  (* SM forces *)
  fQuark    = FQuark[alphaS, r];
  fNeutrino = FNeutrino[mNu, r];
  fALP      = FALP[gaGG, EdotB];
  fDark     = FDark[M, r, SSq];
  fLED      = FLED[M, r, 22];
  fNeutron  = FNeutron[sigmaN, nDensity];
  fTorque   = FTorque[muS, BField, r];
  fDE       = FDE[Lambda, r];

  (* Total *)
  total = (Ug1 + Ug2 + Ug3 + Ug4) + (Ubi1 + Ubi2 + Ubi3 + Ubi4) + Um + Atrace +
          fLENR + fAct + fRes + fQuark + fNeutrino + fALP + fDark + fLED +
          fNeutron + fTorque + fDE;

  <|"F_U_Bi_i" -> total,
    "Ug_sum" -> Ug1 + Ug2 + Ug3 + Ug4,
    "Ubi_sum" -> Ubi1 + Ubi2 + Ubi3 + Ubi4,
    "Um" -> Um,
    "A_trace" -> Atrace,
    "F_external" -> fLENR + fAct + fRes + fQuark + fNeutrino + fALP + fDark + fLED + fNeutron + fTorque + fDE,
    "channels" -> <|
      "Ug1" -> Ug1, "Ug2" -> Ug2, "Ug3" -> Ug3, "Ug4" -> Ug4,
      "Ubi1" -> Ubi1, "Ubi2" -> Ubi2, "Ubi3" -> Ubi3, "Ubi4" -> Ubi4
    |>
  |>
];
"""


def _wl_gw_damping() -> str:
    """GW damping functions from PAPER_001-009."""
    return """\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §6  GRAVITATIONAL WAVE DAMPING (PAPER_001-009)                         *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* D_total = D_Aether * D_SCm * D_TRZ * D_String *)
DTotal[dA_: 1.0, dS_: 1.0, dT_: 0.900, dStr_: 0.370] := dA * dS * dT * dStr;

(* h_UQFF = h_GR * D_total *)
HUQFF[hGR_, Dtotal_: 0.333] := hGR * Dtotal;

(* Chirp mass *)
ChirpMass[m1_, m2_] := (m1 * m2)^(3/5) / (m1 + m2)^(1/5);

(* GR strain at distance *)
HGR[Mc_, f_, DL_] := Module[{},
  4 / DL * (GNewton * Mc / cLight^2)^(5/3) * (piConst * f / cLight)^(2/3)
];

(* Apparent distance bias from UQFF damping *)
ApparentDistanceBias[DL_, Dtotal_] := DL / Dtotal;  (* > DL always *)

(* Phase lag *)
PhaseLag[Dtotal_, fGW_, tMerger_] := (1 - Dtotal) * 2 piConst * fGW * tMerger;

(* Pre-computed GW events *)
GW170817Params = <|"m1" -> 1.46 MSun, "m2" -> 1.27 MSun, "DL_Mpc" -> 40,
  "D_total" -> 0.333, "h_GR" -> 5.4176*^-22, "h_UQFF" -> 1.804*^-22|>;
GW150914Params = <|"m1" -> 36 MSun, "m2" -> 29 MSun, "DL_Mpc" -> 410,
  "D_total" -> 0.810, "h_GR" -> 1.2499*^-21, "h_UQFF" -> 4.1622*^-22|>;
GW190425Params = <|"m1" -> 1.7 MSun, "m2" -> 1.5 MSun, "DL_Mpc" -> 159,
  "D_total" -> 0.530|>;
"""


def _wl_footer() -> str:
    """Package footer."""
    return """\

End[];  (* `Private` *)
EndPackage[];  (* UQFFLagrangian` *)

(* ═══════════════════════════════════════════════════════════════════════ *)
(* USAGE EXAMPLE:                                                         *)
(*   << "uqff_lagrangian_unified.wl"                                      *)
(*   result = LUQFF[<|"M" -> 4.3*^6 * 1.989*^30, "r" -> 2.44*^20|>]     *)
(*   Print["L_UQFF total = ", result["L_UQFF_total"]]                    *)
(*   forces = FUBii[<|"M" -> 4.3*^6 * 1.989*^30, "r" -> 2.44*^20|>]     *)
(*   Print["F_U_Bi_i = ", forces["F_U_Bi_i"]]                            *)
(* ═══════════════════════════════════════════════════════════════════════ *)
"""


# ── §3  EXPORT ENGINE ──────────────────────────────────────────────────────

class WSTPSymbolicExporter:
    """Export unified 9-sector UQFF Lagrangian to Wolfram Language .wl package."""

    def __init__(self, output_path: str = "uqff_lagrangian_unified.wl"):
        self.output_path = output_path

    def generate_wl(self) -> str:
        """Assemble full Wolfram Language package from all sectors."""
        sections = [
            _wl_header(),
            _wl_constants(),
            _wl_sector_1(),
            _wl_sector_2(),
            _wl_sector_3(),
            _wl_sector_4(),
            _wl_sector_5(),
            _wl_sector_6(),
            _wl_sector_7(),
            _wl_sector_8(),
            _wl_sector_9(),
            _wl_unified_lagrangian(),
            _wl_force_assembly(),
            _wl_gw_damping(),
            _wl_footer(),
        ]
        return "\n".join(sections)

    def export(self) -> str:
        """Write .wl package to disk. Returns output path."""
        wl_content = self.generate_wl()
        with open(self.output_path, "w", encoding="utf-8") as f:
            f.write(wl_content)
        return self.output_path

    def export_json_manifest(self, json_path: str = "wstp_export_manifest.json") -> str:
        """Export manifest documenting all exported symbols."""
        manifest = {
            "package": "UQFFLagrangian`",
            "generated": datetime.utcnow().isoformat() + "Z",
            "source": "wstp_symbolic_exporter.py",
            "session": 204,
            "output_file": self.output_path,
            "sectors": [
                {"id": 1, "name": "Einstein-Hilbert", "symbol": "LEH", "fields": ["g_munu"]},
                {"id": 2, "name": "Yang-Mills", "symbol": "LYM", "fields": ["A_mu_a", "B_j"]},
                {"id": 3, "name": "Dirac", "symbol": "LDirac", "fields": ["psi", "psi_bar", "N_R"]},
                {"id": 4, "name": "Scalar-Higgs-Vacuum", "symbol": "LPhi", "fields": ["phi_H", "phi_4"]},
                {"id": 5, "name": "Magnetic-Dipole", "symbol": "LMag", "fields": ["A_SCm", "mu_s", "R_b"]},
                {"id": 6, "name": "Buoyancy-Archimedes", "symbol": "LBuoy", "fields": ["Omega_g", "beta_i", "mu_j"]},
                {"id": 7, "name": "Aether-Tensor", "symbol": "LAether", "fields": ["rho_A", "v_UA", "eta"]},
                {"id": 8, "name": "LENR-Resonance", "symbol": "LLENR", "fields": ["chi", "omega_LENR", "omega_act"]},
                {"id": 9, "name": "Kaluza-Klein-26D", "symbol": "LKK", "fields": ["g_mn_22D", "a_ALP"]},
            ],
            "force_functions": [
                "FUBii", "Ug1Func", "Ug2Func", "Ug3Func", "Ug4Func",
                "UbiFunc", "UmFunc", "FGravityBaseline", "FQuark",
                "FNeutrino", "FNeutron", "FDark", "FDE", "FTorque",
                "FLENR", "FAct", "FRes", "FALP", "FLED",
                "FAetherTrace", "AetherTrace",
            ],
            "gw_functions": ["DTotal", "HUQFF", "ChirpMass", "HGR", "ApparentDistanceBias", "PhaseLag"],
            "gw_events": ["GW170817Params", "GW150914Params", "GW190425Params"],
            "constants_function": "UQFFConstants",
            "calibrated_constants": {
                "kappa": KAPPA, "SSq": SSQ, "H_SCm": H_SCM,
                "beta_i": BETA_I, "U_UA": U_UA, "rho_UA": RHO_UA,
                "rho_SCm": RHO_SCM, "eta_aether": ETA_AETHER,
                "D_TRZ": D_TRZ, "D_String": D_STRING, "F_combined": F_COMBINED,
            },
        }
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump(manifest, f, indent=2)
        return json_path


# ── §4  MAIN ────────────────────────────────────────────────────────────────

def main():
    """Generate unified UQFF Wolfram Language package."""
    print("=" * 72)
    print("WSTP Symbolic Exporter — 9-Sector UQFF Lagrangian → Wolfram .wl")
    print("=" * 72)

    exporter = WSTPSymbolicExporter()

    # Generate and write .wl package
    wl_path = exporter.export()
    wl_size = os.path.getsize(wl_path)
    print(f"\n[OK] Wolfram package written: {wl_path} ({wl_size:,} bytes)")

    # Count exported symbols
    wl_content = exporter.generate_wl()
    func_count = wl_content.count("::usage")
    print(f"[OK] Public symbols exported: {func_count}")
    print(f"[OK] Sectors: 9 (EH, YM, Dirac, Scalar, Mag, Buoy, Aether, LENR, KK)")
    print(f"[OK] GW events: 3 (GW170817, GW150914, GW190425)")

    # Export manifest
    manifest_path = exporter.export_json_manifest()
    print(f"[OK] JSON manifest: {manifest_path}")

    # Verification
    print(f"\n── Verification ──")
    print(f"  D_total = D_Aether × D_SCm × D_TRZ × D_String")
    print(f"         = {D_AETHER} × {D_SCM} × {D_TRZ} × {D_STRING}")
    print(f"         = {F_COMBINED:.3f}  (66.7% strain reduction)")
    print(f"  κ = {KAPPA:.3e} s⁻¹  (0.0005/day)")
    print(f"  [SSq] = {SSQ}")
    print(f"  β_i = {BETA_I}")

    print(f"\n── Usage ──")
    print(f'  Mathematica:  << "{wl_path}"')
    print(f'  WSTP C++:     WSTPutFunction(link, "Get", 1);')
    print(f'                WSTPutString(link, "{wl_path}");')

    print(f"\n{'=' * 72}")
    print(f"EXPORT COMPLETE")
    print(f"{'=' * 72}")


if __name__ == "__main__":
    main()
