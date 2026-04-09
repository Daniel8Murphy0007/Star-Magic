"""
Kozima WSTP Kernel — Parameterized Mathematica Export with Lattice Coupling

Session 204 | Daniel Murphy
PURPOSE: The existing wstp_symbolic_exporter.py (Sector 8, L_LENR) exports a
         GENERIC neutron interaction:
           FNeutron[sigmaN_, nDensity_] := sigmaN * nDensity * cLight;

         This has NO Kozima-specific structure: no SCm lattice coupling,
         no VDS 26-level phonon resonance, no density scaling with [SSq].

         This module generates a Kozima-PARAMETERIZED Mathematica .wl export:
           1. sigmaNSCm[ω, n] with Gaussian + VDS resonance structure
           2. FKozima[...] with full lattice coupling k_Kozima
           3. VDS layer integration across n=0..26
           4. Density-scaled form σ_n^SCm(ρ)
           5. Buoyancy-coupled neutron force F_neutron^SCm
           6. Complete UQFFKozima constants association

ARCHITECTURE: Pure exporter. Generates .wl file content from parameterized
              Kozima physics (kozima_scm_cross_section.py constants).
"""

import math
import json
import os
from datetime import datetime, timezone
from typing import Dict, Optional

# ── §0  CONSTANTS (from kozima_scm_cross_section.py) ─────────────────────

PI = math.pi
C = 2.998e8
HBAR = 1.055e-34
G = 6.674e-11
M_SUN = 1.989e30

# UQFF calibrated
SSQ = 0.57
KAPPA = 5.787e-9
BETA_I = 0.603
H_SCM = 0.99
U_UA = 1e-4

# Vacuum densities
RHO_SCM = 7.09e-37
RHO_UA = 7.09e-36
RHO_VAC_SCM = 9.47e-27
RHO_VAC_UA = 5e-27
ETA_AETHER = 1e-22
K_ETA = 1e-113
E_REACT = 1e46
F_TRZ = 0.1

# LENR / Kozima parameters
F_LENR_THZ = 1.25e12
OMEGA_SCM = 2 * PI * F_LENR_THZ
GAMMA_DEFAULT = 0.1e12 * 2 * PI
SIGMA_0 = 1e-4
N_LEVELS = 26
RHO_CRIT_DEFAULT = 1e-20

# GW damping
D_AETHER = 1.000
D_SCM = 1.000
D_TRZ = 0.900
D_STRING = 0.370
F_COMBINED = D_AETHER * D_SCM * D_TRZ * D_STRING


# ── §1  WOLFRAM LANGUAGE GENERATOR ───────────────────────────────────────

class KozimaWSTPKernel:
    """
    Generates Kozima-parameterized Wolfram Language (.wl) package content.

    Replaces the generic FNeutron[sigmaN_, nDensity_] with full
    Kozima SCm-modulated neutron physics:
      - Frequency-dependent cross-section with VDS structure
      - Density-scaled cross-section with [SSq] exponent
      - Buoyancy-coupled neutron production force
      - Lattice coupling constants

    Output is a complete .wl package importable via:
      << "uqff_kozima_kernel.wl"
    """

    def __init__(self):
        self.timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    # ── 1a. Package header ───────────────────────────────────────────────

    def _wl_header(self) -> str:
        return f"""\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* UQFF Kozima Neutron-Drop Kernel — Wolfram Language Package            *)
(* Generated: {self.timestamp}                                           *)
(* Source: kozima_wstp_kernel.py (Star-Magic Session 204)                *)
(* Author: Daniel Murphy                                                  *)
(* Usage: << "uqff_kozima_kernel.wl"                                     *)
(* Extends: wstp_symbolic_exporter.py Sector 8 (L_LENR)                  *)
(* ═══════════════════════════════════════════════════════════════════════ *)

BeginPackage["UQFFKozima`"];

(* ── Public Symbols ── *)
SigmaNSCm::usage       = "SigmaNSCm[omega, n] returns SCm-modulated Kozima cross-section with VDS.";
SigmaNBase::usage       = "SigmaNBase[omega] returns base CP2 cross-section for comparison.";
SigmaNDensity::usage    = "SigmaNDensity[rho] returns density-scaled SCm cross-section.";
FKozima::usage          = "FKozima[omega, n, Nn, PhiPhonon, FUBi, FU] returns Kozima neutron production force.";
FKozimaDamping::usage   = "FKozimaDamping[FNeutron, FGW] returns D_Kozima damping factor for GW inspiral.";
VDSIntegral::usage      = "VDSIntegral[omega] returns integral of SigmaNSCm over all 26 VDS levels.";
ReversalCondition::usage = "ReversalCondition[FU, Nn, sigma, Phi] returns F_{{U,Bi}} threshold for buoyancy reversal.";
KozimaConstants::usage  = "KozimaConstants[] returns Association of all Kozima-specific constants.";
KozimaValidation::usage = "KozimaValidation[] runs self-test and returns validation report.";

Begin["`Private`"];
"""

    # ── 1b. Constants block ──────────────────────────────────────────────

    def _wl_constants(self) -> str:
        return f"""\
(* ── §2 Kozima-Specific Physical Constants ── *)
(* Inherited from UQFFLagrangian` package: cLight, GNewton, hbarConst, etc. *)
(* These are Kozima-specific additions:                                     *)

cLightKz       = {C};                    (* m/s *)
GNewtonKz      = {G};                    (* m^3 kg^-1 s^-2 *)
hbarKz         = {HBAR};                 (* J s *)
piKz           = Pi;
MSunKz         = {M_SUN};               (* kg *)

(* UQFF calibrated *)
SSqKz          = {SSQ};                   (* string squeezing *)
kappaKz        = {KAPPA};                (* s^-1 *)
betaIKz        = {BETA_I};               (* buoyancy reversal coefficient *)
HScmKz         = {H_SCM};                (* SCm metric *)
uUAKz          = {U_UA};                 (* aether velocity fraction *)

(* Vacuum densities *)
rhoSCmKz       = {RHO_SCM};              (* kg/m^3 SCm vacuum *)
rhoUAKz        = {RHO_UA};               (* kg/m^3 UA vacuum *)
rhoVacSCmKz    = {RHO_VAC_SCM};          (* kg/m^3 VDS SCm *)
rhoVacUAKz     = {RHO_VAC_UA};           (* kg/m^3 VDS UA *)
rhoCritKz      = {RHO_CRIT_DEFAULT};     (* kg/m^3 SCm activation threshold *)

(* LENR / Kozima parameters *)
fLENRTHz       = {F_LENR_THZ};           (* Hz — SCm phonon resonance *)
omegaSCmKz     = 2 piKz fLENRTHz;        (* rad/s angular frequency *)
GammaKz        = 2 piKz {0.1e12};        (* rad/s resonance width ~0.1 THz *)
sigma0Kz       = {SIGMA_0};              (* base cross-section *)
nLevelsKz      = {N_LEVELS};             (* VDS dimensionality *)

(* GW damping factors (PAPER_001-009) *)
dAetherKz      = {D_AETHER};
dSCmKz         = {D_SCM};
dTRZKz         = {D_TRZ};
dStringKz      = {D_STRING};
fCombinedKz    = {F_COMBINED};            (* 4-channel product *)

KozimaConstants[] := <|
  "sigma_0" -> sigma0Kz,
  "omega_SCm" -> omegaSCmKz,
  "f_LENR_THz" -> fLENRTHz / 10^12,
  "Gamma" -> GammaKz,
  "SSq" -> SSqKz,
  "beta_i" -> betaIKz,
  "n_levels" -> nLevelsKz,
  "rho_SCm" -> rhoSCmKz,
  "rho_UA" -> rhoUAKz,
  "rho_vac_SCm" -> rhoVacSCmKz,
  "rho_vac_UA" -> rhoVacUAKz,
  "rho_crit" -> rhoCritKz,
  "D_Aether" -> dAetherKz,
  "D_SCm" -> dSCmKz,
  "D_TRZ" -> dTRZKz,
  "D_String" -> dStringKz,
  "F_combined_4ch" -> fCombinedKz
|>;
"""

    # ── 1c. Cross-section functions ──────────────────────────────────────

    def _wl_cross_sections(self) -> str:
        return """\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §3  KOZIMA SCm-MODULATED CROSS-SECTIONS                               *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* §3.1 Base CP2 cross-section (for comparison) *)
(* sigma_n(omega) = sigma_0 * (omega/omega_LENR)^2 * Exp[-(omega - omega_LENR)^2/(2 Domega^2)] *)
SigmaNBase[omega_] := Module[{freqRatio, gaussian, Domega},
  Domega = GammaKz;
  freqRatio = (omega / omegaSCmKz)^2;
  gaussian = Exp[-(omega - omegaSCmKz)^2 / (2 Domega^2)];
  sigma0Kz * freqRatio * gaussian
];

(* §3.2 SCm-modulated cross-section with VDS 26-level structure *)
(* sigma_n^SCm(omega, n) = sigma_0 * Exp[-(omega - omega_SCm)^2/(2 Gamma^2)] * (1 + [SSq]*n/26) *)
SigmaNSCm[omega_, n_] := Module[{gaussian, vdsFactor},
  gaussian = Exp[-(omega - omegaSCmKz)^2 / (2 GammaKz^2)];
  vdsFactor = 1 + SSqKz * n / nLevelsKz;
  sigma0Kz * gaussian * vdsFactor
];

(* §3.3 Density-scaled cross-section with VDS power-law and cutoff *)
(* sigma_n^SCm(rho) = sigma_0 * (rho_SCm/rho_UA)^[SSq] * Exp(-[SSq]*rho/rho_crit) *)
SigmaNDensity[rho_] := Module[{vdsPower, expCutoff},
  vdsPower = (rhoVacSCmKz / rhoVacUAKz)^SSqKz;
  expCutoff = Exp[-SSqKz * rho / rhoCritKz];
  sigma0Kz * vdsPower * expCutoff
];

(* §3.4 VDS integral over all 26 levels *)
(* Sum_{n=0}^{26} sigma_n^SCm(omega, n) *)
VDSIntegral[omega_] := Sum[SigmaNSCm[omega, n], {n, 0, nLevelsKz}];

(* §3.5 Symbolic VDS integral (closed form) *)
(* Sum_{n=0}^{26} (1 + SSq*n/26) = 27 + SSq * 26*27/(2*26) = 27(1 + SSq/2) *)
VDSIntegralClosed[omega_] := Module[{gaussian, closedSum},
  gaussian = Exp[-(omega - omegaSCmKz)^2 / (2 GammaKz^2)];
  closedSum = (nLevelsKz + 1) * (1 + SSqKz / 2);
  sigma0Kz * gaussian * closedSum
];
"""

    # ── 1d. Force functions ──────────────────────────────────────────────

    def _wl_forces(self) -> str:
        return """\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §4  KOZIMA NEUTRON-DROP FORCE FUNCTIONS                                *)
(* ═══════════════════════════════════════════════════════════════════════ *)

(* §4.1 Full Kozima neutron production force *)
(* F_neutron^SCm = N_n * sigma_n^SCm(omega, n) * Phi_phonon * (F_{U,Bi}/F_U - 1) *)
FKozima[omega_, n_, Nn_, PhiPhonon_, FUBi_, FU_] := Module[
  {sigmaSCm, buoyancyRev, Fneutron},
  sigmaSCm = SigmaNSCm[omega, n];
  buoyancyRev = If[FU == 0, 0, FUBi / FU - 1];
  Fneutron = Nn * sigmaSCm * PhiPhonon * buoyancyRev;
  <|
    "F_neutron" -> Fneutron,
    "sigma_n_SCm" -> sigmaSCm,
    "buoyancy_reversal" -> buoyancyRev,
    "N_n" -> Nn,
    "Phi_phonon" -> PhiPhonon,
    "VDS_level" -> n
  |>
];

(* §4.2 Kozima GW damping factor (5th channel) *)
(* D_Kozima = 1 / (1 + |F_neutron| / F_GW) *)
FKozimaDamping[FNeutron_, FGW_] := Module[{ratio, DKozima, Dtotal5},
  ratio = If[FGW == 0, 0, Abs[FNeutron] / FGW];
  DKozima = 1 / (1 + ratio);
  Dtotal5 = fCombinedKz * DKozima;
  <|
    "D_Kozima" -> DKozima,
    "D_total_4ch" -> fCombinedKz,
    "D_total_5ch" -> Dtotal5,
    "ratio" -> ratio,
    "additional_reduction_pct" -> (fCombinedKz - Dtotal5) / fCombinedKz * 100
  |>
];

(* §4.3 Variational reversal condition *)
(* F_{U,Bi} = F_U * (1 + N_n * sigma_n * Phi / F_U)  *)
ReversalCondition[FU_, Nn_, sigma_, Phi_] := Module[{LambdaK, FUBi, ratio},
  LambdaK = Nn * sigma * Phi;
  FUBi = If[FU == 0, LambdaK, FU * (1 + LambdaK / FU)];
  ratio = If[FU == 0, Infinity, FUBi / FU];
  <|
    "F_U_Bi_threshold" -> FUBi,
    "F_U" -> FU,
    "Lambda_K" -> LambdaK,
    "reversal_ratio" -> ratio,
    "condition" -> "F_{U,Bi} = F_U (1 + N_n sigma_n Phi / F_U)"
  |>
];

(* §4.4 Effective mass from buoyancy sector *)
(* m_eff^2 = beta_i * Ug_sum * Omega_g * M / (d_g * c^2 * hbar^2) * [UA] *)
MEffBuoyancy[UgList_, OmegaG_, M_, dg_] := Module[{UgSum, mEffSq, mEff, lambdaBuoy},
  UgSum = Total[Abs /@ UgList];
  mEffSq = betaIKz * UgSum * OmegaG * M / (dg * cLightKz^2 * hbarKz^2) * uUAKz;
  mEff = Sqrt[Abs[mEffSq]];
  lambdaBuoy = If[mEff > 0, hbarKz / (mEff * cLightKz), Infinity];
  <|
    "m_eff_sq" -> mEffSq,
    "m_eff" -> mEff,
    "lambda_buoy" -> lambdaBuoy,
    "Ug_sum" -> UgSum
  |>
];
"""

    # ── 1e. Validation ───────────────────────────────────────────────────

    def _wl_validation(self) -> str:
        return f"""\
(* ═══════════════════════════════════════════════════════════════════════ *)
(* §5  SELF-VALIDATION                                                    *)
(* ═══════════════════════════════════════════════════════════════════════ *)

KozimaValidation[] := Module[
  {{sigmaBase, sigmaSCm, sigmaDens, fk, dk, rc, vi, viClosed, allPass}},

  (* Test 1: Base cross-section at resonance *)
  sigmaBase = SigmaNBase[omegaSCmKz];

  (* Test 2: SCm cross-section at resonance, n=13 *)
  sigmaSCm = SigmaNSCm[omegaSCmKz, 13];

  (* Test 3: Density cross-section *)
  sigmaDens = SigmaNDensity[rhoVacSCmKz];

  (* Test 4: Kozima force *)
  fk = FKozima[omegaSCmKz, 13, 10^28, 10^16, 1.1, 1.0];

  (* Test 5: Damping factor *)
  dk = FKozimaDamping[fk["F_neutron"], 10^40];

  (* Test 6: Reversal condition *)
  rc = ReversalCondition[1.0, 10^28, sigmaSCm, 10^16];

  (* Test 7: VDS integral consistency *)
  vi = VDSIntegral[omegaSCmKz];
  viClosed = VDSIntegralClosed[omegaSCmKz];

  allPass = And[
    sigmaBase > 0,
    sigmaSCm > sigmaBase,          (* SCm should exceed base at resonance *)
    sigmaDens > 0,
    fk["F_neutron"] > 0,           (* Positive when F_UBi > F_U *)
    0 < dk["D_Kozima"] <= 1,       (* Damping in [0,1] *)
    rc["F_U_Bi_threshold"] >= 1.0, (* Threshold >= F_U *)
    Abs[vi - viClosed] / vi < 0.01 (* Closed form matches sum to 1% *)
  ];

  <|
    "status" -> If[allPass, "VALIDATED", "FAILED"],
    "timestamp" -> "{self.timestamp}",
    "tests" -> <|
      "sigma_base_resonance" -> sigmaBase,
      "sigma_SCm_resonance_n13" -> sigmaSCm,
      "sigma_density" -> sigmaDens,
      "F_neutron" -> fk["F_neutron"],
      "D_Kozima" -> dk["D_Kozima"],
      "D_total_5ch" -> dk["D_total_5ch"],
      "F_UBi_threshold" -> rc["F_U_Bi_threshold"],
      "VDS_integral" -> vi,
      "VDS_closed_form" -> viClosed,
      "VDS_match_pct" -> 100 (1 - Abs[vi - viClosed] / vi)
    |>,
    "all_passed" -> allPass
  |>
];
"""

    # ── 1f. Footer ───────────────────────────────────────────────────────

    def _wl_footer(self) -> str:
        return """\

End[];  (* `Private` *)
EndPackage[];

(* ═══════════════════════════════════════════════════════════════════════ *)
(* END — UQFFKozima` package                                              *)
(* To use: << "uqff_kozima_kernel.wl"                                     *)
(* Then:   KozimaValidation[]  to self-test                               *)
(* ═══════════════════════════════════════════════════════════════════════ *)
"""

    # ── 1g. Full package assembly ────────────────────────────────────────

    def generate_wl_package(self) -> str:
        """
        Assemble the complete Kozima .wl package from all sections.
        Returns the full Wolfram Language source as a string.
        """
        sections = [
            self._wl_header(),
            self._wl_constants(),
            self._wl_cross_sections(),
            self._wl_forces(),
            self._wl_validation(),
            self._wl_footer(),
        ]
        return "\n".join(sections)

    # ── 1h. Export to file ───────────────────────────────────────────────

    def export_wl_file(self, filepath: str = "uqff_kozima_kernel.wl") -> Dict:
        """
        Generate and write the .wl package to disk.
        """
        content = self.generate_wl_package()
        line_count = content.count("\n") + 1

        with open(filepath, "w", encoding="utf-8") as f:
            f.write(content)

        return {
            "filepath": os.path.abspath(filepath),
            "lines": line_count,
            "size_bytes": len(content.encode("utf-8")),
            "sections": [
                "Header + Public symbols",
                "Kozima constants",
                "Cross-section functions (5: Base, SCm, Density, VDS integral, VDS closed)",
                "Force functions (4: FKozima, FKozimaDamping, ReversalCondition, MEffBuoyancy)",
                "Self-validation",
                "Footer",
            ],
            "timestamp": self.timestamp,
        }

    # ── 1i. Summary report ───────────────────────────────────────────────

    def generate_summary(self) -> Dict:
        """
        Return a structured summary of what this kernel provides
        vs. the generic FNeutron in wstp_symbolic_exporter.py.
        """
        return {
            "replaces": "FNeutron[sigmaN_, nDensity_] := sigmaN * nDensity * cLight  (generic, Sector 3)",
            "provides": [
                "SigmaNSCm[omega, n] — SCm Gaussian + VDS 26-level cross-section",
                "SigmaNBase[omega] — CP2 comparison form",
                "SigmaNDensity[rho] — density-scaled with [SSq] exponent",
                "FKozima[omega, n, Nn, Phi, FUBi, FU] — full neutron production force",
                "FKozimaDamping[FN, FGW] — 5th GW damping channel",
                "ReversalCondition[FU, Nn, sigma, Phi] — buoyancy reversal threshold",
                "VDSIntegral[omega] — sum over 26 VDS levels",
                "VDSIntegralClosed[omega] — analytical closed form",
                "MEffBuoyancy[UgList, Omega, M, dg] — effective boson mass",
                "KozimaConstants[] — all Kozima parameters as Association",
                "KozimaValidation[] — built-in self-test",
            ],
            "architecture": "Standalone BeginPackage[UQFFKozima`] — imports alongside UQFFLagrangian`",
            "lattice_coupling": "σ_n^SCm includes phonon resonance at 1.25 THz, VDS (1+[SSq]n/26) structure, "
                                "density scaling (ρ_SCm/ρ_UA)^[SSq], and buoyancy reversal factor",
        }


# ── §2  SELF-TEST ─────────────────────────────────────────────────────────

def main():
    """Validate Kozima WSTP kernel generator."""
    print("=" * 72)
    print("KOZIMA WSTP KERNEL — MATHEMATICA EXPORT VALIDATION")
    print("=" * 72)

    kernel = KozimaWSTPKernel()

    # Generate the .wl content
    wl_content = kernel.generate_wl_package()
    lines = wl_content.count("\n") + 1

    print(f"Generated Wolfram Language package: {lines} lines")
    print()

    # Verify key symbols are present
    required_symbols = [
        "BeginPackage[\"UQFFKozima`\"]",
        "SigmaNSCm[omega_, n_]",
        "SigmaNBase[omega_]",
        "SigmaNDensity[rho_]",
        "FKozima[omega_, n_, Nn_, PhiPhonon_, FUBi_, FU_]",
        "FKozimaDamping[FNeutron_, FGW_]",
        "ReversalCondition[FU_, Nn_, sigma_, Phi_]",
        "VDSIntegral[omega_]",
        "VDSIntegralClosed[omega_]",
        "MEffBuoyancy[UgList_, OmegaG_, M_, dg_]",
        "KozimaConstants[]",
        "KozimaValidation[]",
        "EndPackage[]",
    ]

    all_found = True
    for sym in required_symbols:
        if sym in wl_content:
            print(f"  ✓ {sym}")
        else:
            print(f"  ✗ MISSING: {sym}")
            all_found = False

    assert all_found, "Missing required Mathematica symbols"

    # Verify constants are embedded
    const_checks = [
        (str(SSQ), "SSq"),
        (str(BETA_I), "beta_i"),
        (str(SIGMA_0), "sigma_0"),
        (str(F_LENR_THZ), "f_LENR_THz"),
        (str(N_LEVELS), "n_levels"),
    ]
    print()
    for val, name in const_checks:
        assert val in wl_content, f"Constant {name}={val} not found in .wl"
        print(f"  ✓ {name} = {val}")

    # Summary
    summary = kernel.generate_summary()
    print()
    print(f"Replaces: {summary['replaces']}")
    print(f"Provides: {len(summary['provides'])} functions/symbols")
    for p in summary["provides"]:
        print(f"  • {p}")

    print()
    print(f"Total lines: {lines}")
    print()
    print("ALL ASSERTIONS PASSED")
    print(json.dumps({
        "module": "kozima_wstp_kernel",
        "status": "VALIDATED",
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "wl_lines": lines,
        "symbols_verified": len(required_symbols),
    }, indent=2))


if __name__ == "__main__":
    main()
