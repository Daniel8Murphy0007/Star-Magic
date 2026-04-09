"""
Mock Theta WSTP Kernel — f₂₆(q), φ₂₆(q), ψ₂₆(q) + oneOverPiUQFF Wolfram Export

Session 204 | Daniel Murphy
PURPOSE: The existing s26_wstp_kernel.py exports 8 polylogarithm symbols
         to Wolfram Language (S26, R26, NaiveLi, S26VDS, etc.).
         NONE include mock theta functions or Ramanujan 1/π series.

         This module generates a companion .wl package exporting:
           1. qPochhammer[a, q, n]        — q-Pochhammer symbol (a;q)_n
           2. f26[q]                       — third-order mock theta f₂₆(q)
           3. phi26[q]                     — third-order mock theta φ₂₆(q)
           4. psi26[q]                     — third-order mock theta ψ₂₆(q)
           5. thetaCoupled26[q, ssq, kap]  — 26-state Θ₂₆ amplitude
           6. ramanujanR[n]                — R_n coefficient
           7. oneOverPiClassical[nTerms]   — Ramanujan 1/π
           8. oneOverPiUQFF[nTerms, ssq, kap] — UQFF-modified 1/π
           9. pi26DHypergeometric[nTerms]  — 26D hypergeometric 1/π

         Import via: << "uqff_mock_theta_pi_kernel.wl"

ARCHITECTURE: Pure .wl generator. No hardcoded astrophysical systems.
"""

import json
import math
import os
from datetime import datetime, timezone
from typing import Dict, Optional

# ── §0  CONSTANTS ─────────────────────────────────────────────────────────

PI = math.pi
SSQ = 0.57
KAPPA = 5.787e-9
N_LEVELS = 26

# Ramanujan series
RAMANUJAN_A = 1103
RAMANUJAN_B = 26390
RAMANUJAN_C = 396


# ── §1  WOLFRAM PACKAGE TEMPLATE ─────────────────────────────────────────

WL_PACKAGE = r"""(* ::Package:: *)
(* UQFFMockThetaPi` — Mock Theta Functions & Ramanujan 1/π for UQFF *)
(* Session 204 | Daniel Murphy *)
(* Generated: {timestamp} *)
(* Import: << "uqff_mock_theta_pi_kernel.wl" *)

BeginPackage["UQFFMockThetaPi`"];

(* ── Public symbols ─────────────────────────────────────────────── *)
qPochhammer::usage =
  "qPochhammer[a, q, n] computes the q-Pochhammer symbol \
   (a; q)_n = Product[1 - a q^k, {{k, 0, n-1}}].";

f26::usage =
  "f26[q] computes the third-order Ramanujan mock theta function \
   f_26(q) = Sum[q^(n^2) / QPochhammer[-q, q, n]^2, {{n, 0, 25}}].";

phi26::usage =
  "phi26[q] computes the third-order mock theta function \
   phi_26(q) = Sum[q^(n^2) / QPochhammer[-q^2, q^2, n], {{n, 0, 25}}].";

psi26::usage =
  "psi26[q] computes the third-order mock theta function \
   psi_26(q) = Sum[q^(n^2) / QPochhammer[q, q^2, n], {{n, 1, 26}}].";

thetaCoupled26::usage =
  "thetaCoupled26[q, ssq, kap] computes the 26-state UQFF coupled \
   theta amplitude Theta_26 = Sum[A_i (f26[q_i]+phi26[q_i]+psi26[q_i])].";

ramanujanR::usage =
  "ramanujanR[n] computes the Ramanujan coefficient \
   R_n = (4n)! / ((n!)^4 * 396^(4n)).";

oneOverPiClassical::usage =
  "oneOverPiClassical[nTerms] computes Ramanujan's 1/pi series \
   (2 Sqrt[2]/9801) Sum[R_n (1103 + 26390 n), {{n, 0, nTerms-1}}].";

oneOverPiUQFF::usage =
  "oneOverPiUQFF[nTerms, ssq, kap] computes UQFF-modified 1/pi \
   with 26-state weight W_26(n) = Product[1+ssq Exp[-kap i n/26]].";

pi26DHypergeometric::usage =
  "pi26DHypergeometric[nTerms] computes 26D-generalized Ramanujan 1/pi \
   with a_26 and b_26 modified coefficients, C_26 normalized.";

Begin["`Private`"];

(* ── SSQ = {ssq}, kappa = {kappa}, N_levels = {n_levels} ──── *)
ssqDefault = {ssq};
kappaDefault = {kappa};
nLevels = {n_levels};

(* ── §1 q-Pochhammer ──────────────────────────────────────── *)
qPochhammer[a_, q_, 0] := 1;
qPochhammer[a_, q_, n_Integer /; n > 0] :=
  Product[1 - a q^k, {{k, 0, n - 1}}];

(* ── §2 Mock theta functions ──────────────────────────────── *)
f26[q_] := Module[{{total = 0, denom}},
  Do[
    denom = qPochhammer[-q, q, n]^2;
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {{n, 0, nLevels - 1}}];
  total];

phi26[q_] := Module[{{total = 0, denom}},
  Do[
    denom = qPochhammer[-q^2, q^2, n];
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {{n, 0, nLevels - 1}}];
  total];

psi26[q_] := Module[{{total = 0, denom}},
  Do[
    denom = qPochhammer[q, q^2, n];
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {{n, 1, nLevels}}];
  total];

(* ── §3 UQFF coupled theta ───────────────────────────────── *)
vdsAmplitude[i_] := (2 Pi)^(i/6) * ({rho_scm}/{rho_ua});

thetaCoupled26[q_, ssq_:ssqDefault, kap_:kappaDefault] :=
  Module[{{total = 0, qi, ai}},
    Do[
      qi = ssq * Exp[-kap * i * 0 / nLevels];  (* t=0 default *)
      ai = vdsAmplitude[i];
      total += ai * (f26[qi] + phi26[qi] + psi26[qi]),
      {{i, 1, nLevels}}];
    total];

(* ── §4 Ramanujan R_n coefficient ─────────────────────────── *)
ramanujanR[n_Integer] := (4 n)! / ((n!)^4 * 396^(4 n));

(* Numerical version for large n *)
ramanujanRN[n_Integer] :=
  Exp[LogGamma[4 n + 1] - 4 LogGamma[n + 1] - 4 n Log[396]];

(* ── §5 Classical Ramanujan 1/pi ──────────────────────────── *)
oneOverPiClassical[nTerms_Integer:26] :=
  (2 Sqrt[2] / 9801) *
  Sum[ramanujanR[n] (1103 + 26390 n), {{n, 0, nTerms - 1}}];

(* ── §6 UQFF-modified 1/pi ───────────────────────────────── *)
weight26[n_, ssq_, kap_] :=
  Product[1 + ssq Exp[-kap i n / nLevels], {{i, 1, nLevels}}];

c26Norm[ssq_] := (1 + ssq)^nLevels;

oneOverPiUQFF[nTerms_Integer:26, ssq_:ssqDefault, kap_:kappaDefault] :=
  (2 Sqrt[2] / 9801) / c26Norm[ssq] *
  Sum[ramanujanR[n] (1103 + 26390 n) weight26[n, ssq, kap],
      {{n, 0, nTerms - 1}}];

(* ── §7 26D hypergeometric 1/pi ───────────────────────────── *)
h26Alt = Sum[(-1)^(k + 1) / k, {{k, 1, nLevels}}];
a26 = 1103 * h26Alt;
b26 = 26390 * (nLevels / 13);
c26HyperNorm = h26Alt * (nLevels / 13);

pi26DHypergeometric[nTerms_Integer:26] :=
  1 / ((2 Sqrt[2] / 9801) / c26HyperNorm *
       Sum[ramanujanR[n] (a26 + b26 n), {{n, 0, nTerms - 1}}]);

End[];
EndPackage[];

(* ── Self-test ────────────────────────────────────────────── *)
Print["UQFFMockThetaPi` loaded."];
Print["  f26[0.57]   = ", N[f26[0.57], 15]];
Print["  phi26[0.57] = ", N[phi26[0.57], 15]];
Print["  psi26[0.57] = ", N[psi26[0.57], 15]];
Print["  1/Pi class  = ", N[oneOverPiClassical[4], 40]];
Print["  Pi approx   = ", N[1/oneOverPiClassical[4], 40]];
Print["  Pi ref      = ", N[Pi, 40]];
"""


# ── §2  PACKAGE GENERATOR ────────────────────────────────────────────────

class MockThetaPiWSTPKernel:
    """
    Generates the .wl Wolfram Language package for mock theta functions
    and Ramanujan 1/π series with UQFF 26-state modifications.
    """

    def __init__(self, output_dir: str = "."):
        self.output_dir = output_dir
        self.filename = "uqff_mock_theta_pi_kernel.wl"

    def generate(self) -> Dict[str, str]:
        """Generate the .wl package string."""
        now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

        content = WL_PACKAGE.format(
            timestamp=now,
            ssq=SSQ,
            kappa=KAPPA,
            n_levels=N_LEVELS,
            rho_scm="7.09*^-37",
            rho_ua="7.09*^-36",
        )

        return {
            "filename": self.filename,
            "content": content,
            "symbols": [
                "qPochhammer", "f26", "phi26", "psi26",
                "thetaCoupled26", "ramanujanR", "oneOverPiClassical",
                "oneOverPiUQFF", "pi26DHypergeometric",
            ],
            "n_symbols": 9,
            "package": "UQFFMockThetaPi`",
        }

    def write(self) -> str:
        """Write .wl file to disk and return path."""
        result = self.generate()
        path = os.path.join(self.output_dir, self.filename)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(result["content"])
        return path


# ── §3  SELF-TEST ─────────────────────────────────────────────────────────

def self_test() -> Dict:
    """Validate WSTP kernel generation."""
    results = {}

    gen = MockThetaPiWSTPKernel()
    pkg = gen.generate()

    # Check structure
    assert pkg["n_symbols"] == 9, f"Expected 9 symbols, got {pkg['n_symbols']}"
    results["symbol_count"] = "PASS"

    content = pkg["content"]

    # Check all 9 symbols are present in the package
    for sym in pkg["symbols"]:
        assert sym in content, f"Symbol {sym} missing from package"
    results["all_symbols_present"] = "PASS"

    # Check package structure
    assert "BeginPackage" in content, "Missing BeginPackage"
    assert "EndPackage" in content, "Missing EndPackage"
    assert 'UQFFMockThetaPi`' in content, "Wrong package name"
    results["package_structure"] = "PASS"

    # Check key formulas are present
    assert "9801" in content, "Missing 9801 denominator"
    assert "1103" in content, "Missing 1103 constant"
    assert "26390" in content, "Missing 26390 constant"
    assert "396" in content, "Missing 396 base"
    assert "LogGamma" in content, "Missing LogGamma for overflow safety"
    results["formulas_present"] = "PASS"

    # Check line count is reasonable
    lines = content.strip().split("\n")
    assert len(lines) >= 80, f"Only {len(lines)} lines — too short"
    assert len(lines) <= 300, f"{len(lines)} lines — suspiciously long"
    results["line_count"] = len(lines)

    # Write test
    test_path = gen.write()
    assert os.path.exists(test_path), f"File not written: {test_path}"
    file_size = os.path.getsize(test_path)
    assert file_size > 2000, f"File too small: {file_size} bytes"
    results["file_written"] = {"path": test_path, "size_bytes": file_size}

    results["overall"] = "PASS"
    return results


if __name__ == "__main__":
    print("=" * 72)
    print("Mock Theta Pi WSTP Kernel — Self-Test")
    print("=" * 72)
    r = self_test()
    print(json.dumps(r, indent=2, default=str))
    print("=" * 72)
    print(f"Overall: {r['overall']}")
