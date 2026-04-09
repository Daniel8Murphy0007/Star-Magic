(* ::Package:: *)
(* UQFFMockThetaPi` — Mock Theta Functions & Ramanujan 1/π for UQFF *)
(* Session 204 | Daniel Murphy *)
(* Generated: 2026-04-09 16:30:23 UTC *)
(* Import: << "uqff_mock_theta_pi_kernel.wl" *)

BeginPackage["UQFFMockThetaPi`"];

(* ── Public symbols ─────────────────────────────────────────────── *)
qPochhammer::usage =
  "qPochhammer[a, q, n] computes the q-Pochhammer symbol \
   (a; q)_n = Product[1 - a q^k, {k, 0, n-1}].";

f26::usage =
  "f26[q] computes the third-order Ramanujan mock theta function \
   f_26(q) = Sum[q^(n^2) / QPochhammer[-q, q, n]^2, {n, 0, 25}].";

phi26::usage =
  "phi26[q] computes the third-order mock theta function \
   phi_26(q) = Sum[q^(n^2) / QPochhammer[-q^2, q^2, n], {n, 0, 25}].";

psi26::usage =
  "psi26[q] computes the third-order mock theta function \
   psi_26(q) = Sum[q^(n^2) / QPochhammer[q, q^2, n], {n, 1, 26}].";

thetaCoupled26::usage =
  "thetaCoupled26[q, ssq, kap] computes the 26-state UQFF coupled \
   theta amplitude Theta_26 = Sum[A_i (f26[q_i]+phi26[q_i]+psi26[q_i])].";

ramanujanR::usage =
  "ramanujanR[n] computes the Ramanujan coefficient \
   R_n = (4n)! / ((n!)^4 * 396^(4n)).";

oneOverPiClassical::usage =
  "oneOverPiClassical[nTerms] computes Ramanujan's 1/pi series \
   (2 Sqrt[2]/9801) Sum[R_n (1103 + 26390 n), {n, 0, nTerms-1}].";

oneOverPiUQFF::usage =
  "oneOverPiUQFF[nTerms, ssq, kap] computes UQFF-modified 1/pi \
   with 26-state weight W_26(n) = Product[1+ssq Exp[-kap i n/26]].";

pi26DHypergeometric::usage =
  "pi26DHypergeometric[nTerms] computes 26D-generalized Ramanujan 1/pi \
   with a_26 and b_26 modified coefficients, C_26 normalized.";

Begin["`Private`"];

(* ── SSQ = 0.57, kappa = 5.787e-09, N_levels = 26 ──── *)
ssqDefault = 0.57;
kappaDefault = 5.787e-09;
nLevels = 26;

(* ── §1 q-Pochhammer ──────────────────────────────────────── *)
qPochhammer[a_, q_, 0] := 1;
qPochhammer[a_, q_, n_Integer /; n > 0] :=
  Product[1 - a q^k, {k, 0, n - 1}];

(* ── §2 Mock theta functions ──────────────────────────────── *)
f26[q_] := Module[{total = 0, denom},
  Do[
    denom = qPochhammer[-q, q, n]^2;
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {n, 0, nLevels - 1}];
  total];

phi26[q_] := Module[{total = 0, denom},
  Do[
    denom = qPochhammer[-q^2, q^2, n];
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {n, 0, nLevels - 1}];
  total];

psi26[q_] := Module[{total = 0, denom},
  Do[
    denom = qPochhammer[q, q^2, n];
    If[Abs[denom] > 10^(-300),
      total += q^(n^2) / denom],
    {n, 1, nLevels}];
  total];

(* ── §3 UQFF coupled theta ───────────────────────────────── *)
vdsAmplitude[i_] := (2 Pi)^(i/6) * (7.09*^-37/7.09*^-36);

thetaCoupled26[q_, ssq_:ssqDefault, kap_:kappaDefault] :=
  Module[{total = 0, qi, ai},
    Do[
      qi = ssq * Exp[-kap * i * 0 / nLevels];  (* t=0 default *)
      ai = vdsAmplitude[i];
      total += ai * (f26[qi] + phi26[qi] + psi26[qi]),
      {i, 1, nLevels}];
    total];

(* ── §4 Ramanujan R_n coefficient ─────────────────────────── *)
ramanujanR[n_Integer] := (4 n)! / ((n!)^4 * 396^(4 n));

(* Numerical version for large n *)
ramanujanRN[n_Integer] :=
  Exp[LogGamma[4 n + 1] - 4 LogGamma[n + 1] - 4 n Log[396]];

(* ── §5 Classical Ramanujan 1/pi ──────────────────────────── *)
oneOverPiClassical[nTerms_Integer:26] :=
  (2 Sqrt[2] / 9801) *
  Sum[ramanujanR[n] (1103 + 26390 n), {n, 0, nTerms - 1}];

(* ── §6 UQFF-modified 1/pi ───────────────────────────────── *)
weight26[n_, ssq_, kap_] :=
  Product[1 + ssq Exp[-kap i n / nLevels], {i, 1, nLevels}];

c26Norm[ssq_] := (1 + ssq)^nLevels;

oneOverPiUQFF[nTerms_Integer:26, ssq_:ssqDefault, kap_:kappaDefault] :=
  (2 Sqrt[2] / 9801) / c26Norm[ssq] *
  Sum[ramanujanR[n] (1103 + 26390 n) weight26[n, ssq, kap],
      {n, 0, nTerms - 1}];

(* ── §7 26D hypergeometric 1/pi ───────────────────────────── *)
h26Alt = Sum[(-1)^(k + 1) / k, {k, 1, nLevels}];
a26 = 1103 * h26Alt;
b26 = 26390 * (nLevels / 13);
c26HyperNorm = h26Alt * (nLevels / 13);

pi26DHypergeometric[nTerms_Integer:26] :=
  1 / ((2 Sqrt[2] / 9801) / c26HyperNorm *
       Sum[ramanujanR[n] (a26 + b26 n), {n, 0, nTerms - 1}]);

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
