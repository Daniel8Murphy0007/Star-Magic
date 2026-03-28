#!/usr/bin/env wolframscript
(* ============================================================
   derive_uqff_constants.wls  v2
   UQFF Calibrated Constant Derivation via Wolfram Engine 14.3
   Session 151H — March 28, 2026

   NOTE: No underscores in variable names (Wolfram treats _ as Pattern)
   Derives/constrains: [SSq], kappa, Z26, [UA]
   ============================================================ *)

Print["==================================================="];
Print["  UQFF Constant Derivation - Wolfram Engine 14.3  "];
Print["===================================================\n"];

(* ---------------------------------------------------------- *)
(* PART 1: Z26 = PolyLog[26, SSq]                             *)
(* ---------------------------------------------------------- *)
Print["--- PART 1: Z26 = Li26(SSq) ---"];

SSqnom = 0.57;
Z26 = N[PolyLog[26, SSqnom], 15];
Print["Z26 = PolyLog[26, 0.57] = ", Z26];

Z26at4over7 = N[PolyLog[26, 4/7], 15];
Print["Z26 = PolyLog[26, 4/7] = ", Z26at4over7];
Print["Diff |0.57 - 4/7| = ", N[Abs[0.57 - 4/7], 8]];
Print["Diff |Z26(0.57) - Z26(4/7)| = ", N[Abs[Z26 - Z26at4over7], 10]];

(* Why Z26 ~ SSq: for large n, Li_n(x) ~ x because the k=1 term dominates *)
Print["\nKey insight: For large degree n, Li_n(x) = x + x^2/2^n + x^3/3^n + ..."];
Print["At n=26, second term = 0.57^2/2^26 = ", N[0.57^2/2^26, 4]];
Print["Hence Z26 ~ SSq is expected, not a special property of 0.57"];

(* ---------------------------------------------------------- *)
(* PART 2: [SSq] = 0.57 — where does it come from?           *)
(* ---------------------------------------------------------- *)
Print["\n--- PART 2: Origin of SSq = 0.57 ---"];

(* The triple convergence: CMB/Kepler/ALMA all point to 0.57 *)
(* From ARCHITECTURE_FLOW_DIAGRAM.md line 610: [SSq]=0.57=4/7 approx *)
(* From atomic_uqff_framework.py: calibrated from NS merger Ye~0.1 *)

(* Test: is the connection to Ye via the 26-layer formula? *)
(* beta scaling: beta_i * (1 - Ye)^(26/i) at i=1 gives (1-0.1)^26 *)
YeNS = 0.1;
betaFromYe = N[(1 - YeNS)^(1/26), 8];
Print["beta = (1-Ye)^(1/26) at Ye=0.1: ", betaFromYe];
Print["  -> 0.9^(1/26) = ", N[0.9^(1/26), 8], " not 0.57"];

(* Correct connection: SSq is the squeezing parameter in the density formula *)
(* rho_SCm(r) = rho0 * SSq^(r/r0) — multiplicative suppression per scale *)
(* At 26 layers: total suppression = SSq^26 *)
Print["SSq^26 = ", N[0.57^26, 6], " (26-layer total suppression)"];
Print["SSq^13 = ", N[0.57^13, 6], " (13-layer, half-shell suppression)"];

(* The 4/7 rational test: *)
Print["\n4/7 = ", N[4/7, 10]];
Print["Diff 0.57 - 4/7 = ", N[0.57 - 4/7, 6], " (0.57 is slightly less)"];
Print["FractionBox: SSq = 57/100 exactly (terminating decimal)"];
Print["The physics discriminates 57/100 vs 4/7 only at 4th decimal place"];

(* ---------------------------------------------------------- *)
(* PART 3: kappa = 1/tau derivation                           *)
(* ---------------------------------------------------------- *)
Print["\n--- PART 3: kappa from JWST quasar tau=2000 days ---"];

tauJWST = 2000;  (* days, from bodies.csv JWST_Quasar_Kappa row *)
kappaDerived = N[1/tauJWST, 6];
kappaTarget = 5 * 10^(-4);
Print["kappa = 1/tau = 1/", tauJWST, " days = ", ScientificForm[kappaDerived, 4]];
Print["Target kappa = 5e-4 = ", ScientificForm[kappaTarget, 4]];
Print["Match: kappa derived = kappa target? ",
      If[N[Abs[kappaDerived - kappaTarget], 10] < 10^(-15), "EXACT",
         "Diff = " <> ToString[N[Abs[kappaDerived - kappaTarget], 6]]]];

(* E-folding: at t=tau, E drops to 1/e *)
Print["E-fold check: exp(-kappa * tau) = exp(-1) = ", N[Exp[-1], 6]];

(* ---------------------------------------------------------- *)
(* PART 4: FindFit for SSq from log-spaced density data       *)
(* ---------------------------------------------------------- *)
Print["\n--- PART 4: FindFit SSq from Solar system SCm densities ---"];

(* bodies.csv-inspired data: Sun/Jupiter/Earth/Neptune rho_SCm *)
(* log10 values at n=0,1,2,3 spacing *)
logData = {{0, 15.0}, {1, 13.0}, {2, 12.0}, {3, 11.0}};
linFit = FindFit[logData, {aa + bb*nn, bb < 0}, {{aa, 15}, {bb, -1}}, nn];
slopeLog = N[bb /. linFit, 8];
intercLog = N[aa /. linFit, 8];
Print["Log10 linear fit: intercept=", intercLog, " slope=", slopeLog];
Print["slope = -log10(base) where base is the ratio per step"];
basePerStep = N[10^(-slopeLog), 6];
Print["Density ratio per step = 10^|slope| = ", basePerStep];

(* Natural-log exponential fit *)
expData = {{0, 1*^15}, {1, 1*^13}, {2, 1*^12}, {3, 1*^11}};
expFit = FindFit[expData, rhoA * Exp[-alphaE * nE], {{rhoA, 1*^15}, {alphaE, 1}}, nE];
alphaFit = N[alphaE /. expFit, 6];
rhoAFit = N[rhoA /. expFit, 4];
Print["Exp fit: alpha = ", alphaFit, " (compare SSq=0.57, but note alpha is per body-step)"];
Print["         rho0 = ", ScientificForm[rhoAFit]];
Print["Note: alpha here is per 1-unit n step, not directly SSq"];
Print["SSq enters as: rho(n) = rho0 * Exp[-SSq * n * C] where C = ln(rho_step)"];

(* ---------------------------------------------------------- *)
(* PART 5: [UA] constraint from H_SCm = 0.99                  *)
(* ---------------------------------------------------------- *)
Print["\n--- PART 5: [UA] from H_SCm = 0.99 constraint ---"];

HScm = 0.99;
(* H_SCm = rhoSCm / (rhoSCm + rhoUA) = 0.99 *)
(* => rhoUA = rhoSCm * (1 - HScm) / HScm *)
ratioUAtoSCm = N[(1 - HScm)/HScm, 6];
Print["rhoUA/rhoSCm = (1 - H_SCm)/H_SCm = ", ratioUAtoSCm];
Print["[UA] enters as rhoUA = [UA] * rho_vac_lambda"];
Print["[UA] = 1e-4 from bodies.csv column U_UA (all rows)"];
Print["Interpreted as: Boson-aether density = 1e-4 of total vacuum density"];

(* ---------------------------------------------------------- *)
(* PART 6: Cross-check — all constants together               *)
(* ---------------------------------------------------------- *)
Print["\n--- PART 6: Cross-verification ---"];

(* UQFF Buoyancy force: F_Ubi = [UA] * SSq * kappa * Z26 *)
FUbi = N[1*^-4 * 0.57 * 5*^-4 * Z26, 6];
Print["F_Ubi = [UA]*SSq*kappa*Z26 = ", ScientificForm[FUbi]];
Print["This is the characteristic buoyancy coupling strength"];

(* The 26-layer exponential at SSq=0.57: full shell suppression *)
fullSuppression = N[0.57^26, 8];
Print["SSq^26 = ", fullSuppression, " (full 26-layer suppression factor)"];
Print["sqrt(SSq^26) = SSq^13 = ", N[0.57^13, 6], " (half-stack)"];

(* kappa Z26 product: the decay-weighted polylog *)
kappaZ26 = N[5*^-4 * Z26, 8];
Print["kappa * Z26 = ", ScientificForm[kappaZ26], " day^-1"];

Print["\n=== FINAL RESULTS ==="];
Print["[SSq]  = 0.57 = 57/100  [empirical; ~4/7 but not exact; triple convergence]"];
Print["Z26    = ", Z26, "  [= PolyLog[26, 0.57]; ~SSq by Li_n large-n property]"];
Print["kappa  = ", ScientificForm[kappaDerived, 3], " day^-1  [= 1/2000 exactly; JWST quasar tau]"];
Print["[UA]   = 1e-4           [vacuum density fraction; H_SCm constraint gives ratio only]"];

Print["\n==================================================="];
Print["  DERIVATION COMPLETE — v2                        "];
Print["==================================================="];
