"""
Chain trace -- Fix #5 / Fix #6 / Fix #7 / Fix #9 / Fix #10
Exit 0 if all checks pass, non-zero otherwise.
"""
import sys
import math
import dpm_vacuum_manifold as dpm

PASS = 0
FAIL = 0

def check(label: str, cond: bool, detail: str = "") -> None:
    global PASS, FAIL
    status = "[PASS]" if cond else "[FAIL]"
    print(f"  {status}  {label}  {detail}")
    if cond:
        PASS += 1
    else:
        FAIL += 1
        raise AssertionError(f"Check failed: {label}  {detail}")

# ---------------------------------------------------------------------------
print("=" * 78)
print("  SECTION 1 — Fix #5: ACP denominator proof (S5h)")
print("=" * 78)

r5 = dpm.chain_acp_denominator_proof()
print(f"  rho_SCm  = {r5['rho_SCm']:.2e} kg/m^3")
print(f"  rho_UA   = {r5['rho_UA']:.2e} kg/m^3")
print(f"  DPM ratio = rho_UA/rho_SCm = {r5['dpm_ratio_from_vacuum']}")
print(f"  Denominator canonical = {r5['denominator_canonical']}")
print(f"  M fraction at 1 cycle = {r5['mass_at_1_cycle_frac']:.4f}")
print(f"  Half-mass n_grad = {r5['half_mass_n_grad']:.4f}")
print(f"  Full ladder (10 cycles) = {r5['full_ladder_10_cycles']:.6f}")

check("DPM ratio == 10",
      abs(r5["dpm_ratio_from_vacuum"] - 10.0) < 1e-10,
      f"ratio={r5['dpm_ratio_from_vacuum']}")

check("denominator_canonical == DPM_RATIO",
      abs(r5["error_denominator_vs_ratio"]) < 1e-10,
      f"err={r5['error_denominator_vs_ratio']}")

check("mass_at_1_cycle_frac == 1-1/e",
      abs(r5["mass_at_1_cycle_frac"] - (1.0 - 1.0/math.e)) < 1e-6,
      f"frac={r5['mass_at_1_cycle_frac']:.6f}")

check("half_mass_n_grad == DPM_RATIO*ln2",
      abs(r5["half_mass_n_grad"] - 10.0 * math.log(2.0)) < 1e-6,
      f"n={r5['half_mass_n_grad']:.4f}")

check("full_ladder_10_cycles == DPM_RATIO (within 1e-6)",
      abs(r5["full_ladder_10_cycles"] - 10.0) < 1e-6,
      f"ladder={r5['full_ladder_10_cycles']:.6f}")

check("saturation_table has 20 entries",
      len(r5["saturation_table"]) == 20,
      f"len={len(r5['saturation_table'])}")

check("saturation monotone increasing",
      all(r5["saturation_table"][i]["M_frac"] < r5["saturation_table"][i+1]["M_frac"]
          for i in range(19)),
      "")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  SECTION 2 — Fix #6: DPM ratio scale invariance (S5i)")
print("=" * 78)

r6 = dpm.chain_dpm_ratio_scale_invariance()
print(f"  DPM_RATIO = {r6['dpm_ratio']}")
print(f"  Layer base == DPM_RATIO: {r6['layer_base_check']['base_equals_dpm_ratio']}")
print()
print(f"  {'Scale':<12} {'r [m]':>12}  {'DPM ratio':>12}  {'E_react J/m3':>14}")
print(f"  {'-'*56}")
for sc in r6["scale_checks"]:
    print(f"  {sc['scale']:<12} {sc['r_m']:>12.3e}  {sc['dpm_ratio_at_scale']:>12.1f}  {sc['E_react_J_m3']:>14.3e}")

check("DPM_RATIO == 10 at all 4 scales",
      all(abs(sc["dpm_ratio_at_scale"] - 10.0) < 1e-10 for sc in r6["scale_checks"]),
      "")

check("All scale ratio errors == 0",
      all(sc["ratio_error"] < 1e-10 for sc in r6["scale_checks"]),
      "")

check("Layer base 10 equals DPM_RATIO",
      r6["layer_base_check"]["base_equals_dpm_ratio"],
      f"base={r6['layer_base_check']['E_n_formula_base']}")

check("4 scale checks present",
      len(r6["scale_checks"]) == 4,
      f"n={len(r6['scale_checks'])}")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  SECTION 3 — Fix #7: Falsifiable predictions (S5j)")
print("=" * 78)

r7 = dpm.chain_falsifiable_predictions()
print()
for key, pred in r7.items():
    if "predicted_value" in pred:
        print(f"  {key}: {pred['predicted_value']:.4e} {pred.get('unit','')}")
    elif "predicted_r_cross_m" in pred:
        print(f"  {key}: r_cross = {pred['predicted_r_cross_m']:.4e} m  E_thr = {pred['E_threshold_keV']:.2f} keV")
    elif "predicted_E_MeV" in pred:
        print(f"  {key}: E = {pred['predicted_E_MeV']:.1f} MeV")
    elif "predicted_E_gap_MeV" in pred:
        print(f"  {key}: E_gap = {pred['predicted_E_gap_MeV']:.1f} MeV")
    elif "predicted_scaling" in pred:
        print(f"  {key}: {pred['predicted_scaling']}  err={pred['scaling_error_pct']:.4f}%")
    elif "predicted_value_MeV" in pred:
        print(f"  {key}: {pred['predicted_value_MeV']:.5f} MeV  (PDG: {pred['pdg_value_MeV']:.5f} MeV)")

check("6 predictions present",
      len(r7) == 6,
      f"n={len(r7)}")

check("P1 electron radius in (1e-14, 1e-12) m",
      1e-14 < r7["P1_electron_confinement_radius"]["predicted_value"] < 1e-12,
      f"r={r7['P1_electron_confinement_radius']['predicted_value']:.4e}")

check("P1 test threshold > 0.1 MeV/c (sub-MeV electron scattering)",
      r7["P1_electron_confinement_radius"]["test_threshold_MeV_c"] > 0.1,
      f"q={r7['P1_electron_confinement_radius']['test_threshold_MeV_c']:.3f}")

check("P3 r_cross positive",
      r7["P3_r_cross_Z1_resonance"]["predicted_r_cross_m"] > 0,
      f"r={r7['P3_r_cross_Z1_resonance']['predicted_r_cross_m']:.4e}")

check("P3 E_threshold > 0 keV",
      r7["P3_r_cross_Z1_resonance"]["E_threshold_keV"] > 0,
      f"E={r7['P3_r_cross_Z1_resonance']['E_threshold_keV']:.2f}")

check("P4 Z^(-2/3) scaling error < 1%",
      r7["P4_r_cross_Z_scaling"]["scaling_error_pct"] < 1.0,
      f"err={r7['P4_r_cross_Z_scaling']['scaling_error_pct']:.4f}%")

check("P5 layer-13 energy > 100 GeV (electroweak/LHC scale)",
      r7["P5_layer13_threshold"]["predicted_E_MeV"] > 1e5,
      f"E={r7['P5_layer13_threshold']['predicted_E_MeV']:.1f} MeV = {r7['P5_layer13_threshold']['predicted_E_MeV']/1e3:.1f} GeV")

check("P6 E_crack > 0 (positive mass gap -- Yang-Mills criterion)",
      r7["P6_yang_mills_mass_gap"]["predicted_E_gap_MeV"] > 0.0,
      f"E_gap={r7['P6_yang_mills_mass_gap']['predicted_E_gap_MeV']:.4e} MeV = {r7['P6_yang_mills_mass_gap']['E_crack_J']:.3e} J")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  SECTION 4 — Fix #9: Derive rho_A (S5k)")
print("=" * 78)

r9 = dpm.chain_rhoA_derivation()
print(f"  rho_SCm              = {r9['rho_SCm']:.2e} kg/m^3")
print(f"  DPM_RATIO            = {r9['DPM_RATIO']}")
print(f"  Layer midpoint       = {r9['layer_midpoint']}")
print(f"  rho_A raw (layer 13) = {r9['rho_A_raw_derived']:.4e} kg/m^3")
print(f"  rho_A [SSq]-gated    = {r9['rho_A_ssq_gate']:.4e} kg/m^3")
print(f"  rho_A canonical      = {r9['rho_A_canonical']:.1e} kg/m^3")
print(f"  error (raw)          = {r9['error_raw_pct']:+.1f}%")
print(f"  error ([SSq]-gated)  = {r9['error_ssq_gate_pct']:+.1f}%")
print(f"  quantum coupling     = {r9['quantum_coupling']:.4e}")
print(f"  quantum viscosity    = {r9['quantum_viscosity_Pa_s']:.4e} Pa*s")

check("rho_A raw is rho_SCm * 10^13",
      abs(r9["rho_A_raw_derived"] - r9["rho_SCm"] * 10**13) / r9["rho_A_raw_derived"] < 1e-10,
      f"rho_A={r9['rho_A_raw_derived']:.4e}")

check("rho_A [SSq]-gated within 50% of canonical 1e-23",
      abs(r9["error_ssq_gate_pct"]) < 50.0,
      f"error={r9['error_ssq_gate_pct']:+.1f}%")

check("rho_A_ssq_gate in (1e-24, 1e-22) kg/m^3",
      1e-24 < r9["rho_A_ssq_gate"] < 1e-22,
      f"rho={r9['rho_A_ssq_gate']:.4e}")

check("Layer midpoint == 13",
      r9["layer_midpoint"] == 13,
      f"mid={r9['layer_midpoint']}")

check("26 layer densities computed",
      len(r9["layer_densities"]) == 26,
      f"n={len(r9['layer_densities'])}")

check("Layer densities strictly increasing",
      all(r9["layer_densities"][i]["rho_kg_m3"] < r9["layer_densities"][i+1]["rho_kg_m3"]
          for i in range(25)),
      "")

check("quantum viscosity < 1e-40 Pa*s (quasi-inviscid)",
      r9["quantum_viscosity_Pa_s"] < 1e-40,
      f"mu={r9['quantum_viscosity_Pa_s']:.4e}")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  SECTION 5 — Fix #10: B0_i confinement correction (S5l)")
print("=" * 78)

for Z_test in [1, 2, 26]:
    r10 = dpm.chain_b0_confinement_correction(Z_test)
    print(f"\n  Z={Z_test:3d}  R_nuc={r10['R_nuc_m']:.3e} m  i_crit={r10['i_crit']:.3f}"
          f"  n_conf={r10['n_layers_in_confinement']}  "
          f"A_26_std={r10['A_26_standard']:.4e}  A_26_corr={r10['A_26_corrected']:.4e}"
          f"  f_corr={r10['correction_factor']:.4f}")

# Z=1 checks
r10_1 = dpm.chain_b0_confinement_correction(1)
check("Z=1: A_26_standard == sum(i^6, i=1..26)",
      abs(r10_1["A_26_standard"] - sum(i**6 for i in range(1, 27))) < 1.0,
      f"A={r10_1['A_26_standard']:.0f}")

check("Z=1: ALL 26 layers in confinement (R_nuc < R_C_UP)",
      r10_1["n_layers_in_confinement"] == 26,
      f"n_conf={r10_1['n_layers_in_confinement']}")

check("Z=1: correction_factor < 1 (confinement reduces amplification)",
      r10_1["correction_factor"] < 1.0,
      f"f={r10_1['correction_factor']:.4f}")

# Z=26 checks
r10_26 = dpm.chain_b0_confinement_correction(26)
check("Z=26: at least 2 dipole-regime layers",
      r10_26["n_layers_in_dipole"] >= 2,
      f"n_dip={r10_26['n_layers_in_dipole']}")

check("Z=26: A_26_corrected > 0",
      r10_26["A_26_corrected"] > 0,
      f"A={r10_26['A_26_corrected']:.4e}")

check("Z=26: correction_factor < 1",
      r10_26["correction_factor"] < 1.0,
      f"f={r10_26['correction_factor']:.4f}")

# Z-dependence: i_crit grows with Z
r10_2 = dpm.chain_b0_confinement_correction(2)
check("i_crit(Z) increasing: i_crit(1) < i_crit(2) < i_crit(26)",
      r10_1["i_crit"] < r10_2["i_crit"] < r10_26["i_crit"],
      f"i_crit: {r10_1['i_crit']:.3f} < {r10_2['i_crit']:.3f} < {r10_26['i_crit']:.3f}")

check("correction_factor(Z=1) < correction_factor(Z=26) (more confinement at Z=1)",
      r10_1["correction_factor"] < r10_26["correction_factor"],
      f"f(Z=1)={r10_1['correction_factor']:.4f} vs f(Z=26)={r10_26['correction_factor']:.4f}")

# ---------------------------------------------------------------------------
print()
print("=" * 78)
print("  SUMMARY — Fix #5/#6/#7/#9/#10 chain trace")
print("=" * 78)

print(f"  Fix #5  ACP denominator = DPM_RATIO = {r5['dpm_ratio_from_vacuum']}  (error: {r5['error_denominator_vs_ratio']:.0e})")
print(f"  Fix #6  DPM ratio invariant across 4 scales (nuclear/atomic/stellar/galactic)")
print(f"  Fix #7  {len(r7)} falsifiable predictions generated")
print(f"  Fix #9  rho_A derived = {r9['rho_A_ssq_gate']:.3e} kg/m^3  ({r9['error_ssq_gate_pct']:+.1f}% vs 1e-23)")
print(f"  Fix #10 B0_i confinement: Z=1 f_corr={r10_1['correction_factor']:.4f}, Z=26 f_corr={r10_26['correction_factor']:.4f}")
print()
print(f"  All checks PASSED.  Exit code 0.")
sys.exit(0)
