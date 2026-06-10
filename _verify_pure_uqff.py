"""Verify pure-UQFF state of uqff_pure_calculator.py."""
import uqff_pure_calculator as u

print("=" * 70)
print("PURE-UQFF VERIFICATION")
print("=" * 70)

# 1. Confirm contamination is GONE.
removed = [
    '_l96_uqff_E_hartree_derived',
    '_l96_uqff_A_s_slow_roll',
    '_l96_uqff_A_s_PAPER1092_derived',
    '_l96_tier0_closure_probe',
    '_l96_uqff_higgs_self_coupling_lambda',
    '_l96_uqff_K_HIGGS_bridge',
    '_l96_uqff_m_higgs_PAPER639_GeV',
    '_l96_uqff_m_W_PAPER641_GeV',
    '_l96_uqff_sin2_theta_W_PAPER641',
    '_l96_uqff_CKM_Vus_PAPER634',
    '_l96_uqff_CKM_Vcb_PAPER634',
    '_l96_EW_higgs_pack_probe',
    '_l96_uqff_Lambda_cosmological_PAPER1156',
    '_l96_uqff_H0_dual_anchor',
    '_l96_uqff_t0_age_universe_Gyr',
    '_l96_uqff_T_neutrino_K',
    '_l96_uqff_T_CMB_K',
    '_l96_uqff_BAO_rd_Mpc',
    '_l96_cosmology_pack_probe',
    '_l96_uqff_eht_shadow_PAPER816_uas',
    '_l96_uqff_gw150914_ringdown_PAPER510_Hz',
    '_l96_uqff_gw170817_chirp_freq_PAPER675_Hz',
    '_l96_uqff_jades_z14_mass_PAPER1071_Msun',
    '_l96_multimessenger_pack_probe',
]
print("\nREMOVED (contamination):")
for f in removed:
    status = 'GONE' if not hasattr(u, f) else 'STILL PRESENT!!'
    flag = '  OK ' if status == 'GONE' else '  XX '
    print(f"{flag}{f:<55} {status}")

# 2. Confirm pure-UQFF Tier 2A survives.
kept = [
    '_l96_uqff_m_W_PAPER1209HH_GeV',
    '_l96_uqff_m_Z_PAPER1209HH_GeV',
    '_l96_uqff_m_t_PAPER1209HH_GeV',
    '_l96_uqff_m_H_PAPER1209HH_GeV',
    '_l96_uqff_m_b_PAPER1209HH_GeV',
    '_l96_uqff_m_c_PAPER1209HH_GeV',
    '_l96_uqff_m_tau_PAPER1209HH_GeV',
    '_l96_uqff_m_mu_PAPER1209HH_GeV',
    '_l96_uqff_m_s_PAPER1209HH_GeV',
    '_l96_uqff_m_e_PAPER1209HH_GeV',
    '_l96_uqff_proton_mass_PAPER1155_kg',
    '_l96_uqff_neutron_mass_PAPER1155_kg',
    '_l96_particle_mass_pack_probe',
    '_l96_uqff_omega_lambda_PAPER1156',
    '_l96_omega_lambda_probe',
]
print("\nKEPT (pure UQFF):")
for f in kept:
    status = 'OK' if hasattr(u, f) else 'MISSING!!'
    flag = '  OK ' if status == 'OK' else '  XX '
    print(f"{flag}{f:<55} {status}")

# 3. Run probes.
print("\n" + "=" * 70)
print("Particle-mass pack probe (PAPER_1209HH + PAPER_1155):")
print("=" * 70)
r = u._l96_particle_mass_pack_probe()
print(r['summary'])
print('all_within_5pct =', r['all_within_5pct'])

print("\n" + "=" * 70)
print("Omega_Lambda pure UQFF probe (PAPER_1156 §3):")
print("=" * 70)
r2 = u._l96_omega_lambda_probe()
print(r2['summary'])
