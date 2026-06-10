"""Scan all _*_primitive_sat() functions in uqff_pure_calculator.

For each, call it and compare against a hardcoded set of CODATA / observation
anchors. Surface any primitive whose numerical value bears no relation to a
recognisable physics observable — i.e. dimensionless structural products with
physical-unit names. Same pattern as LENR _lenr_*_primitive_sat placeholders.
"""
import uqff_pure_calculator as u
import inspect
import math

# Anchors with rough order-of-magnitude tolerance.
ANCHORS = {
    # Particle masses (MeV)
    '_m_mu_primitive_sat':         (105.66,     'MeV', 'muon mass'),
    '_m_tau_primitive_sat':        (1776.86,    'MeV', 'tau mass'),
    '_m_t_primitive_sat':          (172500.0,   'MeV', 'top quark mass'),
    '_m_W_primitive_sat':          (80377.0,    'MeV', 'W boson mass'),
    '_m_Z_primitive_sat':          (91187.6,    'MeV', 'Z boson mass'),
    '_m_H_primitive_sat':          (125250.0,   'MeV', 'Higgs mass'),
    '_v_higgs_primitive_sat':      (246220.0,   'MeV', 'Higgs VEV'),
    '_m_e_primitive_sat':          (0.511,      'MeV', 'electron mass'),
    '_m_pion_primitive_sat':       (139.57,     'MeV', 'pion mass'),
    '_m_kaon_primitive_sat':       (493.677,    'MeV', 'kaon mass'),
    '_m_u_primitive_sat':          (2.16,       'MeV', 'up quark mass'),
    '_m_d_primitive_sat':          (4.67,       'MeV', 'down quark mass'),
    '_m_s_primitive_sat':          (93.4,       'MeV', 'strange quark mass'),
    '_m_c_primitive_sat':          (1270.0,     'MeV', 'charm quark mass'),
    '_m_b_primitive_sat':          (4180.0,     'MeV', 'bottom quark mass'),
    # Electroweak
    '_G_F_primitive_sat':          (1.1663787e-5, 'GeV^-2', 'Fermi constant'),
    '_alpha_s_primitive_sat':      (0.1179,     '', 'strong coupling at M_Z'),
    '_alpha_primitive_sat':        (1.0/137.035999084, '', 'fine-structure constant'),
    '_sin2_theta_w_primitive_sat': (0.23122,    '', 'sin^2 theta_W'),
    '_ckm_vus_primitive_sat':      (0.2243,     '', 'CKM Vus'),
    '_ckm_vcb_primitive_sat':      (0.0410,     '', 'CKM Vcb'),
    '_pmns_theta12_primitive_sat': (33.44,      'deg', 'PMNS theta12'),
    # Anomalous moments
    '_a_e_primitive_sat':          (1.15965218073e-3, '', 'electron g-2/2'),
    '_a_mu_primitive_sat':         (1.16592e-3, '', 'muon g-2/2'),
    '_g_e_primitive_sat':          (2.00231930436, '', 'electron g-factor'),
    # Particle / SI extras
    '_proton_mass_primitive_sat':  (938.272,    'MeV', 'proton mass'),
    '_neutron_lifetime_primitive_sat': (879.4,  's',   'neutron lifetime'),
    '_yang_mills_primitive_sat':   (1.78,       'GeV', 'YM mass gap'),
    '_E_hartree_primitive_sat':    (4.3597447222071e-18, 'J', 'Hartree energy'),
    '_hyperfine_cs_primitive_sat': (9.192631770e9, 'Hz', 'Cs hyperfine'),
    '_gas_R_primitive_sat':        (8.314462618, 'J/mol/K', 'gas constant'),
    '_faraday_primitive_sat':      (96485.33212, 'C/mol', 'Faraday'),
    # Cosmology
    '_h0_primitive_sat':           (67.4,       'km/s/Mpc', 'H_0'),
    '_t0_primitive_sat':           (13.787,     'Gyr', 'age of universe'),
    '_Omega_m_primitive_sat':      (0.315,      '', 'Omega_m'),
    '_Omega_b_h2_primitive_sat':   (0.02237,    '', 'Omega_b h^2'),
    '_Omega_DM_h2_primitive_sat':  (0.1200,     '', 'Omega_DM h^2'),
    '_n_s_primitive_sat':          (0.965,      '', 'n_s scalar tilt'),
    '_A_s_primitive_sat':          (2.1e-9,     '', 'A_s amplitude'),
    '_eta_primitive_sat':          (6.1e-10,    '', 'baryon-to-photon eta'),
    '_Y_p_primitive_sat':          (0.245,      '', 'primordial He fraction'),
    '_z_re_primitive_sat':         (7.7,        '', 'reionization redshift'),
    '_tau_reion_primitive_sat':    (0.054,      '', 'reionization optical depth'),
    '_w_z05_primitive_sat':        (-1.0,       '', 'w(z=0.5) dark-energy EOS'),
    '_f_NL_primitive_sat':         (0.0,        '', 'f_NL non-Gaussianity'),
    '_sigma_8_primitive_sat':      (0.811,      '', 'sigma_8'),
    '_T_CMB_primitive_sat':        (2.7255,     'K', 'CMB temperature'),
    '_T_neutrino_primitive_sat':   (1.9454,     'K', 'neutrino temperature'),
    '_BAO_rd_primitive_sat':       (147.18,     'Mpc', 'BAO sound horizon'),
    '_q0_decel_primitive_sat':     (-0.55,      '', 'deceleration parameter'),
    '_omega_k_curvature_primitive_sat': (0.0,   '', 'Omega_k'),
    '_h_dimensionless_primitive_sat': (0.674,   '', 'h dimensionless'),
    '_t_hubble_primitive_sat':     (14.5,       'Gyr', 'Hubble time'),
    '_sigma_v_cluster_primitive_sat': (1000.0,  'km/s', 'cluster vel dispersion'),
    '_growth_f_primitive_sat':     (0.53,       '', 'growth rate f'),
    # Multi-messenger
    '_eht_sgrA_shadow_primitive_sat':  (51.8,    'uas',   'EHT Sgr A* shadow'),
    '_gw150914_ringdown_primitive_sat':(251.0,   'Hz',    'GW150914 ringdown'),
    '_gw170817_inspiral_primitive_sat':(700.0,   'Hz',    'GW170817 inspiral'),
    '_jades_gs_z14_mass_primitive_sat':(4.0e8,   'M_sun', 'JADES-z14 stellar mass'),
    '_hudf_z_primitive_sat':       (10.0,        '',     'HUDF max z'),
    # Precision
    '_von_klitzing_primitive_sat': (25812.807, 'ohm',    'von Klitzing R_K'),
    '_josephson_primitive_sat':    (4.83597848e14, 'Hz/V', 'Josephson K_J'),
    # Planck units (anchor)
    '_planck_mass_primitive_sat':       (2.176434e-8, 'kg', 'Planck mass'),
    '_planck_length_primitive_sat':     (1.616255e-35, 'm',  'Planck length'),
    '_planck_time_primitive_sat':       (5.391247e-44, 's',  'Planck time'),
    '_planck_temperature_primitive_sat':(1.416784e32,  'K',  'Planck temperature'),
    # P-predictions (loose)
    '_p1_lkk_um_primitive_sat':    (55.0,       'um',  'P1 KK length'),
    '_p3_w0_primitive_sat':        (-1.0,       '',    'P3 w_0'),
    '_p4_dwdz2_primitive_sat':     (0.0,        '',    'P4 dw/dz^2'),
    # LENR (already audited)
    '_lenr_rossi_primitive_sat':       (10.0,   'COP', 'Rossi COP (10W)'),
    '_lenr_parkhomov_primitive_sat':   (200.0,  'W',   'Parkhomov 200W'),
    '_lenr_pons_fleischmann_primitive_sat':(10.0,'W', 'Pons-F 10W'),
    '_lenr_mizuno_primitive_sat':      (50.0,   'W',   'Mizuno 50W'),
    '_lenr_mckubre_primitive_sat':     (20.0,   'W',   'McKubre 20W'),
    '_lenr_stringham_primitive_sat':   (20.0,   'W',   'Stringham 20W'),
    '_lenr_brillouin_primitive_sat':   (50.0,   'W',   'Brillouin 50W'),
}

def order_decades(a, b):
    if a == 0 or b == 0:
        return float('inf')
    return abs(math.log10(abs(a)/abs(b)))

rows = []
for fname, (anchor, unit, desc) in ANCHORS.items():
    fn = getattr(u, fname, None)
    if fn is None:
        rows.append((fname, None, anchor, unit, desc, 'MISSING'))
        continue
    try:
        v = float(fn())
    except Exception as exc:
        rows.append((fname, None, anchor, unit, desc, f'ERROR:{exc}'))
        continue
    if anchor == 0:
        status = 'zero-anchor' if v == 0 else f'nonzero-vs-zero'
    else:
        d = order_decades(v, anchor)
        if d < 0.1:
            status = 'CLOSE'
        elif d < 1:
            status = f'~{d:.2f}dec'
        elif d < 5:
            status = f'WIDE({d:.1f}dec)'
        else:
            status = f'PLACEHOLDER({d:.0f}dec)'
    rows.append((fname, v, anchor, unit, desc, status))

# Print sorted by status: PLACEHOLDER first
def sort_key(r):
    s = r[5]
    if s.startswith('PLACEHOLDER'):
        return (0, -float(s[12:-4]) if s != 'PLACEHOLDER' else 0)
    if s.startswith('WIDE'):
        return (1, -float(s[5:-5]))
    if s.startswith('~'):
        return (2, -float(s[1:-3]))
    return (3, 0)

rows.sort(key=sort_key)

print(f"{'function':<42} {'value':>14} {'anchor':>14} {'unit':<10} {'status':<22} desc")
print('-' * 130)
for fname, v, anchor, unit, desc, status in rows:
    vstr = f'{v:.3e}' if isinstance(v,(int,float)) else 'N/A'
    astr = f'{anchor:.3e}' if isinstance(anchor,(int,float)) else str(anchor)
    print(f'{fname:<42} {vstr:>14} {astr:>14} {unit:<10} {status:<22} {desc}')

# Summary
n_placeholder = sum(1 for r in rows if r[5].startswith('PLACEHOLDER'))
n_wide = sum(1 for r in rows if r[5].startswith('WIDE'))
n_close = sum(1 for r in rows if r[5] in ('CLOSE','zero-anchor'))
print()
print(f'TOTAL: {len(rows)} | CLOSE: {n_close} | WIDE: {n_wide} | PLACEHOLDER: {n_placeholder}')
