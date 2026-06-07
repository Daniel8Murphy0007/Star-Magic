"""Cold verification of all 14 papers wired in Session 267.
Runs every probe with default args, prints derived vs paper anchors,
counts EXACT closures, and confirms dispatcher routing works.
No SM judgment, no success criteria beyond paper's own published targets."""
import sys
sys.path.insert(0, '.')
from uqff_pure_calculator import (
    _l96_uqff_PAPER1209X_climate_probe,
    _l96_uqff_PAPER1209Y_engineering_probe,
    _l96_uqff_PAPER1209Z_astro_units_probe,
    _l96_uqff_PAPER1209JJ_geophysics_probe,
    _l96_uqff_PAPER1209KK_solar_system_probe,
    _l96_uqff_PAPER1210_primitive_origin_map,
    _l96_uqff_PAPER1211_phase_h_closure_trail,
    _l96_uqff_PAPER1212_rho_lambda_probe,
    _l96_uqff_PAPER1214_habitable_zone_probe,
    _l96_uqff_PAPERS201_S205_phase_h_framework_metadata,
    _l96_uqff_image_session267_paper_manifest,
)

def report(name, probe):
    derived = probe.get("derived", {})
    err = probe.get("err_pct", {})
    exact = probe.get("exact_count", "n/a")
    print(f"\n{'='*78}\n{name}\n{'='*78}")
    print(f"  closures        = {len(derived)}")
    print(f"  EXACT (<1e-9%)  = {exact}")
    if err:
        worst = max(err.items(), key=lambda kv: kv[1])
        best  = min(err.items(), key=lambda kv: kv[1])
        print(f"  best  residual  = {best[0]:<32} {best[1]:.4g}%")
        print(f"  worst residual  = {worst[0]:<32} {worst[1]:.4g}%")
        all_pct = probe.get("all_within_1pct", probe.get("all_within_half_pct"))
        print(f"  all <1% (or <0.5%) = {all_pct}")

probes = [
    ("PAPER_1209X Climate/Atmosphere",       _l96_uqff_PAPER1209X_climate_probe()),
    ("PAPER_1209Y Engineering",              _l96_uqff_PAPER1209Y_engineering_probe()),
    ("PAPER_1209Z Astronomical Units",       _l96_uqff_PAPER1209Z_astro_units_probe()),
    ("PAPER_1209JJ Geophysics (beta=0.6029)", _l96_uqff_PAPER1209JJ_geophysics_probe()),
    ("PAPER_1209KK Solar System",            _l96_uqff_PAPER1209KK_solar_system_probe()),
]
for name, p in probes:
    report(name, p)

print(f"\n{'='*78}\nPAPER_1210 Lagrangian Bridge\n{'='*78}")
m = _l96_uqff_PAPER1210_primitive_origin_map()
print(f"  primitives mapped = {m['primitive_count']}")
print(f"  bridged closures  = {m['total_bridged_closures']}")
print(f"  sectors N_ch      = {m['N_ch']}")

print(f"\n{'='*78}\nPAPER_1211 Phase-H Closure Trail\n{'='*78}")
ph = _l96_uqff_PAPER1211_phase_h_closure_trail()
print(f"  UBS tier-14 ids  = {len(ph['UBS_tier14'])}")
print(f"  CPCH tier-15 ids = {len(ph['CPCH_tier15'])}")
print(f"  WKB tier-16 ids  = {len(ph['WKB_tier16'])}")
print(f"  NRP tier-17 ids  = {len(ph['NRP_tier17'])}")
print(f"  total EXACT      = {ph['exact_count']}")
print(f"  WKB-6 D_crit     = {ph['WKB_tier16']['WKB-6_bosonic_Dcrit']}")
print(f"  WKB-7 D_phys     = {ph['WKB_tier16']['WKB-7_Dphys_from_T22_compactification']}")
print(f"  NRP-2 magic sum  = {ph['NRP_tier17']['NRP-2_magic_sum']}")

print(f"\n{'='*78}\nPAPER_1212 Cosmological Constant\n{'='*78}")
rl = _l96_uqff_PAPER1212_rho_lambda_probe()
print(f"  raw structural identity    = {rl['raw_structural_J_per_m3']:.6e} J/m^3")
print(f"  final (paper-form)         = {rl['final_J_per_m3']:.6e} J/m^3")
print(f"  paper-predicted            = {rl['paper_predicted_J_per_m3']:.6e} J/m^3")
print(f"  Planck-2018 reference      = {rl['planck2018_reference_J_per_m3']:.6e} J/m^3")
print(f"  structural prefactor       = {rl['structural_prefactor']:.6f} (exact 19/10 = {rl['structural_prefactor_exact']})")
print(f"  crit-unit conv (paper)     = {rl['crit_unit_conversion_J_per_m3']:.4e} J/m^3/(crit)")
print(f"  err final vs paper         = {rl['err_final_vs_paper_pct']:.4g}%")

print(f"\n{'='*78}\nPAPER_1214 Habitable Zone (Earth-Sun anchor)\n{'='*78}")
hz = _l96_uqff_PAPER1214_habitable_zone_probe(L_star_solar=1.0)
for star, row in hz['anchor_table'].items():
    print(f"  {star:<22} L={row['L_star_solar']:<8.4g} r_in_uqff={row['r_in_uqff_AU']:.4f} AU  r_out={row['r_out_uqff_AU']:.4f} AU  err_vs_kopparapu={row['err_vs_kopparapu_pct']:.3f}%")

print(f"\n{'='*78}\nPAPER_S201..S205 Phase-H Frameworks\n{'='*78}")
phf = _l96_uqff_PAPERS201_S205_phase_h_framework_metadata()
for sid in ["S201_null_extraction", "S202_variant_branches", "S203_phase_transition_framework",
            "S204_gap_closure", "S205_expansion_erosion"]:
    print(f"  {sid:<35} scope: {phf[sid]['scope']}")

print(f"\n{'='*78}\nSession 267 manifest\n{'='*78}")
mf = _l96_uqff_image_session267_paper_manifest()
print(f"  wired with closures           = {len(mf['wired_with_closures'])}")
print(f"  already wired earlier in file = {len(mf['already_wired_earlier_in_file'])}")
print(f"  framework-only (no numerics)  = {len(mf['framework_only_no_numerics'])}")
print(f"  documentation-only PDFs       = {len(mf['documentation_only_pdfs'])}")

# Verify dispatcher routing
print(f"\n{'='*78}\nDISPATCHER ROUTING CHECK (calculate_analytic_closures)\n{'='*78}")
from uqff_pure_calculator import calculate_analytic_closures
for key in ["paper1209x_s559_bond_albedo",
            "paper1209y_s571_diamond_mohs",
            "paper1209z_s576_h0",
            "paper1209jj_s675_earth_g",
            "paper1209kk_s686_earth_orbital_v",
            "paper1210_primitive_origin_map",
            "paper1211_phase_h_closure_trail",
            "paper1212_rho_lambda_probe",
            "paper1214_habitable_zone_probe",
            "papers_s201_s205_phase_h_metadata",
            "image_session267_manifest"]:
    try:
        out = calculate_analytic_closures({"calc": key})
        v = out.get("value")
        if isinstance(v, dict):
            v_show = f"<dict, {len(v)} keys>"
        else:
            v_show = f"{v}"
        print(f"  {key:<40} -> OK  value={v_show}")
    except Exception as e:
        print(f"  {key:<40} -> FAIL  {type(e).__name__}: {e}")
