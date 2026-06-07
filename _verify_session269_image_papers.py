"""Cold verification of Session 269 image-batch wiring (26 PDFs: PAPER_1186..PAPER_1200).
10 closure-bearing + 15 framework-only + 1 SCm-dupe.
All 50 closures use form=None (paper bodies don't transcribe forms or use primitives
outside the locked vacuum-ledger). Paper-stated predicted values exposed verbatim
-- same honest pattern as S626 Faraday (Session 268)."""
import sys
sys.path.insert(0, '.')
from uqff_pure_calculator import (
    _PA_S269_CLOSURE_REGISTRY,
    _l96_uqff_closure_S269_lookup,
    # 10 closure probes
    _l96_uqff_PAPER1186_high_z_quasar_accretion_probe,
    _l96_uqff_PAPER1187_cooling_flow_mass_accretion_probe,
    _l96_uqff_PAPER1188_magnetar_thermal_conductivity_probe,
    _l96_uqff_PAPER1189_photoevaporation_hz_orion_probe,
    _l96_uqff_PAPER1190_alma_molecular_gas_probe,
    _l96_uqff_PAPER1191_gw190425_mass_gap_probe,
    _l96_uqff_PAPER1192_snr_shock_velocity_probe,
    _l96_uqff_PAPER1193_pvsnp_conjecture_probe,
    _l96_uqff_PAPER1196_plasma_fusion_probe,
    _l96_uqff_PAPER1199_information_math_constants_probe,
    # 15 framework probes
    _l96_uqff_PAPER1186_standard_model_unified_probe,
    _l96_uqff_PAPER1187_cosmological_tensions_unified_probe,
    _l96_uqff_PAPER1188_number_theory_frontier_probe,
    _l96_uqff_PAPER1189_chemistry_atomic_unified_probe,
    _l96_uqff_PAPER1190_mathematical_constants_unified_probe,
    _l96_uqff_PAPER1191_cosmology_deepset_unified_probe,
    _l96_uqff_PAPER1192_standard_model_deepcuts_unified_probe,
    _l96_uqff_PAPER1193_astrophysics_unified_probe,
    _l96_uqff_PAPER1194_tde_rate_mass_probe,
    _l96_uqff_PAPER1194_condensed_matter_unified_probe,
    _l96_uqff_PAPER1195_biology_unified_probe,
    _l96_uqff_PAPER1197_universal_buoyancy_solver_probe,
    _l96_uqff_PAPER1197_geophysics_atmospheric_unified_probe,
    _l96_uqff_PAPER1198_rhovacscm_derivation_probe,
    _l96_uqff_PAPER1198_particle_physics_unified_probe,
    _l96_uqff_session269_image_paper_manifest,
    calculate_analytic_closures,
)


def report(name, probe):
    derived = probe.get("derived", {})
    exact   = probe.get("exact_count", 0)
    print(f"\n{'='*78}\n{name}\n{'='*78}")
    print(f"  closures   = {len(derived)}")
    print(f"  EXACT      = {exact}")
    print(f"  all <1%PRED= {probe.get('all_within_1pct_predicted')}")
    print(f"  all <1%OBS = {probe.get('all_within_1pct_observed')}")


probes = [
    ("PAPER_1186 High-z Quasar Accretion (S673)",       _l96_uqff_PAPER1186_high_z_quasar_accretion_probe()),
    ("PAPER_1187 Cooling Flow (S674-S677)",             _l96_uqff_PAPER1187_cooling_flow_mass_accretion_probe()),
    ("PAPER_1188 Magnetar Thermal Cond. (S678-S681)",   _l96_uqff_PAPER1188_magnetar_thermal_conductivity_probe()),
    ("PAPER_1189 Photoevap HZ Orion (S682-S686)",       _l96_uqff_PAPER1189_photoevaporation_hz_orion_probe()),
    ("PAPER_1190 ALMA Molecular Gas (S687-S690)",       _l96_uqff_PAPER1190_alma_molecular_gas_probe()),
    ("PAPER_1191 GW190425 Mass Gap (S691-S694)",        _l96_uqff_PAPER1191_gw190425_mass_gap_probe()),
    ("PAPER_1192 SNR Shock Velocity (S695-S698)",       _l96_uqff_PAPER1192_snr_shock_velocity_probe()),
    ("PAPER_1193 P vs NP Conjecture (S699-S702)",       _l96_uqff_PAPER1193_pvsnp_conjecture_probe()),
    ("PAPER_1196 Plasma Fusion (S703-S712)",            _l96_uqff_PAPER1196_plasma_fusion_probe()),
    ("PAPER_1199 Information Math (S713-S722)",         _l96_uqff_PAPER1199_information_math_constants_probe()),
]
total = 0
for name, p in probes:
    report(name, p)
    total += p.get('closures', 0)

print(f"\n{'='*78}\nTOTAL closures wired: {total}  registry size: {len(_PA_S269_CLOSURE_REGISTRY)}")

# Framework-only probes
print(f"\n{'='*78}\nFRAMEWORK-ONLY PAPERS (15)\n{'='*78}")
for name, fn in [
    ("PAPER_1186 Standard Model Unified",         _l96_uqff_PAPER1186_standard_model_unified_probe),
    ("PAPER_1187 Cosmological Tensions",          _l96_uqff_PAPER1187_cosmological_tensions_unified_probe),
    ("PAPER_1188 Number Theory Frontier",         _l96_uqff_PAPER1188_number_theory_frontier_probe),
    ("PAPER_1189 Chemistry/Atomic Unified",       _l96_uqff_PAPER1189_chemistry_atomic_unified_probe),
    ("PAPER_1190 Mathematical Constants",         _l96_uqff_PAPER1190_mathematical_constants_unified_probe),
    ("PAPER_1191 Cosmology Deepset",              _l96_uqff_PAPER1191_cosmology_deepset_unified_probe),
    ("PAPER_1192 SM Deepcuts",                    _l96_uqff_PAPER1192_standard_model_deepcuts_unified_probe),
    ("PAPER_1193 Astrophysics Unified",           _l96_uqff_PAPER1193_astrophysics_unified_probe),
    ("PAPER_1194 TDE Rate Mass",                  _l96_uqff_PAPER1194_tde_rate_mass_probe),
    ("PAPER_1194 Condensed Matter",               _l96_uqff_PAPER1194_condensed_matter_unified_probe),
    ("PAPER_1195 Biology",                        _l96_uqff_PAPER1195_biology_unified_probe),
    ("PAPER_1197 Universal Buoyancy Solver",      _l96_uqff_PAPER1197_universal_buoyancy_solver_probe),
    ("PAPER_1197 Geophysics Atmospheric",         _l96_uqff_PAPER1197_geophysics_atmospheric_unified_probe),
    ("PAPER_1198 RhoVacSCm Derivation",           _l96_uqff_PAPER1198_rhovacscm_derivation_probe),
    ("PAPER_1198 Particle Physics",               _l96_uqff_PAPER1198_particle_physics_unified_probe),
]:
    p = fn()
    print(f"  {name:<48} type='{p['type']}'")

# Manifest
print(f"\n{'='*78}\nSession 269 manifest\n{'='*78}")
mf = _l96_uqff_session269_image_paper_manifest()
print(f"  wired with closures              = {len(mf['wired_with_closures_this_session'])}")
print(f"  framework-only                   = {len(mf['framework_only_metadata_this_session'])}")
print(f"  SCm dupes (existing papers)      = {len(mf['scm_dupes_of_existing_papers'])}")
print(f"  total new closures wired         = {mf['total_new_closures_wired']}")
print(f"  total new probe functions        = {mf['total_new_probe_functions']}")

# Universal lookup spot-check
print(f"\n{'='*78}\nUniversal S### lookup spot-check\n{'='*78}")
for sid in ["S673", "S677", "S686", "S694", "S702", "S707", "S722"]:
    info = _l96_uqff_closure_S269_lookup(sid)
    pred = info.get('paper_predicted')
    obs  = info.get('observed_anchor')
    pred_s = f"{pred:.4g}" if isinstance(pred, (int, float)) else str(pred)
    obs_s  = f"{obs:.4g}"  if isinstance(obs,  (int, float)) else str(obs)
    print(f"  {sid:<6} {info['paper_tag']:<10} {info['label']:<42} pred={pred_s:<10} obs={obs_s}")

# Dispatcher routing
print(f"\n{'='*78}\nDISPATCHER ROUTING CHECK\n{'='*78}")
for key in ["paper1186_high_z_quasar_accretion_probe",
            "paper1191_gw190425_mass_gap_probe",
            "paper1193_pvsnp_conjecture_probe",
            "paper1196_plasma_fusion_probe",
            "paper1199_information_math_constants_probe",
            "paper1198_rhovacscm_derivation_probe",
            "paper1195_biology_unified_probe",
            "image_session269_manifest"]:
    try:
        out = calculate_analytic_closures({"calc": key})
        v = out.get("value")
        v_show = f"<dict, {len(v)} keys>" if isinstance(v, dict) else f"{v}"
        print(f"  {key:<50} -> OK")
    except Exception as e:
        print(f"  {key:<50} -> FAIL  {type(e).__name__}: {e}")

# Universal lookup via dispatcher
try:
    out = calculate_analytic_closures({"calc": "closure_s269_lookup", "s_id": "S696"})
    print(f"  closure_s269_lookup(S696) -> OK")
except Exception as e:
    print(f"  closure_s269_lookup(S696) -> FAIL {type(e).__name__}: {e}")
