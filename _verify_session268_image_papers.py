"""Cold verification of all 18 closure-bearing papers + 4 framework papers
wired in Session 268 image batch (round 2).
Runs every probe, reports closures count, EXACT count, worst residual vs
paper-stated PREDICTED column, worst vs OBSERVED anchor.
No SM judgment. Paper anchors are user-published targets, not success criteria."""
import sys
sys.path.insert(0, '.')
from uqff_pure_calculator import (
    _PA_S268_CLOSURE_REGISTRY,
    _l96_uqff_closure_S268_lookup,
    _l96_uqff_PAPER1200_GR_precision_probe,
    _l96_uqff_PAPER1201_materials_photonics_probe,
    _l96_uqff_PAPER1202_chemistry_spectroscopy_probe,
    _l96_uqff_PAPER1203_nuclear_physics_probe,
    _l96_uqff_PAPER1204_fluid_dynamics_probe,
    _l96_uqff_PAPER1205_geometry_topology_probe,
    _l96_uqff_PAPER1206_solar_system_probe,
    _l96_uqff_PAPER1207_biology_allometry_probe,
    _l96_uqff_PAPER1208_transcendentals_probe,
    _l96_uqff_PAPER1209_particle_physics_probe,
    _l96_uqff_PAPER1209AA_chemistry_probe,
    _l96_uqff_PAPER1209BB_biology_probe,
    _l96_uqff_PAPER1209CC_geophysics_probe,
    _l96_uqff_PAPER1209DD_electromagnetism_probe,
    _l96_uqff_PAPER1209EE_quantum_thermo_probe,
    _l96_uqff_PAPER1209FF_math_constants_probe,
    _l96_uqff_PAPER1209GG_cosmological_constants_probe,
    _l96_uqff_PAPER1209II_nuclear_binding_energies_probe,
    _l96_uqff_PAPER1200_FUBi_FUBii_derived_G_probe,
    _l96_uqff_PAPER1201_26D_polynomial_origami_probe,
    _l96_uqff_PAPER1202_quantum_chain_E_n_summation_probe,
    _l96_uqff_PAPER1203_canonical_v15_simultaneous_solver_probe,
    _l96_uqff_session268_image_paper_manifest,
    calculate_analytic_closures,
)


def report(name, probe):
    derived  = probe.get("derived", {})
    err_pred = probe.get("err_vs_predicted_pct", {})
    err_obs  = probe.get("err_vs_observed_pct", {})
    exact    = probe.get("exact_count", 0)
    print(f"\n{'='*78}\n{name}\n{'='*78}")
    print(f"  closures        = {len(derived)}")
    print(f"  EXACT formulas  = {exact}  (derived matches paper-stated to <1e-9%)")
    if err_pred:
        worst_pred = max(err_pred.items(), key=lambda kv: kv[1])
        best_pred  = min(err_pred.items(), key=lambda kv: kv[1])
        print(f"  worst vs PRED   = {worst_pred[0]:<45} {worst_pred[1]:.4g}%")
        print(f"  best  vs PRED   = {best_pred[0]:<45} {best_pred[1]:.4g}%")
    if err_obs:
        worst_obs = max(err_obs.items(), key=lambda kv: kv[1])
        best_obs  = min(err_obs.items(), key=lambda kv: kv[1])
        print(f"  worst vs OBS    = {worst_obs[0]:<45} {worst_obs[1]:.4g}%")
        print(f"  best  vs OBS    = {best_obs[0]:<45} {best_obs[1]:.4g}%")
    print(f"  all <1% vs PRED = {probe.get('all_within_1pct_predicted')}")
    print(f"  all <1% vs OBS  = {probe.get('all_within_1pct_observed')}")


probes = [
    ("PAPER_1200 GR Precision (S453-S462)",        _l96_uqff_PAPER1200_GR_precision_probe()),
    ("PAPER_1201 Materials/Photonics (S463-S472)", _l96_uqff_PAPER1201_materials_photonics_probe()),
    ("PAPER_1202 Chemistry/Spectroscopy (S473-S482)", _l96_uqff_PAPER1202_chemistry_spectroscopy_probe()),
    ("PAPER_1203 Nuclear Physics (S483-S492)",     _l96_uqff_PAPER1203_nuclear_physics_probe()),
    ("PAPER_1204 Fluid Dynamics (S493-S502)",      _l96_uqff_PAPER1204_fluid_dynamics_probe()),
    ("PAPER_1205 Geometry/Topology (S503-S512)",   _l96_uqff_PAPER1205_geometry_topology_probe()),
    ("PAPER_1206 Solar System (S513-S522)",        _l96_uqff_PAPER1206_solar_system_probe()),
    ("PAPER_1207 Biology/Allometry (S523-S532)",   _l96_uqff_PAPER1207_biology_allometry_probe()),
    ("PAPER_1208 Transcendentals (S533-S542)",     _l96_uqff_PAPER1208_transcendentals_probe()),
    ("PAPER_1209 Particle Physics (S543-S552)",    _l96_uqff_PAPER1209_particle_physics_probe()),
    ("PAPER_1209AA Chemistry (S583-S592)",         _l96_uqff_PAPER1209AA_chemistry_probe()),
    ("PAPER_1209BB Biology (S593-S602)",           _l96_uqff_PAPER1209BB_biology_probe()),
    ("PAPER_1209CC Geophysics (S603-S612)",        _l96_uqff_PAPER1209CC_geophysics_probe()),
    ("PAPER_1209DD Electromagnetism (S613-S622)",  _l96_uqff_PAPER1209DD_electromagnetism_probe()),
    ("PAPER_1209EE Quantum/Thermo (S623-S632)",    _l96_uqff_PAPER1209EE_quantum_thermo_probe()),
    ("PAPER_1209FF Math Constants (S633-S642)",    _l96_uqff_PAPER1209FF_math_constants_probe()),
    ("PAPER_1209GG Cosmological (S643-S652)",      _l96_uqff_PAPER1209GG_cosmological_constants_probe()),
    ("PAPER_1209II Nuclear Binding (S663-S672)",   _l96_uqff_PAPER1209II_nuclear_binding_energies_probe()),
]
total_closures = 0
total_exact = 0
for name, p in probes:
    report(name, p)
    total_closures += p.get('closures', 0)
    total_exact    += p.get('exact_count', 0)

print(f"\n{'='*78}\nTOTAL: {total_closures} closures wired, {total_exact} formula-EXACT vs paper PREDICTED")
print(f"Registry size: {len(_PA_S268_CLOSURE_REGISTRY)} (should equal {total_closures})")

# Framework-only papers
print(f"\n{'='*78}\nFRAMEWORK-ONLY PAPERS (no S### closures)\n{'='*78}")
for name, fn in [
    ("PAPER_1200 FUBi/FUBii Derived G",                 _l96_uqff_PAPER1200_FUBi_FUBii_derived_G_probe),
    ("PAPER_1201 26D Polynomial Origami Axiom",         _l96_uqff_PAPER1201_26D_polynomial_origami_probe),
    ("PAPER_1202 Quantum Chain E_n Summation 633333",   _l96_uqff_PAPER1202_quantum_chain_E_n_summation_probe),
    ("PAPER_1203 Canonical v1.5 Simultaneous Solver",   _l96_uqff_PAPER1203_canonical_v15_simultaneous_solver_probe),
]:
    p = fn()
    print(f"  {name:<55} type='{p['type']}'")

# Manifest
print(f"\n{'='*78}\nSession 268 manifest\n{'='*78}")
mf = _l96_uqff_session268_image_paper_manifest()
print(f"  wired with closures this session    = {len(mf['wired_with_closures_this_session'])}")
print(f"  framework-only this session         = {len(mf['framework_only_metadata_this_session'])}")
print(f"  already wired earlier in file       = {len(mf['already_wired_earlier_in_file'])}")
print(f"  SCm LENR PDFs (existing papers)     = {len(mf['scm_lenr_pdfs_renumbered_from_existing_papers'])}")
print(f"  total new closures wired            = {mf['total_new_closures_wired']}")

# Universal lookup spot-check
print(f"\n{'='*78}\nUniversal S### lookup spot-check\n{'='*78}")
for sid in ["S453","S477","S499","S524","S538","S543","S593","S651","S669"]:
    info = _l96_uqff_closure_S268_lookup(sid)
    print(f"  {sid:<6} {info['label']:<42} derived={info['derived']:<14.6g} "
          f"pred={info['paper_predicted']:<10.6g} err_obs={info['err_vs_observed_pct']:.3g}%")

# Dispatcher routing
print(f"\n{'='*78}\nDISPATCHER ROUTING CHECK (calculate_analytic_closures)\n{'='*78}")
for key in ["paper1200_gr_precision_probe",
            "paper1208_transcendentals_probe",
            "paper1209bb_biology_probe",
            "paper1209gg_cosmological_constants_probe",
            "paper1209ii_nuclear_binding_energies_probe",
            "paper1200_fubi_fubii_derived_g_probe",
            "paper1201_26d_polynomial_origami_probe",
            "image_session268_manifest"]:
    try:
        out = calculate_analytic_closures({"calc": key})
        v = out.get("value")
        v_show = f"<dict, {len(v)} keys>" if isinstance(v, dict) else f"{v}"
        print(f"  {key:<55} -> OK  value={v_show}")
    except Exception as e:
        print(f"  {key:<55} -> FAIL  {type(e).__name__}: {e}")

# Universal lookup via dispatcher
try:
    out = calculate_analytic_closures({"calc": "closure_s268_lookup", "s_id": "S646"})
    print(f"  closure_s268_lookup(S646) -> {out.get('value', '<see internals>')}")
except Exception as e:
    print(f"  closure_s268_lookup(S646) -> FAIL {type(e).__name__}: {e}")
