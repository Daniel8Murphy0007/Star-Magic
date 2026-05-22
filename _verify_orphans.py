"""Read master_closures.csv and print which of the 30 EXTRA derivation files parsed OK
vs PARSE_FAIL. Also count sidecar files for forensic completeness."""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CSV = ROOT / "master_closures.csv"
AUDIT = ROOT / "_audit_outputs"

EXTRAS = {
    "_chain_trace_26layer.py","_chain_trace_C.py","_chain_trace_C_particles.py",
    "_chain_trace_fix348.py","_chain_trace_fix56_7_910.py","_chain_trace_np_split.py",
    "_chain_trace_SSq.py","_constant_derivation_attempt.py","_constant_derivation_v2.py",
    "_constant_derivation_v3.py","_K_Mex_REAL_derivation.py","_lagrangian_rederivation_outline.py",
    "_PAPER_1065_1066_variational_audit.py","_PAPER_1183_first_principles_derivation.py",
    "_variational_sustainability_solution.py","bsm_bounds_derivation.py",
    "buoyancy_lagrangian_eom.py","et_full_lagrangian.py","first_principles_derivation.py",
    "lagrangian_re_runner.py","qcalcgeom_core_derivation.py","thorne_morris_exotic_derivation.py",
    "uqff_lagrangian_derivation.py","UQFF_UNIFIED_CLOSURE_DERIVATIONS.py",
    "variational_reversal_condition.py","vds_dvp_bsh_symbolic_proofs.py",
    "_six_anchor_closures.py","_matter_density_closures.py",
    "_cosmological_closures.py","_lambda_closure_v1.py",
}

rows = {}
with CSV.open("r", encoding="utf-8", newline="") as f:
    for r in csv.DictReader(f):
        if r["script"] in EXTRAS:
            rows[r["script"]] = r

ok, fails = [], []
for fname in sorted(EXTRAS):
    if fname not in rows:
        fails.append((fname, "NOT IN CSV", "", ""))
        continue
    r = rows[fname]
    status = r["status"]
    if status in ("OK", "OK_JSON"):
        ok.append((fname, status, r["predicted"], r["observed"], r["error_pct"]))
    else:
        fails.append((fname, status, r["label"], r["raw_output"][:80]))

print("=" * 80)
print(f"EXTRA DERIVATION FILES: {len(EXTRAS)} total")
print(f"  OK / OK_JSON: {len(ok)}")
print(f"  PARSE_FAIL:   {len(fails)}")
print("=" * 80)
print()
print(f"--- {len(ok)} ORPHANS THAT WORK (parsed cleanly) ---")
for fname, status, pred, obs, err in ok:
    print(f"  [{status:8s}] {fname:50s}  err={err}%")
print()
print(f"--- {len(fails)} ORPHANS THAT DO NOT PARSE ---")
for tup in fails:
    fname = tup[0]; status = tup[1]
    print(f"  [{status:10s}] {fname}")
print()
print("=" * 80)
print(f"Sidecar capture: _audit_outputs/ exists = {AUDIT.exists()}")
if AUDIT.exists():
    files = list(AUDIT.glob("*.txt"))
    print(f"  total stdout sidecars captured: {len(files)}")
    extras_captured = sum(1 for f in files if f.name.endswith(".txt") and (f.stem + ".py") in EXTRAS)
    print(f"  EXTRA derivation sidecars present: {extras_captured} / {len(EXTRAS)}")
print("=" * 80)
