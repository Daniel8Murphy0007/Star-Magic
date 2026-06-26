#!/usr/bin/env python3
"""scan_paradox_dispatch.py — programmatic walk of PARADOX_TO_CLOSURE in uqff_pure_calculator.py.

Reads Daniel's calculator (READ-ONLY) and extracts every paradox dispatch key + closure function
reference. Reports counts, categories, sample formulas.

Writes a JSON index to MY sandbox.
"""
import re
import json
from pathlib import Path

CALC = Path("/sessions/vibrant-keen-bohr/mnt/Star-Magic/uqff_pure_calculator.py")
OUT  = Path("/sessions/vibrant-keen-bohr/mnt/outputs/uqff_recompute/paradox_dispatch_index.json")

print(f"Reading: {CALC}  (READ-ONLY)")
src = CALC.read_text(encoding="utf-8", errors="replace")
print(f"  size: {len(src):,} bytes, {src.count(chr(10)):,} lines")

# 1. find PARADOX_TO_CLOSURE dict definition (it's a large literal)
m = re.search(r"PARADOX_TO_CLOSURE\s*=\s*\{", src)
if not m:
    print("ERROR: PARADOX_TO_CLOSURE not found")
    raise SystemExit(1)
start = m.end()
# find matching close brace by depth tracking
depth = 1
i = start
while i < len(src) and depth > 0:
    if src[i] == '{': depth += 1
    elif src[i] == '}': depth -= 1
    i += 1
dict_body = src[start:i-1]
print(f"PARADOX_TO_CLOSURE dict body: {len(dict_body):,} chars, lines {src[:start].count(chr(10))+1} → {src[:i].count(chr(10))+1}")

# 2. extract keys and closure-fn names
# Pattern: "key": func_name  OR  "key": lambda ...
entries = []
for em in re.finditer(r'["\']([a-zA-Z0-9_\-\.]+)["\']\s*:\s*([a-zA-Z_][a-zA-Z0-9_]*)', dict_body):
    key = em.group(1)
    fn = em.group(2)
    entries.append({"key": key, "fn": fn})
# lambda entries (different pattern)
for em in re.finditer(r'["\']([a-zA-Z0-9_\-\.]+)["\']\s*:\s*lambda', dict_body):
    entries.append({"key": em.group(1), "fn": "<lambda>"})

print(f"Extracted {len(entries):,} entries from PARADOX_TO_CLOSURE")
print()

# 3. categorize keys
cats = {}
def category(key):
    k = key.lower()
    for tag, words in [
        ("quantum_mechanics", ["epr","ghz","hardy","kochen","bell","wigner","quantum","entangle","decoherence","copenhagen","wavefn"]),
        ("relativity_blackhole", ["firewall","amps","page","hawking","schwarzschild","kerr","horizon","singularity","information_paradox","blackhole","bh_"]),
        ("cosmology", ["hubble","sigma_8","dark_flow","dark_energy","lambda","cdm","cmb","reionization","inflation","cosmic","planck","desi"]),
        ("qft", ["strong_cp","axion","gauge_coupling","vacuum_catastrophe","ckm","pmns","running","beta_function","qcd","qed"]),
        ("sm_identities", ["m_w","m_z","m_t","m_h","m_b","m_c","m_tau","m_mu","m_e","m_p","fermion","neutrino","mass_ratio","mixing"]),
        ("nuclear_atomic", ["magic","be_per_a","alpha_binding","deuteron","rydberg","proton_radius","binding_energy","shell_model"]),
        ("lenr_fusion", ["holmlid","parkhomov","pons","mizuno","rossi","star_magic","reactor","cop","ker","transmutation","ecat"]),
        ("astro", ["psr","crab","eta_car","supernova","imf","galaxy","sgr_a","ngc_","gw1","ligo","virgo","kepler"]),
        ("statistical_info", ["monty","sleeping","simpson","bertrand","entropy","shannon","boltzmann"]),
        ("math", ["hilbert","banach","peano","galois","riemann","poincare","hodge","bsd","pvsnp","p_vs_np","navier","yang_mills"]),
        ("engineering", ["surface_code","qubit","threshold","iter","tokamak","fusion","reactor"]),
        ("biology_other", ["genetic","amino","codon","hayflick","mosquito","dna","aging"]),
    ]:
        for w in words:
            if w in k:
                return tag
    return "other"

for e in entries:
    c = category(e["key"])
    cats.setdefault(c, []).append(e["key"])

print(f"{'category':<25} {'count':>6}  sample keys")
print("-" * 100)
for c in sorted(cats, key=lambda k: -len(cats[k])):
    s = ", ".join(cats[c][:3])
    print(f"{c:<25} {len(cats[c]):>6}  {s[:80]}")
print()

# 4. extract _MILLENNIUM_DERIVE (separate small dict)
mm = re.search(r"_MILLENNIUM_DERIVE\s*=\s*\{", src)
mill_entries = []
if mm:
    s2 = mm.end()
    depth = 1; j = s2
    while j < len(src) and depth > 0:
        if src[j] == '{': depth += 1
        elif src[j] == '}': depth -= 1
        j += 1
    body = src[s2:j-1]
    for em in re.finditer(r'["\']([a-zA-Z0-9_\-]+)["\']\s*:\s*([a-zA-Z_][a-zA-Z0-9_]*)', body):
        mill_entries.append({"key": em.group(1), "fn": em.group(2)})
    print(f"_MILLENNIUM_DERIVE entries: {len(mill_entries)}")
    for me in mill_entries:
        print(f"  {me['key']:<35} -> {me['fn']}")
    print()

# 5. find every _l96_uqff_axiom_*_closure function definition
l96 = re.findall(r"^def\s+(_l96_uqff_axiom_[a-zA-Z0-9_]+_closure)\s*\(", src, flags=re.MULTILINE)
print(f"_l96_uqff_axiom_*_closure functions defined: {len(l96)}")

# 6. find every _millennium_*_derive function definition
mill_fns = re.findall(r"^def\s+(_millennium_[a-zA-Z0-9_]+_derive)\s*\(", src, flags=re.MULTILINE)
print(f"_millennium_*_derive functions defined: {len(mill_fns)}")
for f in mill_fns:
    print(f"  {f}")
print()

# 7. count all public calculate_* surfaces
pubs = re.findall(r"^def\s+(calculate_[a-zA-Z0-9_]+)\s*\(", src, flags=re.MULTILINE)
print(f"public calculate_* surfaces: {len(pubs)}")
for p in pubs:
    print(f"  {p}")
print()

# 8. write JSON index
out_data = {
    "paradox_dispatch_entries": entries,
    "paradox_categories": {c: len(v) for c, v in cats.items()},
    "millennium_derive_entries": mill_entries,
    "l96_closure_functions": l96,
    "millennium_derive_functions": mill_fns,
    "public_calculate_surfaces": pubs,
    "calculator_source_path": str(CALC),
    "calculator_source_size_bytes": len(src),
    "calculator_source_lines": src.count('\n'),
}
OUT.write_text(json.dumps(out_data, indent=2))
print(f"Wrote index to: {OUT}")
print(f"Total entries indexed: PARADOX_TO_CLOSURE={len(entries)} + _MILLENNIUM_DERIVE={len(mill_entries)} + _l96_uqff_axiom_*={len(l96)} + _millennium_*_derive={len(mill_fns)} + calculate_*={len(pubs)}")
