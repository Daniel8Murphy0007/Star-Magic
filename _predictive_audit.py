import re
with open('uqff_pure_calculator.py', encoding='utf-8') as f:
    txt = f.read()
anchors_calls = len(re.findall(r'_sm_literal_anchor\(', txt))
primitives_calls = len(re.findall(r'_master_constant_primitive\(', txt))
canonical_sat = len(re.findall(r'_canonical_sat\b', txt))
codata = len(re.findall(r'CODATA', txt))
print('_sm_literal_anchor calls (anchor-only):', anchors_calls)
print('_master_constant_primitive calls (predictive):', primitives_calls)
print('_canonical_sat references (Layer 45 saturations):', canonical_sat)
print('CODATA mentions:', codata)
anch = re.findall(r"return _sm_literal_anchor\([\"']([^\"']+)", txt)
print('constants routed to anchor:', sorted(set(anch)))
prim = re.findall(r"return _master_constant_primitive\([\"']([^\"']+)", txt)
print('constants routed to primitive:', sorted(set(prim)))
# Anchor table keys
keys = re.findall(r"_SM_LITERAL_ANCHOR_SAT\s*=\s*\{([\s\S]*?)\n\}", txt)
if keys:
    found_keys = re.findall(r"[\"']([a-zA-Z_][\w]*)[\"']\s*:", keys[0])
    print('anchor table entries:', len(found_keys))
    print(sorted(set(found_keys)))
