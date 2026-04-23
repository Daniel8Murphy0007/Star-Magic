import json

with open('_cp1_failures.json') as f:
    fails = json.load(f)
print(f'Total CP1 failures: {len(fails)}')
print()

no_compute = [f for f in fails if 'no compute' in f['error']]
attr_dict   = [f for f in fails if 'dict' in f['error'] and 'has no attribute' in f['error']]
init_miss   = [f for f in fails if '__init__' in f['error'] and 'missing' in f['error']]
comp_miss   = [f for f in fails if 'compute()' in f['error'] and 'missing' in f['error']]
other       = [f for f in fails if f not in no_compute and f not in attr_dict and f not in init_miss and f not in comp_miss]

print(f'  no compute() method:         {len(no_compute)}')
print(f'  dict has no .attr:           {len(attr_dict)}')
print(f'  __init__ missing args:       {len(init_miss)}')
print(f'  compute() missing args:      {len(comp_miss)}')
print(f'  other:                       {len(other)}')
print()

print('--- dict.attr failures (compute() uses dataset.X not dataset[X]) ---')
for f in attr_dict:
    attr = f['error'].split("has no attribute '")[-1].rstrip("'")
    print(f'  {f["class"]:50s}  .{attr}')

print()
print('--- compute() missing positional args (old signature) ---')
for f in comp_miss:
    print(f'  {f["class"]:50s}  {f["error"][:100]}')

print()
print('--- __init__ missing args (dataclass/struct) ---')
for f in init_miss:
    print(f'  {f["class"]:50s}  {f["error"][:100]}')

print()
print('--- other ---')
for f in other:
    print(f'  {f["class"]:50s}  {f["error"][:100]}')

print()
print('--- no compute() (infrastructure classes — list for review) ---')
for f in no_compute:
    print(f'  {f["class"]}')
