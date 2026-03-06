import re, hashlib
with open('CondensedPhysics.py', encoding='utf-8', errors='replace') as f:
    lines = f.read().splitlines()
results = []
for i, line in enumerate(lines):
    if re.match(r'\s+def compute_master_buoyant_equation', line):
        cls = 'UNKNOWN'
        for j in range(i-1,-1,-1):
            m = re.match(r'^class (\w+)', lines[j])
            if m:
                cls = m.group(1); break
        indent = len(line)-len(line.lstrip())
        body = [line]
        for k in range(i+1, min(i+60,len(lines))):
            l2 = lines[k]
            if l2.strip()=='':
                body.append(l2); continue
            li2 = len(l2)-len(l2.lstrip())
            if li2<=indent and re.match(r'\s*(def |class )',l2): break
            body.append(l2)
        norm = '\n'.join(l.rstrip() for l in body)
        h = hashlib.md5(norm.encode()).hexdigest()[:8]
        results.append({'line':i+1,'cls':cls,'hash':h})

v2 = [r for r in results if r['hash']=='928fb736']
v3 = [r for r in results if r['hash']=='27402435']
print('VARIANT 2 classes (Form C-2, TRZ resonant):')
for r in v2:
    print('  line %6d  %s' % (r['line'], r['cls']))
print()
print('VARIANT 3 classes (Form C-1, galactic Omega_g):')
for r in v3:
    print('  line %6d  %s' % (r['line'], r['cls']))
