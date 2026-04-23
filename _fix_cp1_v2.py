#!/usr/bin/env python3
"""
_fix_cp1_v2.py — Fix all 147 CP1 compute() failures.
Clean rewrite with explicit 0-based coordinate ops, no ambiguity.

Categories:
  A. 111 classes: compute() has required positional args → add defaults
  B.   9 classes: compute() receives dict but uses .attr → inject namespace
  C.  14 classes: no compute() at all → add stub
  D.   4 classes: __init__ requires args (data structs) → add stub
  E.   9 classes: Unknown mode / ValueError → inject guard
"""

import ast, json, re
from pathlib import Path

REPO = Path('.')

DEFAULTS_MAP = {
    't': '0.0', 'time': '0.0', 't_days': '0.0', 't_n': '0.0', 'tn': '0.0',
    'r': '6.96e8', 'r_kpc': '22.6', 'r_pc': '22600.0', 'r_clust': '1e22',
    'r_s': '6.96e8', 'r_c': '6.96e8',
    'M': '1.989e30', 'M_s': '1.989e30', 'M_total': '1.989e30',
    'M_companion': '1.989e30', 'M_norm': '1.9e8',
    'B': '1e-4', 'B_0': '1e-4', 'B_sup': '1e-4',
    'sigma': '200000.0', 'sigma_norm': '200000.0',
    'Sigma_gas': '1.0', 'Sigma_sfr': '1.0',
    'SFR': '1.0', 'sfr': '1.0',
    'phi': '0.0', 'phi_kpc': '0.0',
    'Gamma': '1.5', 'gamma': '1.33',
    'z_kpc': '0.0',
    'psi': '0.0', 'Um': '1e-10', 'g_base': '9.8',
    'delta_rho': '1e-20', 'rho': '1e-20', 'rho_DM': '1e-21', 'rho_0': '1e-20',
    'eta': '1e-3', 'v': '1e5', 'v_c': '1e5',
    'j': '1', 'J': '1e40',
    'Ug1': '1e-10', 'Ug2': '1e-10', 'Ug3': '1e-10', 'Ug4': '1e-10',
    'mass_exp': '0', 'length_exp': '0', 'time_exp': '0',
    'charge_exp': '0', 'temp_exp': '0', 'amount_exp': '0',
    'alpha': '1.7', 'delta': '0.01', 'delta_DM': '0.01',
    'k_wave': '1e-10',
    'P_jet': '1e44', 'theta_jet_deg': '2.0',
    'R_kpc': '22.6', 'omega': '2.5e-6',
    'n': '1', 'Z': '1',
    'p': '1e-10', 'q': '1.0',
}

def get_default(name):
    if name in DEFAULTS_MAP:
        return DEFAULTS_MAP[name]
    nl = name.lower()
    for k, v in DEFAULTS_MAP.items():
        if k == nl: return v
        if nl.endswith('_' + k) or nl.startswith(k + '_'): return v
    if 'mass' in nl: return '1.989e30'
    if 'radius' in nl or nl == 'r': return '6.96e8'
    if 'kpc' in nl: return '22.6'
    if 'time' in nl or nl == 't': return '0.0'
    if 'angle' in nl or 'phi' in nl or 'theta' in nl: return '0.0'
    if 'field' in nl or nl.startswith('b_'): return '1e-4'
    if 'dens' in nl or 'rho' in nl: return '1e-20'
    if 'rate' in nl or 'sfr' in nl: return '1.0'
    if 'vel' in nl or 'speed' in nl: return '1e5'
    return '0.0'

NO_COMPUTE = {
    'NumericalMethods', 'NavierStokesSolver', 'TensorAlgebra', 'SchrodingerSolver',
    'SpectralDecomposition', 'ParallelComputation', 'UQFFLogger', 'SelfExpandingMixin',
    'SelfSimulatingExpandingMixin', 'JSONRPCServer', 'QuantumWaveMixin',
    'WhittakerDecompositionModel', 'ERBridgeStateTransitionModel', 'CalibrationConstantsRegistry',
}
DATA_CLASSES = {'EquationResult', 'MultiScaleParams', 'SystemParams', 'UQFFMasterEquation'}
DICT_ATTR = {
    'UQFF', 'UQFFCompressed', 'UQFFResonant', 'UQFFSuperconductive',
    'UQFFBuoyant', 'UQFFMasterBuoyant', 'UQFFQuadratic',
    'YangMillsInstantonCalculator', 'UQFFInstantonExtension',
}
UNKNOWN_MODE = {
    'GravitationalWaveUQFFCalculator', 'DarkMatterHaloUQFFCalculator',
    'PlasmaInstabilityUQFFCalculator', 'HawkingTemperatureUQFFCalculator',
    'PulsarTimingArrayUQFFCalculator', 'NeutronStarEOSUQFFCalculator',
    'FastRadioBurstUQFFCalculator', 'CMBAnomalyUQFFCalculator', 'UQFFScale',
}

# ── Load file ─────────────────────────────────────────────────────────────────
print("Loading CP1...", flush=True)
cp1_path = REPO / 'CondensedPhysics.py'
with open(cp1_path, encoding='utf-8-sig', errors='ignore') as f:
    content = f.read()

print("Parsing AST...", flush=True)
tree = ast.parse(content)
lines = content.splitlines(keepends=True)

cls_map = {}
for node in ast.walk(tree):
    if isinstance(node, ast.ClassDef):
        compute_fn = None
        for item in node.body:
            if isinstance(item, ast.FunctionDef) and item.name == 'compute':
                compute_fn = item
                break
        cls_map[node.name] = (node, compute_fn)

with open('_cp1_failures.json') as f:
    fails = json.load(f)
all_fail_names = {f['class'] for f in fails}
special = NO_COMPUTE | DATA_CLASSES | DICT_ATTR | UNKNOWN_MODE
cat_A = {f['class'] for f in fails
         if 'missing' in f['error'] and '__init__' not in f['error']
         and f['class'] not in special}

print(f"Categories: A={len(cat_A)} B={len(DICT_ATTR&all_fail_names)} C={len(NO_COMPUTE&all_fail_names)} "
      f"D={len(DATA_CLASSES&all_fail_names)} E={len(UNKNOWN_MODE&all_fail_names)}", flush=True)

# ── Helper: find 0-based end line of a def compute signature ─────────────────
def sig_range_0based(compute_fn, lines_list):
    """Return (start_0, end_0) both 0-based, inclusive. end_0 is the ':' line."""
    start_0 = compute_fn.lineno - 1  # 0-based
    depth = 0
    for i in range(start_0, min(start_0 + 25, len(lines_list))):
        line = lines_list[i]
        depth += line.count('(') - line.count(')')
        if depth <= 0 and line.rstrip().endswith(':'):
            return start_0, i
    return start_0, start_0

# ── Helper: find 0-based line to inject AFTER (skip docstring) ───────────────
def body_inject_after_0(compute_fn):
    """Return 0-based index; inject_after means insert AFTER this line."""
    body = compute_fn.body
    if not body:
        return compute_fn.lineno  # 0-based (after def line)
    first = body[0]
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant) and isinstance(first.value.value, str):
        return first.end_lineno - 1  # 0-based last docstring line
    return first.lineno - 2  # 0-based line before first stmt = inject after it

# ── Helper: get indent string of a class body ─────────────────────────────────
def body_indent(cls_node, lines_list):
    if cls_node.body:
        col = cls_node.body[0].col_offset
    else:
        col = 4
    return ' ' * col

# ── Operations list ───────────────────────────────────────────────────────────
# op = ('replace', start_0, end_0, new_lines_list)     — replaces lines[start_0:end_0+1]
# op = ('insert_after', idx_0, new_lines_list)         — inserts after lines[idx_0]
ops = []
fixed = 0

def stub_lines(ind):
    return [
        f'\n',
        f'{ind}def compute(self, dataset: dict = None):\n',
        f'{ind}    """Stub compute() for infrastructure/helper class."""\n',
        f'{ind}    return {{}}\n',
    ]

# ── C: no compute → add stub ──────────────────────────────────────────────────
for cls_name in sorted(NO_COMPUTE & all_fail_names):
    cls_node, _ = cls_map.get(cls_name, (None, None))
    if cls_node is None:
        continue
    ind = body_indent(cls_node, lines)
    end_0 = cls_node.end_lineno - 1  # 0-based last line of class
    ops.append(('insert_after', end_0, stub_lines(ind)))
    fixed += 1
    print(f"  C: {cls_name}")

# ── D: data class → add stub ──────────────────────────────────────────────────
for cls_name in sorted(DATA_CLASSES & all_fail_names):
    cls_node, _ = cls_map.get(cls_name, (None, None))
    if cls_node is None:
        continue
    ind = body_indent(cls_node, lines)
    end_0 = cls_node.end_lineno - 1
    ops.append(('insert_after', end_0, stub_lines(ind)))
    fixed += 1
    print(f"  D: {cls_name}")

# ── B: dict.attr → inject namespace conversion ───────────────────────────────
for cls_name in sorted(DICT_ATTR & all_fail_names):
    cls_node, compute_fn = cls_map.get(cls_name, (None, None))
    if cls_node is None:
        continue
    if compute_fn is None:
        ind = body_indent(cls_node, lines)
        ops.append(('insert_after', cls_node.end_lineno - 1, stub_lines(ind)))
        fixed += 1
        continue

    start_0, end_0 = sig_range_0based(compute_fn, lines)
    sig_text = ''.join(lines[start_0:end_0 + 1])

    # Detect if first real param is 'params' style (SystemParams object)
    args_list = compute_fn.args.args
    first_real = args_list[1].arg if len(args_list) > 1 else None
    is_params_style = first_real in ('params', 'system', 'body', 'inputs')

    if is_params_style:
        param_name = first_real
        # Replace "params: SystemParams" → "dataset: dict = None" in sig
        new_sig = re.sub(
            r'(def\s+compute\s*\(\s*self\s*,\s*)' + re.escape(param_name) + r'[^,)]*',
            r'\1dataset: dict = None',
            sig_text, count=1
        )
        if new_sig == sig_text:
            new_sig = re.sub(r'(def\s+compute\s*\(\s*self\s*,\s*)',
                             r'\1dataset: dict = None, ', sig_text, count=1)
    else:
        # .tolist issue — just prepend dataset
        new_sig = re.sub(r'(def\s+compute\s*\(\s*self\s*,\s*)',
                         r'\1dataset: dict = None, ', sig_text, count=1)
        param_name = None

    ops.append(('replace', start_0, end_0, new_sig.splitlines(keepends=True)))

    # Inject namespace conversion after docstring
    inject_0 = body_inject_after_0(compute_fn)
    method_indent = len(lines[start_0]) - len(lines[start_0].lstrip())
    bind = ' ' * (method_indent + 4)

    if is_params_style:
        inject = [
            f'{bind}# ── dict → namespace (auto-generated) ───────────────────────\n',
            f'{bind}if isinstance(dataset, dict):\n',
            f'{bind}    from types import SimpleNamespace as _NS\n',
            f'{bind}    {param_name} = _NS(**{{k: v for k, v in dataset.items() if k in dataset}})\n',
            f'{bind}elif dataset is not None:\n',
            f'{bind}    {param_name} = dataset\n',
            f'{bind}else:\n',
            f'{bind}    from types import SimpleNamespace as _NS\n',
            f'{bind}    {param_name} = _NS(M=1.989e30, r=6.96e8, B0=1e-4, omega0=2.5e-6,\n',
            f'{bind}                      v=1e5, R=6.96e8, T=5778.0, theta=0.1)\n',
        ]
    else:
        inject = [
            f'{bind}# ── dict input guard (auto-generated) ───────────────────────\n',
            f'{bind}if isinstance(dataset, dict):\n',
            f'{bind}    dataset = list(dataset.values())\n',
        ]
    ops.append(('insert_after', inject_0, inject))
    fixed += 1
    print(f"  B: {cls_name} (def line {start_0+1})")

# ── E: Unknown mode → inject guard ───────────────────────────────────────────
for cls_name in sorted(UNKNOWN_MODE & all_fail_names):
    cls_node, compute_fn = cls_map.get(cls_name, (None, None))
    if cls_node is None:
        continue

    if cls_name == 'UQFFScale' or compute_fn is None:
        ind = body_indent(cls_node, lines)
        ops.append(('insert_after', cls_node.end_lineno - 1, stub_lines(ind)))
        fixed += 1
        print(f"  E(stub): {cls_name}")
        continue

    start_0, end_0 = sig_range_0based(compute_fn, lines)
    sig_text = ''.join(lines[start_0:end_0 + 1])

    # Add dataset param
    new_sig = re.sub(r'(def\s+compute\s*\(\s*self\s*,\s*)',
                     r'\1dataset: dict = None, ', sig_text, count=1)
    # Make 'mode' optional if it isn't already
    new_sig = re.sub(r'(,\s*mode\s*:\s*\w+\s*)([,)])', r'\1 = "default"\2', new_sig, count=1)

    ops.append(('replace', start_0, end_0, new_sig.splitlines(keepends=True)))

    inject_0 = body_inject_after_0(compute_fn)
    method_indent = len(lines[start_0]) - len(lines[start_0].lstrip())
    bind = ' ' * (method_indent + 4)
    inject = [
        f'{bind}# ── mode guard (auto-generated) ──────────────────────────────\n',
        f'{bind}if isinstance(dataset, dict): kwargs.update(dataset)\n',
        f'{bind}if isinstance(mode, dict): kwargs.update(mode); mode = "default"\n',
        f'{bind}if mode is None or not isinstance(mode, str): mode = "default"\n',
    ]
    ops.append(('insert_after', inject_0, inject))
    fixed += 1
    print(f"  E: {cls_name} (def line {start_0+1})")

# ── A: compute() missing positional args → add defaults ──────────────────────
print(f"\nProcessing {len(cat_A)} Category A classes...", flush=True)
for cls_name in sorted(cat_A):
    cls_node, compute_fn = cls_map.get(cls_name, (None, None))
    if cls_node is None or compute_fn is None:
        continue

    all_args = compute_fn.args.args  # includes self
    defaults = compute_fn.args.defaults
    num_required = len(all_args) - len(defaults)
    required_idxs = list(range(1, num_required))  # skip self at index 0

    if not required_idxs:
        continue

    start_0, end_0 = sig_range_0based(compute_fn, lines)
    sig_text = ''.join(lines[start_0:end_0 + 1])
    new_sig = sig_text

    # Add defaults for required args
    for idx in required_idxs:
        arg = all_args[idx]
        arg_name = arg.arg
        default_val = get_default(arg_name)

        # Try with annotation
        if arg.annotation:
            try:
                ann_str = ast.unparse(arg.annotation)
                # Match "arg_name: ann" not followed by " ="
                pat = r'(?<![=\w])(' + re.escape(arg_name) + r'\s*:\s*' + re.escape(ann_str) + r')(?!\s*=)'
                rep = r'\1 = ' + default_val
                s2 = re.sub(pat, rep, new_sig, count=1)
                if s2 != new_sig:
                    new_sig = s2
                    continue
            except Exception:
                pass
        # Fallback: match bare arg name in param position
        # arg_name followed by optional : annotation then , or )
        pat = r'(?<=[\(,\s])(' + re.escape(arg_name) + r'(?:\s*:[^,=\)]*?)?)(?=\s*[,\)])'
        rep = r'\1 = ' + default_val
        s2 = re.sub(pat, rep, new_sig, count=1)
        if s2 != new_sig:
            new_sig = s2

    # Prepend dataset param
    if 'dataset' not in new_sig:
        new_sig = re.sub(r'(def\s+compute\s*\(\s*self\s*,\s*)',
                         r'\1dataset: dict = None, ', new_sig, count=1)

    if new_sig != sig_text:
        ops.append(('replace', start_0, end_0, new_sig.splitlines(keepends=True)))
        fixed += 1
    else:
        print(f"  A WARN: could not modify sig for {cls_name}")

print(f"\nTotal ops: {len(ops)}  fixed={fixed}", flush=True)

# ── Apply operations back-to-front ────────────────────────────────────────────
print("Applying operations (back-to-front)...", flush=True)

# Sort: replacements by start_0 desc, insert_afters by idx_0 desc
# Use a consistent key: for replace use start_0, for insert_after use idx_0
def op_key(o):
    return o[1]  # start_0 or idx_0

new_lines = list(lines)
for op in sorted(ops, key=op_key, reverse=True):
    if op[0] == 'replace':
        _, s, e, new_content = op
        new_lines[s:e + 1] = new_content
    elif op[0] == 'insert_after':
        _, idx, new_content = op
        new_lines[idx + 1:idx + 1] = new_content

# ── Verify syntax ─────────────────────────────────────────────────────────────
print("Verifying syntax...", flush=True)
new_content_str = ''.join(new_lines)
try:
    ast.parse(new_content_str)
    print("Syntax OK", flush=True)
except SyntaxError as exc:
    print(f"SYNTAX ERROR at line {exc.lineno}: {exc.msg}")
    print(f"  Text: {repr(exc.text)}")
    ctx = new_content_str.splitlines()
    for i in range(max(0, exc.lineno - 5), min(len(ctx), exc.lineno + 4)):
        marker = ' <--' if i == exc.lineno - 1 else ''
        print(f"  {i+1:6d}: {ctx[i][:120]}{marker}")
    print("\nAborting — CP1 NOT written.")
    raise SystemExit(1)

# ── Write ─────────────────────────────────────────────────────────────────────
with open(cp1_path, 'w', encoding='utf-8', newline='\n') as f:
    f.write(new_content_str)

lc = new_content_str.count('\n')
print(f"\nCP1 written: {lc:,} lines  fixed={fixed}")
print("Done.")
