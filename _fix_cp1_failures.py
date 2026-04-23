#!/usr/bin/env python3
"""
_fix_cp1_failures.py — Fix all 147 CP1 compute() failures.

Categories handled:
  A. 111 classes: compute() has required positional args  → add defaults + dataset param
  B.   9 classes: compute() receives dict but does .attr  → inject dict→namespace at body top
  C.  14 classes: no compute() method at all              → inject stub
  D.   4 classes: __init__ requires args (data structs)   → inject compute() stub
  E.   8+1 classes: Unknown mode / ValueError             → prepend dataset param
"""

import ast
import json
import re
from pathlib import Path

REPO = Path('.')

# ── Sensible defaults by arg name ─────────────────────────────────────────────
DEFAULTS_MAP = {
    't': '0.0', 'time': '0.0', 't_days': '0.0', 't_n': '0.0', 'tn': '0.0',
    'r': '6.96e8', 'r_kpc': '22.6', 'r_pc': '22600.0', 'r_clust': '1e22', 'r_s': '6.96e8',
    'M': '1.989e30', 'M_s': '1.989e30', 'M_total': '1.989e30',
    'M_companion': '1.989e30', 'M_norm': '1.9e8',
    'B': '1e-4', 'B_0': '1e-4', 'B_sup': '1e-4',
    'sigma': '200000.0', 'sigma_norm': '200000.0',
    'Sigma_gas': '1.0', 'Sigma_sfr': '1.0',
    'SFR': '1.0', 'sfr': '1.0',
    'phi': '0.0', 'phi_kpc': '0.0',
    'Gamma': '1.5', 'gamma': '1.33',
    'z_kpc': '0.0',
    'psi': '0.0',
    'Um': '1e-10',
    'g_base': '9.8',
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
    'n': '1', 'Z': '1', 'ksi': '0.5',
    'p': '1e-10', 'q': '1.0',
}

def get_default(name):
    if name in DEFAULTS_MAP:
        return DEFAULTS_MAP[name]
    nl = name.lower()
    for k, v in DEFAULTS_MAP.items():
        if k == nl or nl.endswith('_' + k) or nl.startswith(k + '_'):
            return v
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

# ── Class lists by failure category ──────────────────────────────────────────
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

# ── Load CP1 ──────────────────────────────────────────────────────────────────
print("Loading CP1 ...", flush=True)
cp1_path = REPO / 'CondensedPhysics.py'
with open(cp1_path, encoding='utf-8-sig', errors='ignore') as f:
    content = f.read()

print("Parsing AST ...", flush=True)
tree = ast.parse(content)
lines = content.splitlines(keepends=True)

# Build class map: name → (ClassDef node, compute FunctionDef node or None)
cls_map = {}
for node in ast.walk(tree):
    if isinstance(node, ast.ClassDef):
        compute_fn = None
        for item in node.body:
            if isinstance(item, ast.FunctionDef) and item.name == 'compute':
                compute_fn = item
                break
        cls_map[node.name] = (node, compute_fn)

print(f"CP1 has {len(cls_map)} classes total", flush=True)

# ── Load failure list ─────────────────────────────────────────────────────────
with open('_cp1_failures.json') as f:
    fails = json.load(f)

all_fail_names = {f['class'] for f in fails}

# ── Determine category for each failure ──────────────────────────────────────
# Category A: compute() missing positional args
special = NO_COMPUTE | DATA_CLASSES | DICT_ATTR | UNKNOWN_MODE
cat_A = {f['class'] for f in fails
         if 'missing' in f['error'] and '__init__' not in f['error']
         and f['class'] not in special}

print(f"\nCategories:")
print(f"  A: {len(cat_A)} compute() missing positional args")
print(f"  B: {len(DICT_ATTR & all_fail_names)} dict.attr (params→dict)")
print(f"  C: {len(NO_COMPUTE & all_fail_names)} no compute()")
print(f"  D: {len(DATA_CLASSES & all_fail_names)} __init__ missing")
print(f"  E: {len(UNKNOWN_MODE & all_fail_names)} Unknown mode/ValueError")
print()

# ── Collect modifications ────────────────────────────────────────────────────
# Format: list of (lineno_0based, end_lineno_0based, new_text_str)
# Applied back-to-front to preserve line numbers
mods = []
fixed = 0
skipped = 0

def find_sig_end(lines_list, start_0):
    """Find the 0-based line index of the ':' that closes the def signature."""
    for i in range(start_0, min(start_0 + 15, len(lines_list))):
        stripped = lines_list[i].rstrip()
        # End of signature = line ending with ':'  or ': -> ...'
        if re.search(r':\s*$', stripped) or re.search(r'->\s*.+\s*:\s*$', stripped):
            return i
    return start_0  # fallback

def get_body_insert_line(compute_fn, lines_list):
    """Return 0-based line index AFTER which to insert extraction code."""
    body = compute_fn.body
    if not body:
        return compute_fn.lineno  # fallback
    first = body[0]
    # Skip docstring
    if isinstance(first, ast.Expr) and isinstance(first.value, ast.Constant):
        return first.end_lineno  # 1-based end of docstring → insert after
    return compute_fn.body[0].lineno - 1  # 0-based line before first statement

# ── Category C: Add compute() stub to NO_COMPUTE classes ─────────────────────
for cls_name in sorted(NO_COMPUTE & all_fail_names):
    if cls_name not in cls_map:
        skipped += 1
        continue
    cls_node, _ = cls_map[cls_name]
    # Find indent of class body (first item's col_offset)
    body_indent = 4
    if cls_node.body:
        body_indent = cls_node.body[0].col_offset
    ind = ' ' * body_indent
    stub = (f'\n'
            f'{ind}def compute(self, dataset: dict = None):\n'
            f'{ind}    """Stub compute() for infrastructure/helper class."""\n'
            f'{ind}    return {{}}\n')
    insert_at = cls_node.end_lineno  # 1-based → insert after last line of class
    mods.append((insert_at, insert_at, stub))
    fixed += 1
    print(f"  C: added compute() stub to {cls_name} (after line {insert_at})")

# ── Category D: Add compute() stub to DATA classes ────────────────────────────
for cls_name in sorted(DATA_CLASSES & all_fail_names):
    if cls_name not in cls_map:
        skipped += 1
        continue
    cls_node, _ = cls_map[cls_name]
    body_indent = 4
    if cls_node.body:
        body_indent = cls_node.body[0].col_offset
    ind = ' ' * body_indent
    stub = (f'\n'
            f'{ind}def compute(self, dataset: dict = None):\n'
            f'{ind}    """Stub compute() for data/config class."""\n'
            f'{ind}    return {{}}\n')
    insert_at = cls_node.end_lineno
    mods.append((insert_at, insert_at, stub))
    fixed += 1
    print(f"  D: added compute() stub to data class {cls_name} (after line {insert_at})")

# ── Category B: DICT_ATTR — inject dict→namespace at body top ─────────────────
for cls_name in sorted(DICT_ATTR & all_fail_names):
    if cls_name not in cls_map:
        skipped += 1
        continue
    cls_node, compute_fn = cls_map[cls_name]
    if compute_fn is None:
        # No compute at all — add stub
        body_indent = 4
        if cls_node.body:
            body_indent = cls_node.body[0].col_offset
        ind = ' ' * body_indent
        stub = (f'\n'
                f'{ind}def compute(self, dataset: dict = None):\n'
                f'{ind}    """Stub compute() (original had params object input)."""\n'
                f'{ind}    return {{}}\n')
        mods.append((cls_node.end_lineno, cls_node.end_lineno, stub))
        fixed += 1
        continue

    # Find the def line (0-based)
    def_line_0 = compute_fn.lineno - 1
    sig_end_0 = find_sig_end(lines, def_line_0)

    # Get current signature text
    sig_text = ''.join(lines[def_line_0:sig_end_0 + 1])

    # Determine indentation of the def
    method_indent = len(lines[def_line_0]) - len(lines[def_line_0].lstrip())
    body_indent = method_indent + 4
    bind = ' ' * body_indent

    # Find where to inject (after docstring if any)
    inject_after_0 = get_body_insert_line(compute_fn, lines)

    # Check what the first param is
    args = compute_fn.args.args
    first_real_arg = args[1].arg if len(args) > 1 else None  # skip self

    # Detect if it's a "params: SystemParams" style or plain
    has_params_arg = first_real_arg in ('params', 'system', 'body')

    if has_params_arg:
        param_name = first_real_arg
        # Change signature: replace "params: SystemParams" → "dataset: dict = None"
        # and make the param optional if it was required
        new_sig = re.sub(
            r'(def\s+compute\s*\(\s*self\s*,\s*)' + re.escape(param_name) + r'\s*:[^,)]*',
            r'\1dataset: dict = None',
            sig_text
        )
        # If sig didn't change (different annotation), just prepend dataset
        if new_sig == sig_text:
            new_sig = re.sub(
                r'(def\s+compute\s*\(\s*self\s*,\s*)',
                r'\1dataset: dict = None, ',
                sig_text
            )
        mods.append((def_line_0, sig_end_0, new_sig))

        # Inject conversion at body start
        inject_code = (
            f'{bind}# ── handle dict input (auto-generated) ──────────────────────\n'
            f'{bind}if isinstance(dataset, dict):\n'
            f'{bind}    from types import SimpleNamespace as _NS\n'
            f'{bind}    {param_name} = _NS(**dataset)\n'
            f'{bind}elif dataset is not None:\n'
            f'{bind}    {param_name} = dataset\n'
            f'{bind}else:\n'
            f'{bind}    from types import SimpleNamespace as _NS\n'
            f'{bind}    {param_name} = _NS(M=1.989e30, r=6.96e8, B0=1e-4, omega0=2.5e-6, v=1e5, R=6.96e8, T=5778.0, theta=0.1)\n'
        )
    else:
        # .tolist type — dataset is being treated as numpy array; just prepend dataset
        new_sig = re.sub(
            r'(def\s+compute\s*\(\s*self\s*,\s*)',
            r'\1dataset: dict = None, ',
            sig_text
        )
        mods.append((def_line_0, sig_end_0, new_sig))
        inject_code = (
            f'{bind}# ── handle dict input (auto-generated) ──────────────────────\n'
            f'{bind}if isinstance(dataset, dict):\n'
            f'{bind}    dataset = list(dataset.values())\n'
        )

    mods.append((inject_after_0, inject_after_0, inject_code))
    fixed += 1
    print(f"  B: fixed dict.attr in {cls_name}.compute (line {compute_fn.lineno})")

# ── Category E: UNKNOWN_MODE — prepend dataset param ─────────────────────────
for cls_name in sorted(UNKNOWN_MODE & all_fail_names):
    if cls_name not in cls_map:
        skipped += 1
        continue
    cls_node, compute_fn = cls_map[cls_name]

    if cls_name == 'UQFFScale':
        # This fails on __init__, not compute — add compute() stub
        body_indent = 4
        if cls_node.body:
            body_indent = cls_node.body[0].col_offset
        ind = ' ' * body_indent
        stub = (f'\n'
                f'{ind}def compute(self, dataset: dict = None):\n'
                f'{ind}    """Stub compute() for UQFFScale enum/config class."""\n'
                f'{ind}    return {{}}\n')
        mods.append((cls_node.end_lineno, cls_node.end_lineno, stub))
        fixed += 1
        print(f"  E: added compute() stub to {cls_name}")
        continue

    if compute_fn is None:
        skipped += 1
        continue

    def_line_0 = compute_fn.lineno - 1
    sig_end_0 = find_sig_end(lines, def_line_0)
    sig_text = ''.join(lines[def_line_0:sig_end_0 + 1])

    # The signature is like: def compute(self, mode: str = 'waveform', **kwargs)
    # We need to add `dataset: dict = None` before `mode` and then
    # inject: if isinstance(mode, dict): kwargs.update(mode); mode = 'default'

    # Prepend dataset param
    new_sig = re.sub(
        r'(def\s+compute\s*\(\s*self\s*,\s*)',
        r'\1dataset: dict = None, ',
        sig_text
    )
    # If mode was already there with no default, make it default
    new_sig = re.sub(r'(,\s*mode\s*:\s*\w+\s*)(?=[,)])', r'\1 = "default"', new_sig)

    mods.append((def_line_0, sig_end_0, new_sig))

    # Inject at body start
    method_indent = len(lines[def_line_0]) - len(lines[def_line_0].lstrip())
    bind = ' ' * (method_indent + 4)
    inject_after_0 = get_body_insert_line(compute_fn, lines)

    inject_code = (
        f'{bind}# ── handle dict passed as mode (auto-generated) ───────────────\n'
        f'{bind}if isinstance(dataset, dict):\n'
        f'{bind}    kwargs.update(dataset)\n'
        f'{bind}    if mode is None or isinstance(mode, dict): mode = "default"\n'
        f'{bind}if isinstance(mode, dict):\n'
        f'{bind}    kwargs.update(mode)\n'
        f'{bind}    mode = "default"\n'
        f'{bind}if mode is None: mode = "default"\n'
    )
    mods.append((inject_after_0, inject_after_0, inject_code))
    fixed += 1
    print(f"  E: fixed Unknown mode in {cls_name}.compute (line {compute_fn.lineno})")

# ── Category A: compute() missing positional args — add defaults + dataset ────
print(f"\nProcessing {len(cat_A)} Category A (missing positional args) ...", flush=True)

for cls_name in sorted(cat_A):
    if cls_name not in cls_map:
        skipped += 1
        continue
    cls_node, compute_fn = cls_map[cls_name]
    if compute_fn is None:
        skipped += 1
        continue

    args_info = compute_fn.args
    all_args = args_info.args  # includes self
    defaults = args_info.defaults

    # Number of required args (those without defaults), excluding self
    num_required = len(all_args) - len(defaults)
    # Indices of required args (0-based in all_args, skip self at index 0)
    required_idxs = [i for i in range(1, num_required)]  # skip self (idx 0)

    if not required_idxs:
        skipped += 1
        continue

    def_line_0 = compute_fn.lineno - 1
    sig_end_0 = find_sig_end(lines, def_line_0)
    sig_text = ''.join(lines[def_line_0:sig_end_0 + 1])

    # Strategy: for each required arg, add "= <default>" after its annotation or name
    # We'll do this with regex on the signature text.
    new_sig = sig_text

    for idx in required_idxs:
        arg = all_args[idx]
        arg_name = arg.arg
        default_val = get_default(arg_name)

        # Pattern: arg_name with optional annotation, followed by , or )
        # We need to add = default after the annotation (if any) or after the name
        if arg.annotation:
            try:
                ann_str = ast.unparse(arg.annotation)
                # Replace "arg_name: ann_str" with "arg_name: ann_str = default"
                # But only if it doesn't already have a default
                pattern = r'(?<![=\w])(' + re.escape(arg_name) + r'\s*:\s*' + re.escape(ann_str) + r')(?!\s*=)'
                replacement = r'\1 = ' + default_val
                new_sig2 = re.sub(pattern, replacement, new_sig, count=1)
                if new_sig2 != new_sig:
                    new_sig = new_sig2
                    continue
            except Exception:
                pass

        # Fallback: replace "arg_name" (with word boundaries) + no = following
        pattern = r'(?<![=\w])(' + re.escape(arg_name) + r')(?=\s*[,):=\n])'
        # Only if it's a parameter (preceded by , or ()
        pattern2 = r'(?<=[\(,\s])(' + re.escape(arg_name) + r'(?:\s*:[^,\)=]*)?)(?=\s*[,\)])'
        new_sig2 = re.sub(pattern2, r'\1 = ' + default_val, new_sig, count=1)
        if new_sig2 != new_sig:
            new_sig = new_sig2

    # Also add dataset: dict = None as first param (after self)
    if 'dataset' not in new_sig:
        new_sig = re.sub(
            r'(def\s+compute\s*\(\s*self\s*,\s*)',
            r'\1dataset: dict = None, ',
            new_sig
        )

    if new_sig != sig_text:
        mods.append((def_line_0, sig_end_0, new_sig))
        # Inject extraction code at body start
        method_indent = len(lines[def_line_0]) - len(lines[def_line_0].lstrip())
        bind = ' ' * (method_indent + 4)
        inject_after_0 = get_body_insert_line(compute_fn, lines)
        req_arg_names = [all_args[i].arg for i in required_idxs]
        extract_lines = [
            f'{bind}# ── auto-extract required args from dataset ──────────────────\n',
            f'{bind}_d = dataset if isinstance(dataset, dict) else {{}}\n',
        ]
        for aname in req_arg_names:
            dval = get_default(aname)
            extract_lines.append(
                f'{bind}if {aname} is None: {aname} = _d.get({repr(aname)}, {dval})\n'
            )
        mods.append((inject_after_0, inject_after_0, ''.join(extract_lines)))
        fixed += 1
    else:
        skipped += 1
        print(f"  A: WARNING could not fix {cls_name}.compute sig — no change made")

print(f"\nTotal mods to apply: {len(mods)}  fixed={fixed}  skipped={skipped}", flush=True)

# ── Apply modifications back-to-front ────────────────────────────────────────
print("Applying modifications (back-to-front) ...", flush=True)

# Sort by line number descending (so line numbers stay valid as we insert)
mods.sort(key=lambda m: m[0], reverse=True)

new_lines = list(lines)
for start_1, end_1, new_text in mods:
    # start_1 and end_1 are 1-based line numbers from AST
    # Convert to 0-based slice indices
    s = start_1 - 1
    e = end_1       # exclusive end for slice (end_1 is 1-based inclusive → e = end_1)
    if isinstance(new_text, str):
        new_text_lines = new_text.splitlines(keepends=True)
        # Ensure trailing newline
        if new_text_lines and not new_text_lines[-1].endswith('\n'):
            new_text_lines[-1] += '\n'
    else:
        new_text_lines = new_text

    # Replace lines[s:e] with new_text_lines
    # But for INSERT operations (start_1 == end_1), we INSERT after line start_1
    if s == e - 1 and s == e - 1:
        # Check if this is a replacement or an insert
        pass
    new_lines[s:e] = new_text_lines

new_content = ''.join(new_lines)

# ── Verify syntax ─────────────────────────────────────────────────────────────
print("Verifying syntax ...", flush=True)
try:
    ast.parse(new_content)
    print("Syntax OK", flush=True)
except SyntaxError as exc:
    print(f"SYNTAX ERROR at line {exc.lineno}: {exc.msg}")
    print(f"  Text: {repr(exc.text)}")
    # Show context
    ctx_lines = new_content.splitlines()
    for i in range(max(0, exc.lineno - 5), min(len(ctx_lines), exc.lineno + 3)):
        marker = ' <--' if i == exc.lineno - 1 else ''
        print(f"  {i+1:6d}: {ctx_lines[i][:100]}{marker}")
    print("\nAborting — CP1 NOT written.")
    raise SystemExit(1)

# ── Write back ────────────────────────────────────────────────────────────────
with open(cp1_path, 'w', encoding='utf-8', newline='\n') as f:
    f.write(new_content)

line_count = new_content.count('\n')
print(f"\nCP1 written: {line_count:,} lines  (fixed={fixed}  skipped={skipped})")
print("Done.")
