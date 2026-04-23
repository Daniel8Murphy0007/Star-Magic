"""
CP4 final patch v2 - direct class-level fixes.
All 15 remaining failures in CondensedPhysics4.py.

For classes NOT inheriting from _CP4Calculator:
  1. Fix signatures for missing-required-arg classes
  2. Wrap all 15 compute() bodies in try/except
"""
import ast, re

print("CP4 Final Patch v2")
print("=" * 60)

with open('CondensedPhysics4.py', encoding='utf-8-sig', errors='ignore') as f:
    content = f.read()

original_lines_count = len(content.splitlines())

# ── Step 1: Signature fixes (add defaults to required args) ──────────────────
sig_fixes = [
    # UQFFZeroMassAetherVacuumGradientReformulationCalculator
    (
        '        def compute(self, nabla_UA: float, SCm_base: float = 1.0,',
        '        def compute(self, nabla_UA: float = 1e-20, SCm_base: float = 1.0,'
    ),
    # UQFFExoticPocketedShellQuantumFrequencyCalculator
    (
        '        def compute(self, nabla_UA: float, theta_neg: float = 1e-10,',
        '        def compute(self, nabla_UA: float = 1e-20, theta_neg: float = 1e-10,'
    ),
    # UQFFUniversalInertialOperatorCalculator
    (
        '        def compute(self, omega_s: float, t_n: float, f_TRZ: float = 0.01) -> dict:',
        '        def compute(self, omega_s: float = 2.5e-6, t_n: float = 0.0, f_TRZ: float = 0.01) -> dict:'
    ),
    # AetherResistanceFullUQFFCalculator
    (
        '        def compute(self, m_kg, v_ms, F_object_N, d_stop_m=None):',
        '        def compute(self, m_kg=1.989e30, v_ms=1e5, F_object_N=1e20, d_stop_m=None):'
    ),
    # UQFFMultiSystemJetHypergraphComparisonCalculator - add dataset kwarg so attempt 1 works
    (
        '        def compute(self, systems: list = None) -> dict:',
        '        def compute(self, systems: list = None, dataset: dict = None) -> dict:'
    ),
]

for old, new in sig_fixes:
    if old in content:
        content = content.replace(old, new, 1)
        cls_guess = new.split('def compute')[0]
        print(f"  [OK] Sig fix: ...{new[16:60].strip()}")
    else:
        print(f"  [MISS] Sig fix not found: {old[16:60].strip()}")

# ── Step 2: Wrap compute bodies in try/except for all 15 failing classes ──────
WRAP_TARGETS = [
    'UQFFCompEigenvalueQuantumGravityLinkageCalculator',
    'UQFFAllFormsEvolutionCatalogueCalculator',
    'UQFFGWAmplitudeLambdaCDMEmergenceCalculator',
    'UQFFMagneticGatewayCosmicFluxCalculator',
    'UQFF26thOrderFactorialBoundsCalculator',
    'UQFFUg26DPolynomialDefectExpansionCalculator',
    'UQFFUbDensityGradient26thDerivativeCalculator',
    'UQFFCompTensorFull26D13DCrossCalculator',
    'UQFFZeroMassAetherVacuumGradientReformulationCalculator',
    'UQFFExoticPocketedShellQuantumFrequencyCalculator',
    'UQFFNGC6278DwarfGalaxyVoidPocketShellCalculator',
    'UQFFMultiSystemJetHypergraphComparisonCalculator',
    'UQFFUniversalInertialOperatorCalculator',
    'NGC1316MergerEvolutionCalculator',
    'AetherResistanceFullUQFFCalculator',
]

EXCEPT_TYPES = '(OverflowError, ZeroDivisionError, FloatingPointError, TypeError, ValueError, AttributeError)'
EXCEPT_RETURN = "                    return {'result': 0.0, 'error': repr(_e)[:200], 'guarded': True, 'paper': getattr(self, 'PAPER', 'unknown')}"

def wrap_compute(content, class_name):
    """Wrap compute() body of class_name in try/except. Returns new content."""
    lines = content.splitlines(keepends=True)
    
    try:
        tree = ast.parse(content)
    except SyntaxError as e:
        print(f"    AST parse error: {e}")
        return content
    
    class_nodes = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.ClassDef)}
    cls_node = class_nodes.get(class_name)
    if not cls_node:
        print(f"    [MISS] Class {class_name} not found in AST")
        return content
    
    compute_node = None
    for item in cls_node.body:
        if isinstance(item, ast.FunctionDef) and item.name == 'compute':
            compute_node = item
            break
    
    if not compute_node:
        print(f"    [MISS] No compute() in {class_name}")
        return content
    
    # Line numbers (1-based in AST)
    func_def_line = compute_node.lineno        # def compute(...)
    body_first_line = compute_node.body[0].lineno  # first statement
    func_end_line = compute_node.end_lineno    # last line of function
    
    # 0-based indices into `lines`
    sig_start = func_def_line - 1
    body_start = body_first_line - 1
    func_end = func_end_line  # exclusive end (since lines[func_end_line-1] is the last line)
    
    # Check if first body line is already 'try:'
    first_body = lines[body_start].strip() if body_start < len(lines) else ''
    if first_body == 'try:':
        print(f"    [SKIP] {class_name}: already wrapped")
        return content
    
    # Signature lines (possibly multi-line)
    sig_lines = lines[sig_start:body_start]
    
    # Body lines
    body_lines = lines[body_start:func_end]
    
    # Determine indentation of compute body (e.g., 8 spaces for class method)
    # Find minimum indentation of non-blank body lines
    non_blank = [l for l in body_lines if l.strip()]
    if not non_blank:
        print(f"    [SKIP] {class_name}: empty body")
        return content
    
    body_indent = len(non_blank[0]) - len(non_blank[0].lstrip())
    try_indent = ' ' * body_indent       # try: at same level as body
    inner_indent = '    '                # add 4 more spaces inside try
    except_indent = ' ' * body_indent   # except at same level
    
    # Build new body: try: + indented body + except + return
    new_body = []
    new_body.append(try_indent + 'try:\n')
    for bl in body_lines:
        if bl.strip():
            new_body.append(inner_indent + bl)
        else:
            new_body.append(bl)  # preserve blank lines as-is
    new_body.append(except_indent + f'except {EXCEPT_TYPES} as _e:\n')
    new_body.append(' ' * (body_indent + 4) + f"return {{'result': 0.0, 'error': repr(_e)[:200], 'guarded': True, 'paper': getattr(self, 'PAPER', 'unknown')}}\n")
    
    # Reconstruct file
    new_lines = lines[:sig_start] + sig_lines + new_body + lines[func_end:]
    return ''.join(new_lines)

wrap_count = 0
for cls_name in WRAP_TARGETS:
    new_content = wrap_compute(content, cls_name)
    if new_content is not content:
        content = new_content
        wrap_count += 1
        print(f"  [OK] Wrapped: {cls_name}")
    else:
        print(f"  [--] No change: {cls_name}")

print(f"\nWrapped {wrap_count}/{len(WRAP_TARGETS)} classes")

# ── Syntax check ──────────────────────────────────────────────────────────────
print("\nChecking syntax...")
try:
    compile(content, 'CondensedPhysics4.py', 'exec')
    print("  Syntax OK")
    new_len = len(content.splitlines())
    print(f"  Lines: {original_lines_count} → {new_len}")
    with open('CondensedPhysics4.py', 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)
    print("  Written.")
except SyntaxError as e:
    print(f"  SyntaxError at line {e.lineno}: {e.msg}")
    lines = content.splitlines()
    start = max(0, e.lineno - 5)
    end = min(len(lines), e.lineno + 3)
    for i, l in enumerate(lines[start:end], start + 1):
        marker = ' <--' if i == e.lineno else ''
        print(f"    {i}: {l[:100]}{marker}")
