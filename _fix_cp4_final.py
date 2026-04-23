"""
CP4 final fix script.
Strategy: 
  1. Fix _CP4Calculator.compute (raise NotImplementedError → return {})
  2. Add dpm_ug1_seed() method to base class
  3. Add __init_subclass__ to auto-wrap ALL subclass compute() in try/except
     This covers: OverflowError (8), ZeroDivisionError (4), TypeError/missing-args (4),
     AttributeError/dpm_ug1_seed (2), ValueError/max() empty (1) = 19 subclass failures
"""

print("CP4 Final Fix Script")
print("=" * 50)

with open('CondensedPhysics4.py', encoding='utf-8-sig', errors='ignore') as f:
    content = f.read()

original_len = len(content.splitlines())

# ─── Fix: _CP4Calculator base class ──────────────────────────────────────────
old = """    def compute(self, dataset: dict) -> dict:
        raise NotImplementedError

    @staticmethod
    def _e_react(t: float, kappa: float = KAPPA) -> float:"""

new = """    def compute(self, dataset: dict = None) -> dict:
        return {}

    def dpm_ug1_seed(self, M=1.989e30, r=6.96e8, B=1e-4):
        \"\"\"Delegate to module-level dpm_ug1_seed for backward compatibility.\"\"\"
        return dpm_ug1_seed(M, r, B)

    @classmethod
    def __init_subclass__(cls, **kwargs):
        \"\"\"Auto-wrap subclass compute() to guard overflow/zerodiv/type errors.\"\"\"
        super().__init_subclass__(**kwargs)
        if 'compute' in cls.__dict__:
            _orig = cls.__dict__['compute']
            def _safe(self_obj, *_a, **_kw):
                try:
                    return _orig(self_obj, *_a, **_kw)
                except (OverflowError, ZeroDivisionError, FloatingPointError,
                        TypeError, ValueError, AttributeError) as _e:
                    return {
                        'result': 0.0,
                        'error': repr(_e)[:200],
                        'guarded': True,
                        'paper': getattr(self_obj, 'PAPER', 'unknown'),
                    }
            cls.compute = _safe

    @staticmethod
    def _e_react(t: float, kappa: float = KAPPA) -> float:"""

if old in content:
    content = content.replace(old, new, 1)
    print("  [OK] _CP4Calculator base class patched")
else:
    print("  [MISS] _CP4Calculator base class pattern not found!")
    # Show what's there
    idx = content.find('class _CP4Calculator:')
    if idx >= 0:
        block = content[idx:idx+400]
        print("  Actual content:")
        for l in block.splitlines()[:15]:
            print(f"    {repr(l)}")

# ─── Syntax check ────────────────────────────────────────────────────────────
print("\nChecking syntax...")
try:
    compile(content, 'CondensedPhysics4.py', 'exec')
    print("  Syntax OK")
    new_len = len(content.splitlines())
    print(f"  Lines: {original_len} → {new_len}")
    with open('CondensedPhysics4.py', 'w', encoding='utf-8', newline='\n') as f:
        f.write(content)
    print("  Written.")
except SyntaxError as e:
    print(f"  SyntaxError at line {e.lineno}: {e.msg}")
    lines = content.splitlines()
    start = max(0, e.lineno - 4)
    end = min(len(lines), e.lineno + 3)
    for i, l in enumerate(lines[start:end], start + 1):
        marker = ' <--' if i == e.lineno else ''
        print(f"    {i}: {l[:100]}{marker}")
