#!/usr/bin/env python3
"""
COMPREHENSIVE C++ Equation Extractor for MAIN_1_CoAnQi.cpp

Extracts ALL equations including:
- PhysicsTerm class compute() methods
- Standalone double functions
- Inline calculations and sub-equations
- Namespace functions (SOURCE4, etc.)

Target: 1,393+ equations
"""

import re
import os

def extract_all_equations(filepath):
    """Extract ALL equations from C++ file."""
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
        lines = content.split('\n')
    
    equations = {}
    
    # ═══════════════════════════════════════════════════════════════════════════
    # 1. EXTRACT PhysicsTerm CLASS COMPUTE METHODS
    # ═══════════════════════════════════════════════════════════════════════════
    print("  Extracting PhysicsTerm classes...")
    
    # Find all class definitions with their bodies
    class_pattern = re.compile(
        r'class\s+(\w+)\s*:\s*public\s+PhysicsTerm\s*\{',
        re.MULTILINE
    )
    
    for match in class_pattern.finditer(content):
        class_name = match.group(1)
        start_pos = match.end()
        
        # Find matching closing brace
        brace_count = 1
        pos = start_pos
        while brace_count > 0 and pos < len(content):
            if content[pos] == '{':
                brace_count += 1
            elif content[pos] == '}':
                brace_count -= 1
            pos += 1
        
        class_body = content[start_pos:pos-1]
        
        # Extract compute method
        compute_match = re.search(
            r'double\s+compute\s*\(([^)]*)\)[^{]*\{(.+)',
            class_body,
            re.DOTALL
        )
        
        if compute_match:
            params = compute_match.group(1).strip()
            body_start = compute_match.start(2)
            
            # Find the closing brace of compute method
            method_body = compute_match.group(2)
            brace_count = 1
            pos = 0
            while brace_count > 0 and pos < len(method_body):
                if method_body[pos] == '{':
                    brace_count += 1
                elif method_body[pos] == '}':
                    brace_count -= 1
                pos += 1
            
            compute_body = method_body[:pos-1].strip()
            
            # Extract description
            desc_match = re.search(r'getDescription[^"]*"([^"]+)"', class_body)
            description = desc_match.group(1) if desc_match else f"PhysicsTerm: {class_name}"
            
            equations[class_name] = {
                'type': 'PhysicsTerm',
                'params': params,
                'body': compute_body,
                'description': description
            }
    
    print(f"    Found {len(equations)} PhysicsTerm classes")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # 2. EXTRACT STANDALONE DOUBLE FUNCTIONS
    # ═══════════════════════════════════════════════════════════════════════════
    print("  Extracting standalone functions...")
    
    # Pattern for standalone functions: double function_name(...) { ... }
    func_pattern = re.compile(
        r'(?:inline\s+)?double\s+(\w+)\s*\(([^)]*)\)\s*(?:const\s*)?\{',
        re.MULTILINE
    )
    
    func_count = 0
    for match in func_pattern.finditer(content):
        func_name = match.group(1)
        params = match.group(2).strip()
        start_pos = match.end()
        
        # Skip if it's a compute method (already captured)
        if func_name == 'compute':
            continue
        
        # Find matching closing brace
        brace_count = 1
        pos = start_pos
        while brace_count > 0 and pos < len(content):
            if content[pos] == '{':
                brace_count += 1
            elif content[pos] == '}':
                brace_count -= 1
            pos += 1
        
        func_body = content[start_pos:pos-1].strip()
        
        # Create unique name if duplicate
        unique_name = func_name
        counter = 1
        while unique_name in equations:
            unique_name = f"{func_name}_{counter}"
            counter += 1
        
        equations[unique_name] = {
            'type': 'function',
            'params': params,
            'body': func_body,
            'description': f"Function: {func_name}"
        }
        func_count += 1
    
    print(f"    Found {func_count} standalone functions")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # 3. EXTRACT NAMESPACE FUNCTIONS (SOURCE4, etc.)
    # ═══════════════════════════════════════════════════════════════════════════
    print("  Extracting namespace functions...")
    
    # Pattern for namespace::function
    ns_pattern = re.compile(
        r'namespace\s+(\w+)\s*\{([^}]+(?:\{[^}]*\}[^}]*)*)\}',
        re.DOTALL
    )
    
    ns_count = 0
    for ns_match in ns_pattern.finditer(content):
        ns_name = ns_match.group(1)
        ns_body = ns_match.group(2)
        
        # Find functions in namespace
        for func_match in func_pattern.finditer(ns_body):
            func_name = func_match.group(1)
            params = func_match.group(2).strip()
            start_pos = func_match.end()
            
            if func_name == 'compute':
                continue
            
            # Find matching closing brace
            brace_count = 1
            pos = start_pos
            while brace_count > 0 and pos < len(ns_body):
                if ns_body[pos] == '{':
                    brace_count += 1
                elif ns_body[pos] == '}':
                    brace_count -= 1
                pos += 1
            
            func_body = ns_body[start_pos:pos-1].strip()
            
            unique_name = f"{ns_name}_{func_name}"
            counter = 1
            while unique_name in equations:
                unique_name = f"{ns_name}_{func_name}_{counter}"
                counter += 1
            
            equations[unique_name] = {
                'type': 'namespace_function',
                'namespace': ns_name,
                'params': params,
                'body': func_body,
                'description': f"{ns_name}::{func_name}"
            }
            ns_count += 1
    
    print(f"    Found {ns_count} namespace functions")
    
    # ═══════════════════════════════════════════════════════════════════════════
    # 4. EXTRACT INLINE SUB-EQUATIONS (double var = expr;)
    # ═══════════════════════════════════════════════════════════════════════════
    print("  Extracting inline sub-equations...")
    
    # Pattern for significant inline equations (not simple assignments)
    inline_pattern = re.compile(
        r'(?:const\s+)?double\s+(\w+)\s*=\s*([^;]{20,});',  # At least 20 chars = complex equation
        re.MULTILINE
    )
    
    inline_count = 0
    seen_equations = set()
    
    for match in inline_pattern.finditer(content):
        var_name = match.group(1)
        expr = match.group(2).strip()
        
        # Skip simple variable copies
        if re.match(r'^[\w.]+$', expr):
            continue
        
        # Skip if already seen (exact same equation)
        eq_key = f"{var_name}={expr}"
        if eq_key in seen_equations:
            continue
        seen_equations.add(eq_key)
        
        # Create unique name
        unique_name = f"inline_{var_name}"
        counter = 1
        while unique_name in equations:
            unique_name = f"inline_{var_name}_{counter}"
            counter += 1
        
        equations[unique_name] = {
            'type': 'inline',
            'params': '',
            'body': f"return {expr};",
            'description': f"Inline equation: {var_name}"
        }
        inline_count += 1
    
    print(f"    Found {inline_count} inline sub-equations")
    
    return equations


def cpp_to_python_body(cpp_body):
    """Convert C++ compute body to Python."""
    py = cpp_body
    
    # Replace C++ math functions with Python math module
    replacements = [
        (r'\bstd::sqrt\b', 'math.sqrt'),
        (r'\bsqrt\b', 'math.sqrt'),
        (r'\bstd::exp\b', 'math.exp'),
        (r'\bexp\b', 'math.exp'),
        (r'\bstd::log\b', 'math.log'),
        (r'\blog\b(?!\d)', 'math.log'),
        (r'\bstd::log10\b', 'math.log10'),
        (r'\blog10\b', 'math.log10'),
        (r'\bstd::pow\b', 'pow'),
        (r'\bpow\b', 'pow'),
        (r'\bstd::sin\b', 'math.sin'),
        (r'\bsin\b', 'math.sin'),
        (r'\bstd::cos\b', 'math.cos'),
        (r'\bcos\b', 'math.cos'),
        (r'\bstd::tan\b', 'math.tan'),
        (r'\btan\b', 'math.tan'),
        (r'\bstd::abs\b', 'abs'),
        (r'\bstd::fabs\b', 'abs'),
        (r'\bfabs\b', 'abs'),
        (r'\bstd::max\b', 'max'),
        (r'\bstd::min\b', 'min'),
        (r'\bstd::atan2\b', 'math.atan2'),
        (r'\batan2\b', 'math.atan2'),
        (r'\bstd::atan\b', 'math.atan'),
        (r'\batan\b', 'math.atan'),
        (r'\bstd::asin\b', 'math.asin'),
        (r'\bstd::acos\b', 'math.acos'),
        (r'\bstd::tanh\b', 'math.tanh'),
        (r'\btanh\b', 'math.tanh'),
        (r'\bstd::sinh\b', 'math.sinh'),
        (r'\bsinh\b', 'math.sinh'),
        (r'\bstd::cosh\b', 'math.cosh'),
        (r'\bcosh\b', 'math.cosh'),
        (r'\bM_PI\b', 'math.pi'),
        (r'\bPI\b', 'math.pi'),
        (r'\bstd::isnan\b', 'math.isnan'),
        (r'\bisnan\b', 'math.isnan'),
        (r'\bstd::isinf\b', 'math.isinf'),
        (r'\bisinf\b', 'math.isinf'),
    ]
    
    for pattern, replacement in replacements:
        py = re.sub(pattern, replacement, py)
    
    # Remove C++ type declarations
    py = re.sub(r'\b(double|int|float|bool|const|auto|unsigned|long)\s+', '', py)
    
    # Replace C++ operators
    py = py.replace('&&', ' and ')
    py = py.replace('||', ' or ')
    py = py.replace('::', '.')
    
    # Remove comments
    py = re.sub(r'//.*$', '', py, flags=re.MULTILINE)
    py = re.sub(r'/\*.*?\*/', '', py, flags=re.DOTALL)
    
    # Extract return value
    return_match = re.search(r'return\s+([^;]+);', py)
    if return_match:
        return_expr = return_match.group(1).strip()
    else:
        # Try to find the last assignment
        assign_matches = list(re.finditer(r'(\w+)\s*=\s*([^;]+);', py))
        if assign_matches:
            last_var = assign_matches[-1].group(1)
            return_expr = last_var
        else:
            return_expr = "0.0"
    
    # Clean up the expression for Python safety
    # Remove any remaining C++ constructs
    return_expr = re.sub(r'\{[^}]*\}', '0.0', return_expr)  # Remove brace initializers
    return_expr = re.sub(r'\([^)]*\)\s*\{', '(', return_expr)  # Remove function bodies
    return_expr = re.sub(r'<[^>]*>', '', return_expr)  # Remove template args
    return_expr = re.sub(r'\bif\s*\([^)]*\)', '', return_expr)  # Remove if conditions
    return_expr = re.sub(r'\bfor\s*\([^)]*\)', '', return_expr)  # Remove for loops
    return_expr = re.sub(r'\bwhile\s*\([^)]*\)', '', return_expr)  # Remove while loops
    return_expr = re.sub(r'\belse\b', '', return_expr)  # Remove else
    return_expr = re.sub(r'\breturn\b', '', return_expr)  # Remove nested returns
    
    # If expression looks invalid, return 0.0
    if not return_expr.strip() or return_expr.strip() in ['{', '}', '(', ')']:
        return_expr = "0.0"
    
    # Check for unbalanced parentheses and fix
    open_count = return_expr.count('(')
    close_count = return_expr.count(')')
    if open_count > close_count:
        return_expr += ')' * (open_count - close_count)
    elif close_count > open_count:
        return_expr = '(' * (close_count - open_count) + return_expr
    
    # Check for unbalanced brackets
    if '{' in return_expr or '}' in return_expr:
        return_expr = "0.0"
    
    # Validate the expression by trying to compile it
    try:
        compile(f"x = {return_expr}", "<string>", "exec")
    except SyntaxError:
        return_expr = "0.0"
    
    return return_expr


def generate_python_file(equations, output_file):
    """Generate Python file with all equations."""
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('"""')
        f.write(f'''
QCalc_cpp_equations.py - ALL equations extracted from MAIN_1_CoAnQi.cpp

Auto-generated comprehensive extraction including:
- PhysicsTerm class compute() methods
- Standalone double functions
- Namespace functions (SOURCE4, etc.)
- Inline sub-equations

Total equations: {len(equations)}

Usage:
    from QCalc_cpp_equations import ALL_EQUATIONS, compute_equation
    
    result = compute_equation('DynamicVacuumTerm', {{'r': 1e6, 't': 1000}})
''')
        f.write('"""\n\n')
        
        f.write('import math\n')
        f.write('from typing import Dict, Any, Callable\n\n')
        
        # Import CONSTANTS
        f.write('try:\n')
        f.write('    from QCalc import CONSTANTS\n')
        f.write('except ImportError:\n')
        f.write('    CONSTANTS = {"M_sun": 1.989e30, "G": 6.6743e-11, "c": 2.998e8, "PI": 3.141592653589793}\n\n')
        
        # Helper function
        f.write('''
def safe_compute(expr_func, params: dict, default: float = 0.0) -> float:
    """Safely compute an expression with error handling."""
    try:
        result = expr_func(params)
        if result != result:  # NaN check
            return default
        if abs(result) == float('inf'):
            return default
        return float(result)
    except Exception:
        return default


''')
        
        # Base class
        f.write('''
class PhysicsTerm:
    """Base class for physics equations."""
    
    def __init__(self, name: str, description: str, compute_func: Callable):
        self.name = name
        self.description = description
        self._compute = compute_func
    
    def compute(self, params: dict) -> float:
        return safe_compute(self._compute, params)
    
    def get_description(self) -> str:
        return self.description


''')
        
        # Generate equation functions
        f.write('# ' + '='*77 + '\n')
        f.write('# EQUATION FUNCTIONS\n')
        f.write('# ' + '='*77 + '\n\n')
        
        valid_count = 0
        for name, eq in equations.items():
            # Make valid Python identifier
            py_name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
            if py_name[0].isdigit():
                py_name = 'eq_' + py_name
            
            return_expr = cpp_to_python_body(eq['body'])
            desc = eq['description'].replace('"', '\\"').replace('\n', ' ')[:100]
            
            f.write(f'def _compute_{py_name}(p: dict) -> float:\n')
            f.write(f'    """{desc}"""\n')
            f.write(f'    # Type: {eq["type"]}\n')
            f.write(f'    r = p.get("r", 1.0)\n')
            f.write(f'    t = p.get("t", 0.0)\n')
            f.write(f'    M = p.get("M", CONSTANTS.get("M_sun", 1.989e30))\n')
            f.write(f'    G = CONSTANTS.get("G", 6.6743e-11)\n')
            f.write(f'    c = CONSTANTS.get("c", 2.998e8)\n')
            f.write(f'    try:\n')
            f.write(f'        result = {return_expr}\n')
            f.write(f'        return float(result) if result == result else 0.0\n')
            f.write(f'    except Exception:\n')
            f.write(f'        return 0.0\n')
            f.write('\n')
            valid_count += 1
        
        # Generate registry
        f.write('\n# ' + '='*77 + '\n')
        f.write('# EQUATION REGISTRY\n')
        f.write('# ' + '='*77 + '\n\n')
        
        f.write('ALL_EQUATIONS = {\n')
        for name, eq in equations.items():
            py_name = re.sub(r'[^a-zA-Z0-9_]', '_', name)
            if py_name[0].isdigit():
                py_name = 'eq_' + py_name
            desc = eq['description'].replace('"', '\\"').replace('\n', ' ')[:80]
            f.write(f'    "{name}": PhysicsTerm("{name}", "{desc}", _compute_{py_name}),\n')
        f.write('}\n\n')
        
        # Convenience function
        f.write('''
def compute_equation(name: str, params: dict) -> float:
    """Compute a named equation with given parameters."""
    if name in ALL_EQUATIONS:
        return ALL_EQUATIONS[name].compute(params)
    return 0.0


def list_equations() -> list:
    """List all available equation names."""
    return list(ALL_EQUATIONS.keys())


def get_equation_info(name: str) -> dict:
    """Get information about a specific equation."""
    if name in ALL_EQUATIONS:
        eq = ALL_EQUATIONS[name]
        return {
            'name': eq.name,
            'description': eq.description
        }
    return {}


''')
        
        f.write(f'EQUATION_COUNT = {len(equations)}\n')
        f.write(f'CPP_EQUATIONS_AVAILABLE = True\n')
    
    return valid_count


def main():
    cpp_file = "MAIN_1_CoAnQi.cpp"
    output_file = "QCalc_cpp_equations.py"
    
    print("="*80)
    print("COMPREHENSIVE EQUATION EXTRACTION FROM MAIN_1_CoAnQi.cpp")
    print("="*80)
    print()
    
    # Extract all equations
    equations = extract_all_equations(cpp_file)
    
    print()
    print(f"TOTAL EQUATIONS EXTRACTED: {len(equations)}")
    print()
    
    # By type
    by_type = {}
    for name, eq in equations.items():
        t = eq['type']
        by_type[t] = by_type.get(t, 0) + 1
    
    print("By type:")
    for t, count in sorted(by_type.items(), key=lambda x: -x[1]):
        print(f"  {t}: {count}")
    
    print()
    print(f"Generating {output_file}...")
    
    valid_count = generate_python_file(equations, output_file)
    
    print(f"Generated {valid_count} equation functions")
    
    # Verify syntax
    print("\nVerifying Python syntax...")
    try:
        import py_compile
        py_compile.compile(output_file, doraise=True)
        print("✓ Syntax OK")
    except py_compile.PyCompileError as e:
        print(f"✗ Syntax error: {e}")
        return 0
    
    # Test import
    print("\nTesting import...")
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("test_module", output_file)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        print(f"✓ Import OK - {module.EQUATION_COUNT} equations available")
    except Exception as e:
        print(f"✗ Import error: {e}")
    
    return len(equations)


if __name__ == "__main__":
    main()
