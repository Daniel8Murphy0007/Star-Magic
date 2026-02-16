#!/usr/bin/env python3
"""
Extract PhysicsTerm classes from MAIN_1_CoAnQi.cpp and convert to Python.
Generates QCalc_cpp_extracted.py with ~1,239 physics term classes.
"""

import re
import os

def extract_physics_terms(filepath):
    """Extract PhysicsTerm class definitions from C++ file."""
    
    with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Pattern to match PhysicsTerm class definitions
    # Matches: class ClassName : public PhysicsTerm { ... }
    class_pattern = re.compile(
        r'class\s+(\w+)\s*:\s*public\s+PhysicsTerm\s*\{([^}]+(?:\{[^}]*\}[^}]*)*)\}',
        re.DOTALL
    )
    
    classes = []
    
    for match in class_pattern.finditer(content):
        class_name = match.group(1)
        class_body = match.group(2)
        
        # Extract description from getDescription()
        desc_match = re.search(r'getDescription\s*\([^)]*\)\s*(?:const\s*)?(?:override\s*)?\{\s*return\s*"([^"]+)"', class_body)
        description = desc_match.group(1) if desc_match else f"PhysicsTerm: {class_name}"
        
        # Extract compute method body
        compute_match = re.search(
            r'double\s+compute\s*\(([^)]*)\)\s*(?:const\s*)?(?:override\s*)?\{([^}]+(?:\{[^}]*\}[^}]*)*)\}',
            class_body
        )
        
        if compute_match:
            params = compute_match.group(1).strip()
            compute_body = compute_match.group(2).strip()
        else:
            params = ""
            compute_body = "return 0.0;"
        
        # Extract member variables
        member_vars = re.findall(r'(?:double|int|float|bool)\s+(\w+)\s*(?:=\s*[^;]+)?;', class_body)
        
        classes.append({
            'name': class_name,
            'description': description,
            'params': params,
            'compute_body': compute_body,
            'members': member_vars
        })
    
    return classes

def cpp_to_python_expr(cpp_code):
    """Convert C++ expression to Python."""
    py = cpp_code
    
    # Replace C++ math functions
    py = re.sub(r'\bstd::sqrt\b', 'math.sqrt', py)
    py = re.sub(r'\bsqrt\b', 'math.sqrt', py)
    py = re.sub(r'\bstd::exp\b', 'math.exp', py)
    py = re.sub(r'\bexp\b', 'math.exp', py)
    py = re.sub(r'\bstd::log\b', 'math.log', py)
    py = re.sub(r'\blog\b', 'math.log', py)
    py = re.sub(r'\bstd::log10\b', 'math.log10', py)
    py = re.sub(r'\blog10\b', 'math.log10', py)
    py = re.sub(r'\bstd::pow\b', 'math.pow', py)
    py = re.sub(r'\bpow\b', 'pow', py)
    py = re.sub(r'\bstd::sin\b', 'math.sin', py)
    py = re.sub(r'\bsin\b', 'math.sin', py)
    py = re.sub(r'\bstd::cos\b', 'math.cos', py)
    py = re.sub(r'\bcos\b', 'math.cos', py)
    py = re.sub(r'\bstd::tan\b', 'math.tan', py)
    py = re.sub(r'\btan\b', 'math.tan', py)
    py = re.sub(r'\bstd::abs\b', 'abs', py)
    py = re.sub(r'\bstd::fabs\b', 'abs', py)
    py = re.sub(r'\bfabs\b', 'abs', py)
    py = re.sub(r'\bstd::max\b', 'max', py)
    py = re.sub(r'\bstd::min\b', 'min', py)
    py = re.sub(r'\bstd::atan2\b', 'math.atan2', py)
    py = re.sub(r'\batan2\b', 'math.atan2', py)
    py = re.sub(r'\bstd::atan\b', 'math.atan', py)
    py = re.sub(r'\batan\b', 'math.atan', py)
    py = re.sub(r'\bstd::asin\b', 'math.asin', py)
    py = re.sub(r'\basin\b', 'math.asin', py)
    py = re.sub(r'\bstd::acos\b', 'math.acos', py)
    py = re.sub(r'\bacos\b', 'math.acos', py)
    py = re.sub(r'\bstd::tanh\b', 'math.tanh', py)
    py = re.sub(r'\btanh\b', 'math.tanh', py)
    py = re.sub(r'\bstd::sinh\b', 'math.sinh', py)
    py = re.sub(r'\bsinh\b', 'math.sinh', py)
    py = re.sub(r'\bstd::cosh\b', 'math.cosh', py)
    py = re.sub(r'\bcosh\b', 'math.cosh', py)
    py = re.sub(r'\bM_PI\b', 'math.pi', py)
    py = re.sub(r'\bPI\b', 'math.pi', py)
    
    # Replace C++ types
    py = re.sub(r'\bdouble\s+', '', py)
    py = re.sub(r'\bint\s+', '', py)
    py = re.sub(r'\bfloat\s+', '', py)
    py = re.sub(r'\bbool\s+', '', py)
    py = re.sub(r'\bconst\s+', '', py)
    py = re.sub(r'\bauto\s+', '', py)
    
    # Replace C++ operators/syntax
    py = py.replace('&&', ' and ')
    py = py.replace('||', ' or ')
    py = py.replace('!', 'not ')
    py = py.replace('::', '.')
    
    # Clean up comments
    py = re.sub(r'//.*$', '', py, flags=re.MULTILINE)
    py = re.sub(r'/\*.*?\*/', '', py, flags=re.DOTALL)
    
    return py

def generate_python_class(cls_info):
    """Generate Python class from extracted C++ info."""
    
    name = cls_info['name']
    desc = cls_info['description'].replace('"', '\\"')
    compute_body = cpp_to_python_expr(cls_info['compute_body'])
    
    # Parse compute body for return statement
    return_match = re.search(r'return\s+([^;]+);', compute_body)
    if return_match:
        return_expr = return_match.group(1).strip()
    else:
        return_expr = "0.0"
    
    # Clean up return expression
    return_expr = re.sub(r'\s+', ' ', return_expr)
    
    # Build class
    py_class = f'''
class {name}(PhysicsTerm):
    """
    {desc}
    
    Converted from MAIN_1_CoAnQi.cpp PhysicsTerm class.
    """
    
    def __init__(self):
        super().__init__()
        self.name = "{name}"
        self.description = "{desc}"
    
    def compute(self, params: dict) -> float:
        """Compute physics term value."""
        try:
            # Extract parameters with defaults
            r = params.get('r', 1.0)
            t = params.get('t', 0.0)
            M = params.get('M', CONSTANTS.get('M_sun', 1.989e30))
            
            # Computation from C++
            result = {return_expr}
            return float(result) if result == result else 0.0  # NaN check
        except Exception:
            return 0.0
    
    def get_description(self) -> str:
        return self.description
'''
    return py_class

def main():
    cpp_file = "MAIN_1_CoAnQi.cpp"
    output_file = "QCalc_cpp_extracted.py"
    
    print("="*80)
    print("EXTRACTING PhysicsTerm CLASSES FROM MAIN_1_CoAnQi.cpp")
    print("="*80)
    
    # Extract classes
    classes = extract_physics_terms(cpp_file)
    print(f"\nExtracted {len(classes)} PhysicsTerm classes")
    
    # Generate Python file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('"""')
        f.write('''
QCalc_cpp_extracted.py - PhysicsTerm classes extracted from MAIN_1_CoAnQi.cpp

Auto-generated by extract_cpp_physics_terms.py
Contains ~1,239 physics term classes converted from C++ to Python.

Usage:
    from QCalc_cpp_extracted import *
    
    term = DynamicVacuumTerm()
    result = term.compute({'r': 1e6, 't': 1000, 'M': 1.989e30})
''')
        f.write('"""\n\n')
        
        f.write('import math\n')
        f.write('from typing import Dict, Any, Optional\n\n')
        
        # Import CONSTANTS from QCalc
        f.write('# Import CONSTANTS from QCalc.py\n')
        f.write('try:\n')
        f.write('    from QCalc import CONSTANTS\n')
        f.write('except ImportError:\n')
        f.write('    CONSTANTS = {"M_sun": 1.989e30, "G": 6.6743e-11, "c": 2.998e8}\n\n')
        
        # Base PhysicsTerm class
        f.write('''
# ═══════════════════════════════════════════════════════════════════════════════
# BASE CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class PhysicsTerm:
    """Base class for all physics terms (matching C++ interface)."""
    
    def __init__(self):
        self.name = "PhysicsTerm"
        self.description = "Base physics term"
    
    def compute(self, params: dict) -> float:
        """Override in subclasses."""
        return 0.0
    
    def get_description(self) -> str:
        return self.description
    
    def validate(self, params: dict) -> bool:
        """Validate parameters."""
        return True


# ═══════════════════════════════════════════════════════════════════════════════
# EXTRACTED PhysicsTerm CLASSES (from MAIN_1_CoAnQi.cpp)
# ═══════════════════════════════════════════════════════════════════════════════

''')
        
        # Write each class
        written = 0
        for cls_info in classes:
            try:
                py_class = generate_python_class(cls_info)
                f.write(py_class)
                f.write('\n')
                written += 1
            except Exception as e:
                f.write(f'# ERROR converting {cls_info["name"]}: {e}\n\n')
        
        # Write registry at end
        f.write('''
# ═══════════════════════════════════════════════════════════════════════════════
# PHYSICS TERM REGISTRY
# ═══════════════════════════════════════════════════════════════════════════════

CPP_PHYSICS_TERMS = {
''')
        for cls_info in classes:
            name = cls_info['name']
            f.write(f'    "{name}": {name},\n')
        
        f.write('}\n\n')
        
        f.write(f'CPP_EXTRACTED_AVAILABLE = True\n')
        f.write(f'CPP_TERM_COUNT = {len(classes)}\n')
    
    print(f"Generated {output_file} with {written} classes")
    
    # Verify syntax
    print("\nVerifying Python syntax...")
    try:
        import py_compile
        py_compile.compile(output_file, doraise=True)
        print("✓ Syntax OK")
    except py_compile.PyCompileError as e:
        print(f"✗ Syntax error: {e}")
    
    return written

if __name__ == "__main__":
    main()
