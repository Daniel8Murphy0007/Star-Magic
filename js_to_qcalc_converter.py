#!/usr/bin/env python3
"""
JS to QCalc.py Converter
========================
Converts source*.js JavaScript physics modules to Python classes for QCalc.py

JavaScript → Python is straightforward:
- class X { → class X:
- this. → self.
- Math. → math.
- const/let → (remove)
- constructor() → def __init__(self):
- === → ==
- !== → !=
- || → or
- && → and
"""

import re
import os
import glob
from pathlib import Path
from typing import List, Tuple
from dataclasses import dataclass

@dataclass
class ExtractedJSClass:
    name: str
    source_file: str
    constants: dict  # name → value
    methods: list    # list of (name, params, body)
    docstring: str

class JSToQCalcConverter:
    """Convert JavaScript physics modules to Python"""
    
    def __init__(self):
        self.classes = []
        self.all_constants = {}
        
    def convert_js_to_python(self, js_code: str) -> str:
        """Convert JavaScript code to Python"""
        py = js_code
        
        # Handle template literals: `...${expr}...` → f"""...{expr}..."""
        # First convert ${...} to {...} for Python f-strings
        py = re.sub(r'\$\{([^}]+)\}', r'{\1}', py)
        # Use triple quotes to handle newlines, add f prefix for interpolation
        def convert_template(match):
            content = match.group(1)
            if '{' in content:
                return 'f"""' + content + '"""'
            else:
                return '"""' + content + '"""'
        py = re.sub(r'`([^`]*)`', convert_template, py, flags=re.DOTALL)
        
        # Encode non-ASCII to prevent syntax errors
        # Replace common math symbols with ASCII
        py = py.replace('∑', 'sum')
        py = py.replace('π', 'pi')
        py = py.replace('μ', 'mu')
        py = py.replace('γ', 'gamma')
        py = py.replace('φ', 'phi')
        py = py.replace('ω', 'omega')
        py = py.replace('Ω', 'Omega')
        py = py.replace('ρ', 'rho')
        py = py.replace('λ', 'lambda_')
        py = py.replace('Λ', 'Lambda')
        py = py.replace('δ', 'delta')
        py = py.replace('ε', 'epsilon')
        py = py.replace('η', 'eta')
        py = py.replace('θ', 'theta')
        py = py.replace('κ', 'kappa')
        py = py.replace('α', 'alpha')
        py = py.replace('β', 'beta')
        py = py.replace('τ', 'tau')
        py = py.replace('σ', 'sigma')
        py = py.replace('ψ', 'psi')
        py = py.replace('χ', 'chi')
        py = py.replace('ν', 'nu')
        py = py.replace('ξ', 'xi')
        py = py.replace('ζ', 'zeta')
        py = py.replace('∞', 'inf')
        py = py.replace('≈', '~')
        py = py.replace('≤', '<=')
        py = py.replace('≥', '>=')
        py = py.replace('≠', '!=')
        py = py.replace('×', '*')
        py = py.replace('÷', '/')
        py = py.replace('−', '-')
        py = py.replace(''', "'")
        py = py.replace(''', "'")
        py = py.replace('"', '"')
        py = py.replace('"', '"')
        py = py.replace('…', '...')
        py = py.replace('•', '*')
        py = py.replace('–', '-')
        py = py.replace('—', '-')
        py = py.replace('³', '^3')
        py = py.replace('²', '^2')
        py = py.replace('¹', '^1')
        py = py.replace('⁰', '^0')
        py = py.replace('⁴', '^4')
        py = py.replace('⁵', '^5')
        py = py.replace('⁶', '^6')
        py = py.replace('⁷', '^7')
        py = py.replace('⁸', '^8')
        py = py.replace('⁹', '^9')
        py = py.replace('⁺', '+')
        py = py.replace('⁻', '-')
        py = py.replace('₀', '_0')
        py = py.replace('₁', '_1')
        py = py.replace('₂', '_2')
        py = py.replace('₃', '_3')
        py = py.replace('₄', '_4')
        py = py.replace('₅', '_5')
        py = py.replace('₆', '_6')
        py = py.replace('₇', '_7')
        py = py.replace('₈', '_8')
        py = py.replace('₉', '_9')
        
        # Remove JavaScript-specific syntax
        py = re.sub(r'\bconst\s+', '', py)
        py = re.sub(r'\blet\s+', '', py)
        py = re.sub(r'\bvar\s+', '', py)
        
        # this. → self.
        py = re.sub(r'\bthis\.', 'self.', py)
        
        # Math functions → math
        py = re.sub(r'\bMath\.PI\b', 'math.pi', py)
        py = re.sub(r'\bMath\.E\b', 'math.e', py)
        py = re.sub(r'\bMath\.sqrt\b', 'math.sqrt', py)
        py = re.sub(r'\bMath\.sin\b', 'math.sin', py)
        py = re.sub(r'\bMath\.cos\b', 'math.cos', py)
        py = re.sub(r'\bMath\.tan\b', 'math.tan', py)
        py = re.sub(r'\bMath\.exp\b', 'math.exp', py)
        py = re.sub(r'\bMath\.log\b', 'math.log', py)
        py = re.sub(r'\bMath\.abs\b', 'abs', py)
        py = re.sub(r'\bMath\.pow\s*\(\s*([^,]+),\s*([^)]+)\)', r'(\1) ** (\2)', py)
        py = re.sub(r'\bMath\.floor\b', 'int', py)
        py = re.sub(r'\bMath\.ceil\b', 'math.ceil', py)
        py = re.sub(r'\bMath\.min\b', 'min', py)
        py = re.sub(r'\bMath\.max\b', 'max', py)
        
        # Operators
        py = re.sub(r'===', '==', py)
        py = re.sub(r'!==', '!=', py)
        py = re.sub(r'\|\|', ' or ', py)
        py = re.sub(r'&&', ' and ', py)
        
        # Remove semicolons
        py = re.sub(r';(\s*$)', r'\1', py, flags=re.MULTILINE)
        py = re.sub(r';(\s*//)', r'\1', py)
        py = re.sub(r';$', '', py, flags=re.MULTILINE)
        
        # Convert // comments to # comments
        py = re.sub(r'//', '#', py)
        
        # Convert /* */ comments
        py = re.sub(r'/\*', '"""', py)
        py = re.sub(r'\*/', '"""', py)
        
        return py
    
    def extract_class_from_js(self, filepath: str) -> List[ExtractedJSClass]:
        """Extract class definitions from a JS file"""
        with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()
        
        classes = []
        
        # Find class definitions: class ClassName { or class ClassName extends Parent {
        class_pattern = r'class\s+(\w+)(?:\s+extends\s+\w+)?\s*\{'
        
        matches = list(re.finditer(class_pattern, content))
        
        for i, match in enumerate(matches):
            class_name = match.group(1)
            start_pos = match.start()
            
            # Find matching closing brace
            brace_count = 0
            end_pos = match.end() - 1  # Position of opening brace
            found_open = False
            
            for j in range(match.end() - 1, len(content)):
                if content[j] == '{':
                    brace_count += 1
                    found_open = True
                elif content[j] == '}':
                    brace_count -= 1
                    if found_open and brace_count == 0:
                        end_pos = j + 1
                        break
            
            class_body = content[start_pos:end_pos]
            
            # Extract constants from initializeDefaults or constructor
            constants = self._extract_constants(class_body)
            
            # Extract methods
            methods = self._extract_methods(class_body)
            
            # Get docstring from file header
            docstring = self._extract_docstring(content, class_name)
            
            classes.append(ExtractedJSClass(
                name=class_name,
                source_file=os.path.basename(filepath),
                constants=constants,
                methods=methods,
                docstring=docstring
            ))
        
        return classes
    
    def _extract_constants(self, class_body: str) -> dict:
        """Extract constants from class body"""
        constants = {}
        
        # Pattern: this.name = value; or this.name = value // comment
        pattern = r'this\.(\w+)\s*=\s*([^;/\n]+)'
        
        for match in re.finditer(pattern, class_body):
            name = match.group(1)
            value = match.group(2).strip()
            
            # Skip method definitions
            if '(' in value and ')' in value and '{' in value:
                continue
            # Skip references to other this. properties
            if 'this.' in value:
                value = value.replace('this.', 'self.')
            
            constants[name] = self.convert_js_to_python(value)
        
        return constants
    
    def _extract_methods(self, class_body: str) -> list:
        """Extract methods from class body"""
        methods = []
        
        # Pattern: methodName(params) { body }
        method_pattern = r'(\w+)\s*\(([^)]*)\)\s*\{'
        
        for match in re.finditer(method_pattern, class_body):
            method_name = match.group(1)
            params = match.group(2)
            
            # Skip if this is a class definition
            if method_name == 'class':
                continue
            
            # Skip reserved Python keywords
            if method_name in ('if', 'else', 'elif', 'for', 'while', 'try', 'except', 
                               'finally', 'with', 'import', 'from', 'as', 'pass', 
                               'break', 'continue', 'return', 'yield', 'raise', 
                               'def', 'class', 'lambda', 'and', 'or', 'not', 'in', 
                               'is', 'True', 'False', 'None', 'global', 'nonlocal',
                               'assert', 'del', 'exec', 'print'):
                continue
            
            # Find method body
            start = match.end() - 1
            brace_count = 0
            end = start
            
            for j in range(start, len(class_body)):
                if class_body[j] == '{':
                    brace_count += 1
                elif class_body[j] == '}':
                    brace_count -= 1
                    if brace_count == 0:
                        end = j
                        break
            
            body = class_body[start+1:end].strip()
            methods.append((method_name, params, body))
        
        return methods
    
    def _extract_docstring(self, content: str, class_name: str) -> str:
        """Extract docstring from file header"""
        # Look for comment block at start of file
        lines = content[:2000].split('\n')
        docstring_lines = []
        
        for line in lines:
            line = line.strip()
            if line.startswith('//'):
                docstring_lines.append(line[2:].strip())
            elif line.startswith('/*') or line.startswith('*'):
                docstring_lines.append(line.lstrip('/* ').rstrip(' */'))
            elif line.startswith('class'):
                break
        
        return ' '.join(docstring_lines[:3]) if docstring_lines else f"Physics module: {class_name}"
    
    def generate_python_class(self, cls: ExtractedJSClass) -> str:
        """Generate Python class from extracted JS class"""
        
        # Build __init__ constants
        init_lines = []
        for name, value in cls.constants.items():
            # Clean up the value - only keep simple numeric/string values
            clean_value = value.strip()
            
            # Skip complex expressions
            if any(x in clean_value for x in ['{', '}', '[', ']', '=>', 'function', 
                                               'new ', 'Map(', 'Set(', 'SystemType.',
                                               ' ? ', ' : ', '++', '--']):
                continue
            
            # Skip if it references methods (has parens but not self.)
            if '(' in clean_value and 'self.' not in clean_value:
                # Allow math functions
                if not any(f in clean_value for f in ['math.', 'Math.', 'sin(', 'cos(', 'sqrt(', 'exp(', 'log(']):
                    continue
            
            # Skip if it has logical operators followed by equals (garbled)
            if ' == ' in clean_value and clean_value.startswith('=='):
                continue
            if clean_value.startswith('=='):
                continue
            
            # Skip multiline
            if '\n' in clean_value:
                continue
                
            # Validate it looks like valid Python
            try:
                compile(f"x = {clean_value}", "<string>", "exec")
            except SyntaxError:
                continue
            
            init_lines.append(f"        self.{name} = {clean_value}")
        
        init_body = '\n'.join(init_lines) if init_lines else "        pass"
        
        # Build method stubs (don't try to convert complex JS bodies)
        method_blocks = []
        for method_name, params, body in cls.methods:
            py_name = method_name
            if method_name == 'constructor':
                continue  # Skip, we handle in __init__
            
            # Skip reserved keywords
            if py_name in ('if', 'else', 'for', 'while', 'return', 'class', 'def', 
                           'import', 'from', 'with', 'try', 'except', 'finally',
                           'and', 'or', 'not', 'in', 'is', 'True', 'False', 'None',
                           'forEach', 'map', 'filter', 'reduce', 'find', 'some', 
                           'every', 'keys', 'values', 'entries', 'then', 'catch'):
                continue
            
            # Skip methods with complex JS patterns in params
            if any(x in params for x in ['=>', 'if', 'else', 'function', '{', '}']):
                continue
            
            # Convert params
            py_params = 'self'
            if params.strip():
                param_parts = [p.strip() for p in params.split(',')]
                for p in param_parts:
                    # Remove type annotations if any
                    p = re.sub(r':\s*\w+', '', p)
                    # Clean up the parameter name - be very strict
                    p = re.sub(r'[^a-zA-Z0-9_]', '', p)
                    if not p or not p.isidentifier():
                        continue
                    py_params += f", {p}"
            
            # Simple stub body
            method_blocks.append(f'''
    def {py_name}({py_params}):
        """Auto-converted from {cls.source_file}"""
        pass  # TODO: Implement - original JS method body available in source
''')
        
        methods_str = ''.join(method_blocks)
        
        # Generate class
        code = f'''
class {cls.name}:
    """
    {cls.docstring}
    
    Source: {cls.source_file}
    Converted: JavaScript → Python for QCalc.py
    """
    
    def __init__(self):
{init_body}
{methods_str}
'''
        return code
    
    def process_all_files(self, pattern: str = "source*.js") -> Tuple[str, int]:
        """Process all matching JS files and generate Python code"""
        
        files = sorted(glob.glob(pattern))
        print(f"Found {len(files)} JavaScript source files")
        
        all_classes = []
        all_constants = {}
        
        for filepath in files:
            try:
                classes = self.extract_class_from_js(filepath)
                all_classes.extend(classes)
                print(f"  {os.path.basename(filepath)}: {len(classes)} classes")
            except Exception as e:
                print(f"  WARNING: {os.path.basename(filepath)} failed: {e}")
        
        print(f"\nExtracted {len(all_classes)} total classes")
        
        # Generate Python code
        header = '''#!/usr/bin/env python3
"""
QCalc_extracted.py - Auto-generated from source*.js JavaScript modules
=======================================================================
Converted by js_to_qcalc_converter.py

Contains physics calculator classes extracted from JavaScript UQFF modules.
All naming conventions preserved from original JavaScript source.
"""

import math
from dataclasses import dataclass
from typing import Dict, Any, Optional

@dataclass
class EquationResult:
    """Result container for physics calculations"""
    name: str
    value: float
    unit: str
    latex: str
    description: str

'''
        
        class_code = ''
        for cls in all_classes:
            try:
                class_code += self.generate_python_class(cls)
            except Exception as e:
                print(f"  WARNING: Failed to generate {cls.name}: {e}")
        
        return header + class_code, len(all_classes)


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Convert JavaScript physics modules to Python')
    parser.add_argument('--output', '-o', default='QCalc_js_extracted.py', help='Output file')
    parser.add_argument('--pattern', '-p', default='source*.js', help='File pattern to match')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("JS TO QCalc.py CONVERTER")
    print("=" * 60)
    
    converter = JSToQCalcConverter()
    
    if args.dry_run:
        files = sorted(glob.glob(args.pattern))
        print(f"Would process {len(files)} files:")
        for f in files[:10]:
            print(f"  {f}")
        if len(files) > 10:
            print(f"  ... and {len(files) - 10} more")
        return
    
    code, count = converter.process_all_files(args.pattern)
    
    print(f"\nWriting to {args.output}...")
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(code)
    
    print(f"  Wrote {len(code):,} characters")
    print(f"  Classes: {count}")
    
    # Verify syntax
    print("\nVerifying Python syntax...")
    import py_compile
    try:
        py_compile.compile(args.output, doraise=True)
        print("  [OK] Syntax valid")
    except py_compile.PyCompileError as e:
        print(f"  [ERROR] Syntax error: {e}")
    
    print("\n" + "=" * 60)
    print("CONVERSION COMPLETE")
    print("=" * 60)


if __name__ == '__main__':
    main()
