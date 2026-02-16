#!/usr/bin/env python3
"""
================================================================================
C++ TO QCalc.py CONVERTER
================================================================================
Purpose: Automated extraction of PhysicsTerm classes from MAIN_1_CoAnQi.cpp
         and generation of equivalent Python calculators for QCalc.py

Source: MAIN_1_CoAnQi.cpp (106,894 lines, 1,240+ PhysicsTerm classes)
        wolfram_physics_classes.cpp (132,550 lines, 5,703 classes)
        
Target: QCalc.py (monolithic calculator architecture)

NAMING CONVENTION: ALL NAMES PRESERVED EXACTLY FROM C++ SOURCE
- Class names: DynamicVacuumTerm → DynamicVacuumTerm (unchanged)
- Constants: G, c_light, hbar → G, c_light, hbar (unchanged)
- Methods: compute() → compute() (unchanged)

================================================================================
USAGE:
    python cpp_to_qcalc_converter.py [--dry-run] [--output QCalc_generated.py]
    
    --dry-run       Preview without writing files
    --output FILE   Write to specific file (default: append to QCalc.py)
    --constants     Extract constants only
    --classes       Extract classes only
    --limit N       Limit to first N classes (for testing)
================================================================================
Author: Automated extraction for Star-Magic UQFF Framework
Date: 2026-02-15
================================================================================
"""

import re
import sys
import argparse
from pathlib import Path
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
import json

# =============================================================================
# CONFIGURATION
# =============================================================================

MAIN_CPP = "MAIN_1_CoAnQi.cpp"
WOLFRAM_CPP = "wolfram_physics_classes.cpp"
QCALC_PY = "QCalc.py"
OUTPUT_BACKUP = "QCalc_backup.py"

# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class ExtractedConstant:
    """Represents a constant extracted from C++"""
    name: str
    value: str
    cpp_type: str  # 'double', 'int', 'constexpr'
    comment: str
    source_file: str
    line_number: int

@dataclass
class ExtractedMethod:
    """Represents a method body extracted from C++"""
    name: str
    return_type: str
    parameters: List[str]
    body: str
    line_start: int
    line_end: int

@dataclass
class ExtractedClass:
    """Represents a PhysicsTerm class extracted from C++"""
    name: str
    base_class: str
    constants: List[ExtractedConstant]
    methods: List[ExtractedMethod]
    compute_body: str
    get_name_value: str
    get_description_value: str
    source_file: str
    line_start: int
    line_end: int
    metadata: Dict[str, str]

# =============================================================================
# C++ PARSER
# =============================================================================

class CppPhysicsTermParser:
    """
    Parses C++ PhysicsTerm classes and extracts:
    - Class name (preserved exactly)
    - Member constants
    - compute() method body
    - getName() return value
    - getDescription() return value
    """
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.content = ""
        self.lines = []
        self.classes: List[ExtractedClass] = []
        self.constants: List[ExtractedConstant] = []
        
    def load(self):
        """Load and preprocess C++ file"""
        print(f"Loading {self.filepath}...")
        with open(self.filepath, 'r', encoding='utf-8', errors='replace') as f:
            self.content = f.read()
        self.lines = self.content.split('\n')
        print(f"  Loaded {len(self.lines):,} lines")
        
    def extract_global_constants(self) -> List[ExtractedConstant]:
        """Extract global const/constexpr definitions"""
        patterns = [
            # const double NAME = VALUE; // comment
            r'^const\s+double\s+(\w+)\s*=\s*([^;]+);?\s*(?://\s*(.*))?$',
            # const int NAME = VALUE;
            r'^const\s+int\s+(\w+)\s*=\s*([^;]+);?\s*(?://\s*(.*))?$',
            # constexpr double NAME = VALUE;
            r'^constexpr\s+double\s+(\w+)\s*=\s*([^;]+);?\s*(?://\s*(.*))?$',
            # constexpr int NAME = VALUE;
            r'^constexpr\s+int\s+(\w+)\s*=\s*([^;]+);?\s*(?://\s*(.*))?$',
            # #define NAME VALUE
            r'^#define\s+(\w+)\s+([^\s/]+)\s*(?://\s*(.*))?$',
        ]
        
        constants = []
        for i, line in enumerate(self.lines):
            line = line.strip()
            for pattern in patterns:
                match = re.match(pattern, line)
                if match:
                    name = match.group(1)
                    value = match.group(2).strip()
                    comment = match.group(3) if match.lastindex >= 3 else ""
                    
                    # Determine type
                    if 'constexpr' in line:
                        cpp_type = 'constexpr'
                    elif 'const int' in line:
                        cpp_type = 'int'
                    elif 'const double' in line:
                        cpp_type = 'double'
                    elif '#define' in line:
                        cpp_type = 'define'
                    else:
                        cpp_type = 'unknown'
                    
                    constants.append(ExtractedConstant(
                        name=name,
                        value=value,
                        cpp_type=cpp_type,
                        comment=comment or "",
                        source_file=self.filepath,
                        line_number=i + 1
                    ))
                    break
        
        self.constants = constants
        print(f"  Extracted {len(constants)} global constants")
        return constants
    
    def find_class_boundaries(self) -> List[Tuple[int, int, str]]:
        """Find start/end line numbers for each PhysicsTerm class"""
        # Pattern: class ClassName : public PhysicsTerm
        class_pattern = r'^class\s+(\w+)\s*:\s*public\s+PhysicsTerm'
        
        boundaries = []
        i = 0
        while i < len(self.lines):
            match = re.match(class_pattern, self.lines[i].strip())
            if match:
                class_name = match.group(1)
                start_line = i
                
                # Find matching closing brace
                brace_count = 0
                found_open = False
                for j in range(i, len(self.lines)):
                    line = self.lines[j]
                    brace_count += line.count('{') - line.count('}')
                    if '{' in line:
                        found_open = True
                    if found_open and brace_count == 0:
                        boundaries.append((start_line, j, class_name))
                        i = j
                        break
            i += 1
        
        print(f"  Found {len(boundaries)} PhysicsTerm class boundaries")
        return boundaries
    
    def extract_class(self, start: int, end: int, name: str) -> ExtractedClass:
        """Extract a single PhysicsTerm class"""
        class_lines = self.lines[start:end+1]
        class_text = '\n'.join(class_lines)
        
        # Extract compute() body
        compute_body = self._extract_method_body(class_text, 'compute')
        
        # Extract getName() return value
        get_name_match = re.search(r'getName\s*\(\s*\)\s*const[^{]*\{\s*return\s*"([^"]+)"', class_text)
        get_name_value = get_name_match.group(1) if get_name_match else name
        
        # Extract getDescription() return value
        get_desc_match = re.search(r'getDescription\s*\(\s*\)\s*const[^{]*\{\s*return\s*"([^"]+)"', class_text)
        get_description_value = get_desc_match.group(1) if get_desc_match else f"Physics term: {name}"
        
        # Extract local constants (member variables)
        local_constants = self._extract_local_constants(class_text, start)
        
        # Extract metadata
        metadata = {}
        meta_matches = re.findall(r'setMetadata\s*\(\s*"(\w+)"\s*,\s*"([^"]+)"\s*\)', class_text)
        for key, value in meta_matches:
            metadata[key] = value
        
        return ExtractedClass(
            name=name,
            base_class="PhysicsTerm",
            constants=local_constants,
            methods=[],
            compute_body=compute_body,
            get_name_value=get_name_value,
            get_description_value=get_description_value,
            source_file=self.filepath,
            line_start=start + 1,
            line_end=end + 1,
            metadata=metadata
        )
    
    def _extract_method_body(self, class_text: str, method_name: str) -> str:
        """Extract the body of a method"""
        # Find method signature and body
        pattern = rf'{method_name}\s*\([^)]*\)\s*(?:const)?\s*(?:override)?\s*\{{'
        match = re.search(pattern, class_text)
        if not match:
            return ""
        
        start_pos = match.end() - 1  # Position of opening brace
        brace_count = 1
        end_pos = start_pos + 1
        
        while end_pos < len(class_text) and brace_count > 0:
            if class_text[end_pos] == '{':
                brace_count += 1
            elif class_text[end_pos] == '}':
                brace_count -= 1
            end_pos += 1
        
        body = class_text[start_pos+1:end_pos-1].strip()
        return body
    
    def _extract_local_constants(self, class_text: str, base_line: int) -> List[ExtractedConstant]:
        """Extract member variable constants from a class"""
        constants = []
        
        # Only extract simple member initializations, NOT constructor params
        # Pattern: in member variable section (before methods)
        # Look for: double name = value; (standalone, not in constructor)
        
        # Find the private/public section with member variables
        member_section = re.search(r'private:\s*(.*?)(?:public:|protected:|\};)', class_text, re.DOTALL)
        if not member_section:
            return constants
            
        member_text = member_section.group(1)
        
        # Only match simple initializations (not function params)
        patterns = [
            r'^\s*double\s+(\w+)\s*;\s*$',  # double name;
            r'^\s*const\s+double\s+(\w+)\s*;\s*$',  # const double name;
        ]
        
        for line in member_text.split('\n'):
            for pattern in patterns:
                match = re.match(pattern, line)
                if match:
                    constants.append(ExtractedConstant(
                        name=match.group(1),
                        value="0.0",  # Default, actual value set in constructor
                        cpp_type='double',
                        comment="",
                        source_file=self.filepath,
                        line_number=base_line
                    ))
        
        return constants
    
    def extract_all_classes(self, limit: Optional[int] = None) -> List[ExtractedClass]:
        """Extract all PhysicsTerm classes"""
        boundaries = self.find_class_boundaries()
        
        if limit:
            boundaries = boundaries[:limit]
        
        classes = []
        for start, end, name in boundaries:
            try:
                cls = self.extract_class(start, end, name)
                classes.append(cls)
            except Exception as e:
                print(f"  WARNING: Failed to extract {name}: {e}")
        
        self.classes = classes
        print(f"  Extracted {len(classes)} complete classes")
        return classes

# =============================================================================
# PYTHON CODE GENERATOR
# =============================================================================

class PythonCalculatorGenerator:
    """
    Generates Python calculator classes from extracted C++ PhysicsTerm classes
    
    NAMING: All names preserved exactly from C++ source
    OUTPUT FORMAT: Matches QCalc.py patterns (EquationResult, dataclass, etc.)
    """
    
    # Python reserved keywords that need renaming WHEN USED AS VARIABLE NAMES
    # Only include keywords that could plausibly be C++ variable names
    # Do NOT include: return, if, else, for, while, and, or, not (Python syntax)
    PYTHON_RESERVED = {
        'lambda': 'lambda_',
        'class': 'class_',
        'global': 'global_',
        'True': 'True_',
        'False': 'False_',
        'None': 'None_',
        'yield': 'yield_',
    }
    
    def __init__(self):
        self.indent = "    "
    
    def sanitize_identifier(self, name: str) -> str:
        """Rename Python reserved keywords in identifiers"""
        if name in self.PYTHON_RESERVED:
            return self.PYTHON_RESERVED[name]
        return name
        
    def cpp_value_to_python(self, value: str) -> str:
        """Convert C++ constant value to Python"""
        # Handle scientific notation
        value = value.strip()
        
        # Sanitize reserved keywords used as identifiers
        # Replace all occurrences of reserved words used as variable names
        for reserved, replacement in self.PYTHON_RESERVED.items():
            # Use word boundary to replace all variable-like usage
            value = re.sub(rf'\b{reserved}\b', replacement, value)
        
        # Remove C++ digit separators (1'000'000 -> 1000000)
        value = re.sub(r"(\d)'(\d)", r'\1\2', value)
        
        # Replace M_PI with math.pi
        value = re.sub(r'\bM_PI\b', 'math.pi', value)
        value = re.sub(r'\bPI\b', 'math.pi', value)
        
        # Replace pow(x, y) with x ** y or math.pow(x, y)
        value = re.sub(r'pow\s*\(\s*([^,]+),\s*([^)]+)\)', r'(\1) ** (\2)', value)
        
        # Replace sqrt with math.sqrt (only if not already prefixed)
        value = re.sub(r'(?<!math\.)\bsqrt\s*\(', 'math.sqrt(', value)
        
        # Replace sin, cos, exp, log - but not if already prefixed with math.
        value = re.sub(r'(?<!math\.)\bsin\s*\(', 'math.sin(', value)
        value = re.sub(r'(?<!math\.)\bcos\s*\(', 'math.cos(', value)
        value = re.sub(r'(?<!math\.)\bexp\s*\(', 'math.exp(', value)
        value = re.sub(r'(?<!math\.)\blog\s*\(', 'math.log(', value)
        value = re.sub(r'(?<!math\.)\babs\s*\(', 'abs(', value)
        
        # Convert C++ logical operators to Python
        value = re.sub(r'\|\|', ' or ', value)
        value = re.sub(r'&&', ' and ', value)
        value = re.sub(r'!(?!=)', ' not ', value)  # ! but not !=
        
        # Replace 1e-10 style (already valid in Python)
        # Replace 1.0e-10 style (already valid in Python)
        
        return value
    
    def cpp_body_to_python(self, body: str) -> str:
        """Convert C++ method body to Python"""
        if not body:
            return "result = 0.0  # TODO: Implement computation"
        
        # Pre-process: merge C++ multi-line if statements
        # Pattern: if (condition)\n  return x; → if (condition) return x;
        body = re.sub(
            r'if\s*\(([^)]+)\)\s*[\r\n]+\s*(return [^;]+;)',
            r'if (\1) \2',
            body
        )
        
        # Pre-process: merge continuation lines (lines ending with operators)
        # C++ allows: x = a + b +\n    c; → x = a + b + c;
        body = re.sub(r'([+\-*/=,])\s*[\r\n]+\s*', r'\1 ', body)
        
        # Pre-process: merge simple for loops with their single-statement body
        # for (...)\n    statement; → for (...) statement;
        body = re.sub(
            r'for\s*\([^)]+\)\s*[\r\n]+\s*([^;{]+;)',
            lambda m: f"// FOR_LOOP: {m.group(0)[:80]}",  # Convert to comment for now
            body
        )
        
        lines = body.split('\n')
        python_lines = []
        
        for line in lines:
            line = line.strip()
            if not line:
                continue
            
            # Handle full-line C++ comments
            if line.startswith('//'):
                # Convert to Python comment, stripping non-ASCII
                comment = line[2:].strip()
                comment = comment.encode('ascii', 'replace').decode('ascii')
                python_lines.append(f"# {comment}")
                continue
            if line.startswith('/*') or line.startswith('*') or line.endswith('*/'):
                continue
            
            # CRITICAL: Strip inline C++ comments before processing
            if '//' in line:
                line = line.split('//')[0].strip()
            
            # Skip constructor initialization list remnants
            if line.startswith(':') or line.startswith('setMetadata'):
                continue
            if 'double ' in line and '=' in line and ')' in line and '(' in line:
                # Constructor parameter - skip
                continue
            
            # Convert variable declarations
            # const double x = ...; → x = ...
            line = re.sub(r'^const\s+double\s+', '', line)
            line = re.sub(r'^double\s+', '', line)
            line = re.sub(r'^const\s+int\s+', '', line)
            line = re.sub(r'^int\s+', '', line)
            
            # Remove semicolons
            line = line.rstrip(';')
            
            # Skip empty lines after processing
            if not line.strip():
                continue
            
            # Skip any remaining for loops (converted to comments in preprocessing)
            if line.startswith('for ') or line.startswith('for('):
                python_lines.append(f"# TODO: Convert C++ for loop")
                continue
            
            # Convert C++ single-line if with body: if (cond) body; → if cond: body
            inline_if = re.match(r'^if\s*\((.+?)\)\s*(.+)$', line)
            if inline_if:
                condition = inline_if.group(1)
                body = inline_if.group(2).strip().rstrip(';')
                condition = self.cpp_value_to_python(condition)
                body = self.cpp_value_to_python(body)
                if body.startswith('return'):
                    line = f"if {condition}: {body}"
                else:
                    line = f"if {condition}:\n            {body}"
            else:
                # Convert C++ if statements to Python (multi-line)
                # if (condition) → if condition:
                if_match = re.match(r'^if\s*\((.+)\)\s*$', line)
                if if_match:
                    condition = if_match.group(1)
                    condition = self.cpp_value_to_python(condition)
                    line = f"if {condition}:"
            
            # Convert C++ else if → elif
            elif_match = re.match(r'^else\s+if\s*\((.+)\)\s*$', line)
            if elif_match:
                condition = elif_match.group(1)
                condition = self.cpp_value_to_python(condition)
                line = f"elif {condition}:"
            
            # Convert else
            if line.strip() == 'else':
                line = 'else:'
            
            # Convert return statements (standalone)
            return_match = re.match(r'^return\s+(.+)$', line)
            if return_match and not line.startswith('if'):
                value = return_match.group(1)
                value = self.cpp_value_to_python(value)
                line = f"return {value}"
            
            # Convert C++ ternary operator with params.count
            # params.count("key") ? params.at("key") : default
            line = re.sub(
                r'params\.count\s*\(\s*"(\w+)"\s*\)\s*\?\s*params\.at\s*\(\s*"\1"\s*\)\s*:\s*([^\s;]+)',
                r'params.get("\1", \2)',
                line
            )
            
            # Convert std::map access
            # params.at("key") → params.get("key", 0.0)
            line = re.sub(r'params\.at\s*\(\s*"(\w+)"\s*\)', r'params.get("\1", 0.0)', line)
            line = re.sub(r'params\["(\w+)"\]', r'params.get("\1", 0.0)', line)
            
            # Convert math functions
            line = self.cpp_value_to_python(line)
            
            # Convert getDynamicParameter
            line = re.sub(r'getDynamicParameter\s*\(\s*"(\w+)"\s*,\s*([^)]+)\)', 
                         r'self.get_param("\1", \2)', line)
            
            # Skip any remaining raw C++ patterns
            if '{' in line or '}' in line:
                continue
            if 'std::' in line:
                continue
            
            python_lines.append(line)
        
        if not python_lines:
            return "result = 0.0  # TODO: Implement computation"
        
        # Ensure we have a result assignment
        result_code = '\n        '.join(python_lines)
        if 'result' not in result_code and 'return' not in result_code:
            result_code += '\n        result = 0.0  # Assign result'
        
        return result_code
    
    def generate_constant(self, const: ExtractedConstant) -> str:
        """Generate Python constant entry for CONSTANTS dict"""
        py_value = self.cpp_value_to_python(const.value)
        comment = f"  # {const.comment}" if const.comment else ""
        return f"    '{const.name}': {py_value},{comment}"
    
    def generate_class(self, cls: ExtractedClass) -> str:
        """Generate Python calculator class from extracted C++ class"""
        
        # Generate local constants as class attributes
        local_consts = ""
        if cls.constants:
            const_lines = []
            for c in cls.constants:
                py_val = self.cpp_value_to_python(c.value)
                safe_name = self.sanitize_identifier(c.name)
                const_lines.append(f"        self.{safe_name} = {py_val}")
            local_consts = '\n'.join(const_lines)
        
        # Generate compute body
        compute_body = self.cpp_body_to_python(cls.compute_body)
        
        # Generate metadata
        metadata_lines = []
        for key, value in cls.metadata.items():
            # Escape single quotes in value to prevent string termination issues
            safe_value = str(value).replace("'", "\\'")
            metadata_lines.append(f"        self.metadata['{key}'] = '{safe_value}'")
        metadata_init = '\n'.join(metadata_lines) if metadata_lines else ""
        
        # Escape strings for Python
        safe_name = str(cls.get_name_value).replace('"', '\\"')
        safe_desc = str(cls.get_description_value).replace('"', '\\"')
        
        # Build class
        code = f'''
class {cls.name}:
    """
    {safe_desc}
    
    Source: {cls.source_file} (lines {cls.line_start}-{cls.line_end})
    Naming: Preserved exactly from C++ PhysicsTerm class
    """
    
    def __init__(self):
        self.name = "{safe_name}"
        self.description = "{safe_desc}"
        self.metadata = {{}}
{metadata_init}
{local_consts}
    
    def get_param(self, key: str, default: float = 0.0) -> float:
        """Get dynamic parameter (compatibility layer)"""
        return default  # TODO: Implement dynamic parameter storage
    
    def compute(self, t: float = 0.0, params: dict = None) -> EquationResult:
        """
        Compute physics term value
        
        Args:
            t: Time parameter (seconds)
            params: Dictionary of input parameters
            
        Returns:
            EquationResult with computed value and LaTeX equation
        """
        if params is None:
            params = {{}}
        
        # Computation (auto-converted from C++)
        {compute_body}
        
        return EquationResult(
            name=self.name,
            value=result if 'result' in dir() else 0.0,
            unit="SI",
            latex=r"\\text{{{cls.name}}}",
            description=self.description
        )
'''
        return code
    
    def generate_constants_block(self, constants: List[ExtractedConstant]) -> str:
        """Generate CONSTANTS dictionary additions"""
        lines = ["# ========== AUTO-EXTRACTED CONSTANTS FROM C++ =========="]
        lines.append("# Source: MAIN_1_CoAnQi.cpp")
        lines.append("# Naming: Preserved exactly from C++ source")
        lines.append("EXTRACTED_CONSTANTS = {")
        
        # Group by source/type
        for const in constants:
            lines.append(self.generate_constant(const))
        
        lines.append("}")
        lines.append("")
        lines.append("# Merge into main CONSTANTS")
        lines.append("CONSTANTS.update(EXTRACTED_CONSTANTS)")
        lines.append("")
        
        return '\n'.join(lines)
    
    def generate_all(self, constants: List[ExtractedConstant], 
                     classes: List[ExtractedClass]) -> str:
        """Generate complete Python code block"""
        
        output = []
        output.append("""
# ==============================================================================
# AUTO-GENERATED PHYSICS CALCULATORS
# ==============================================================================
# Source: MAIN_1_CoAnQi.cpp (C++ PhysicsTerm classes)
# Generator: cpp_to_qcalc_converter.py
# Date: 2026-02-15
# 
# NAMING CONVENTION: ALL NAMES PRESERVED EXACTLY FROM C++ SOURCE
# - Class names unchanged
# - Constant names unchanged  
# - Method signatures preserved
# ==============================================================================
""")
        
        # Generate constants
        if constants:
            output.append(self.generate_constants_block(constants))
        
        # Generate classes
        output.append("# ========== EXTRACTED CALCULATOR CLASSES ==========")
        for cls in classes:
            try:
                output.append(self.generate_class(cls))
            except Exception as e:
                output.append(f"# ERROR generating {cls.name}: {e}")
        
        # Generate registry
        output.append(self._generate_registry(classes))
        
        return '\n'.join(output)
    
    def _generate_registry(self, classes: List[ExtractedClass]) -> str:
        """Generate calculator registry"""
        lines = [
            "",
            "# ========== CALCULATOR REGISTRY ==========",
            "EXTRACTED_CALCULATORS = {",
        ]
        for cls in classes:
            lines.append(f"    '{cls.name}': {cls.name},")
        lines.append("}")
        lines.append("")
        return '\n'.join(lines)

# =============================================================================
# MAIN CONVERTER
# =============================================================================

def run_converter(args):
    """Main converter execution"""
    
    print("=" * 60)
    print("C++ TO QCalc.py CONVERTER")
    print("=" * 60)
    
    # Parse MAIN_1_CoAnQi.cpp
    parser = CppPhysicsTermParser(MAIN_CPP)
    parser.load()
    
    # Extract constants
    if not args.classes_only:
        print("\nExtracting global constants...")
        constants = parser.extract_global_constants()
    else:
        constants = []
    
    # Extract classes
    if not args.constants_only:
        print("\nExtracting PhysicsTerm classes...")
        classes = parser.extract_all_classes(limit=args.limit)
    else:
        classes = []
    
    # Generate Python code
    print("\nGenerating Python code...")
    generator = PythonCalculatorGenerator()
    python_code = generator.generate_all(constants, classes)
    
    # Output
    if args.dry_run:
        print("\n[DRY RUN] Generated code preview (first 2000 chars):")
        print("-" * 60)
        print(python_code[:2000])
        print("-" * 60)
        print(f"\nTotal generated: {len(python_code):,} characters")
        print(f"Constants: {len(constants)}")
        print(f"Classes: {len(classes)}")
    else:
        output_file = args.output or "QCalc_generated.py"
        print(f"\nWriting to {output_file}...")
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(python_code)
        print(f"  Wrote {len(python_code):,} characters")
        print(f"  Constants: {len(constants)}")
        print(f"  Classes: {len(classes)}")
    
    # Generate summary
    summary = {
        "source_file": MAIN_CPP,
        "source_lines": len(parser.lines),
        "constants_extracted": len(constants),
        "classes_extracted": len(classes),
        "output_file": args.output or "QCalc_generated.py",
        "output_size": len(python_code)
    }
    
    with open("converter_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n" + "=" * 60)
    print("CONVERSION COMPLETE")
    print("=" * 60)
    print(f"Summary saved to: converter_summary.json")
    
    return summary

# =============================================================================
# CLI
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert C++ PhysicsTerm classes to Python QCalc.py format"
    )
    parser.add_argument('--dry-run', action='store_true',
                       help="Preview without writing files")
    parser.add_argument('--output', '-o', type=str,
                       help="Output file (default: QCalc_generated.py)")
    parser.add_argument('--constants-only', action='store_true',
                       help="Extract constants only")
    parser.add_argument('--classes-only', action='store_true',
                       help="Extract classes only")
    parser.add_argument('--limit', '-n', type=int,
                       help="Limit number of classes to extract")
    
    args = parser.parse_args()
    run_converter(args)
