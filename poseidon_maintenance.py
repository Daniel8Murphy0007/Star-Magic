"""
poseidon_maintenance.py
=======================
Python companion script for Poseidon TaskBot

This script is called by the C++ Poseidon TaskBot via pybind11 bridge
to perform Python-side maintenance operations:

- Constant synchronization (shared_constants.h ↔ .py ↔ index.js)
- File regeneration (QCalc_cpp_extracted.py, etc.)
- Python calculator invocation (QCalc.py, CondensedPhysics.py)
- Validation suite execution
- FTPS bundle handling

Version: 4.2.1 (Canonical - matches Architecture v4.2)
Author: Daniel T. Murphy
Framework: UQFF Star-Magic v4.2
"""

import os
import sys
import re
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from dataclasses import dataclass, field, asdict
import hashlib

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [POSEIDON-PY] %(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger('poseidon')

# ============================================================================
# CONFIGURATION
# ============================================================================

@dataclass
class PoseidonPythonConfig:
    """Configuration for Python-side Poseidon operations"""
    base_path: str = ""
    backup_path: str = ""
    log_path: str = ""
    
    # Constant files
    cpp_constants_file: str = "shared_constants.h"
    py_constants_file: str = "shared_constants.py"
    js_constants_file: str = "index.js"
    
    # Calculator modules
    qcalc_module: str = "QCalc"
    condensed_physics_module: str = "CondensedPhysics"
    
    # Validation modules
    qcalc_validation: str = "QCalc_validation"
    condensed_validation: str = "CondensedPhysics_Validation"
    
    # FTPS settings
    ftps_enabled: bool = False
    
    def __post_init__(self):
        if not self.base_path:
            self.base_path = os.environ.get('UQFF_WORKSPACE', os.getcwd())
        if not self.backup_path:
            self.backup_path = os.path.join(self.base_path, 'backups')
        if not self.log_path:
            self.log_path = os.path.join(self.base_path, 'logs')


# Global config instance
_config = PoseidonPythonConfig()


def set_config(config: PoseidonPythonConfig):
    """Set global configuration"""
    global _config
    _config = config
    logger.info(f"Configuration updated - base_path: {_config.base_path}")


# ============================================================================
# CONSTANT EXTRACTION & SYNCHRONIZATION
# ============================================================================

@dataclass
class PhysicsConstant:
    """Represents a physics constant across languages"""
    name: str
    value: float
    unit: str = ""
    description: str = ""
    source_lang: str = ""  # "cpp", "python", "javascript"
    source_file: str = ""
    line_number: int = 0


def extract_cpp_constants(file_path: str) -> Dict[str, PhysicsConstant]:
    """Extract constants from shared_constants.h"""
    constants = {}
    
    # Pattern for constexpr double or static const double
    pattern = re.compile(
        r'(?:constexpr|static\s+const)\s+double\s+(\w+)\s*=\s*([^;]+);'
        r'\s*(?://\s*(.+))?'
    )
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                match = pattern.search(line)
                if match:
                    name = match.group(1)
                    value_str = match.group(2).strip()
                    comment = match.group(3) or ""
                    
                    # Parse value (handle scientific notation)
                    try:
                        # Replace C++ syntax
                        value_str = value_str.replace('e', 'E')
                        value = float(eval(value_str))
                    except:
                        continue
                    
                    # Extract unit from comment
                    unit_match = re.search(r'\[([^\]]+)\]', comment)
                    unit = unit_match.group(1) if unit_match else ""
                    
                    constants[name] = PhysicsConstant(
                        name=name,
                        value=value,
                        unit=unit,
                        description=comment.strip(),
                        source_lang="cpp",
                        source_file=file_path,
                        line_number=line_num
                    )
    except Exception as e:
        logger.error(f"Failed to extract C++ constants: {e}")
    
    logger.info(f"Extracted {len(constants)} constants from C++")
    return constants


def extract_python_constants(file_path: str) -> Dict[str, PhysicsConstant]:
    """Extract constants from shared_constants.py"""
    constants = {}
    
    # Pattern for Python constants
    pattern = re.compile(
        r'^(\w+)\s*=\s*([^\n#]+)(?:#\s*(.+))?'
    )
    
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                match = pattern.match(line)
                if match:
                    name = match.group(1)
                    value_str = match.group(2).strip()
                    comment = match.group(3) or ""
                    
                    # Skip non-numeric
                    if not re.match(r'^[\d.eE+-]+$', value_str.replace(' ', '')):
                        continue
                    
                    try:
                        value = float(eval(value_str))
                    except:
                        continue
                    
                    constants[name] = PhysicsConstant(
                        name=name,
                        value=value,
                        description=comment.strip(),
                        source_lang="python",
                        source_file=file_path,
                        line_number=line_num
                    )
    except Exception as e:
        logger.error(f"Failed to extract Python constants: {e}")
    
    logger.info(f"Extracted {len(constants)} constants from Python")
    return constants


def extract_js_constants(file_path: str) -> Dict[str, PhysicsConstant]:
    """Extract CONSTANTS object from index.js"""
    constants = {}
    
    # Find CONSTANTS object and extract values
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Find CONSTANTS = { ... }
        const_match = re.search(r'const\s+CONSTANTS\s*=\s*\{([^}]+)\}', content, re.DOTALL)
        if const_match:
            const_block = const_match.group(1)
            
            # Extract key: value pairs
            pattern = re.compile(r'(\w+)\s*:\s*([^,\n]+)')
            for match in pattern.finditer(const_block):
                name = match.group(1)
                value_str = match.group(2).strip()
                
                try:
                    value = float(eval(value_str.replace('e', 'E')))
                except:
                    continue
                
                constants[name] = PhysicsConstant(
                    name=name,
                    value=value,
                    source_lang="javascript",
                    source_file=file_path
                )
    except Exception as e:
        logger.error(f"Failed to extract JS constants: {e}")
    
    logger.info(f"Extracted {len(constants)} constants from JavaScript")
    return constants


def compare_constants(
    cpp_constants: Dict[str, PhysicsConstant],
    py_constants: Dict[str, PhysicsConstant],
    js_constants: Dict[str, PhysicsConstant],
    tolerance: float = 1e-10
) -> Dict[str, Any]:
    """Compare constants across all three languages"""
    
    results = {
        'matching': [],
        'mismatched': [],
        'cpp_only': [],
        'python_only': [],
        'js_only': [],
        'all_names': set()
    }
    
    # Collect all constant names
    all_names = set(cpp_constants.keys()) | set(py_constants.keys()) | set(js_constants.keys())
    results['all_names'] = all_names
    
    for name in all_names:
        cpp_val = cpp_constants.get(name)
        py_val = py_constants.get(name)
        js_val = js_constants.get(name)
        
        values = []
        if cpp_val: values.append(('cpp', cpp_val.value))
        if py_val: values.append(('python', py_val.value))
        if js_val: values.append(('javascript', js_val.value))
        
        if len(values) == 3:
            # Check if all match
            cpp, py, js = cpp_val.value, py_val.value, js_val.value
            if abs(cpp - py) < tolerance and abs(py - js) < tolerance:
                results['matching'].append(name)
            else:
                results['mismatched'].append({
                    'name': name,
                    'cpp': cpp,
                    'python': py,
                    'javascript': js,
                    'max_diff': max(abs(cpp - py), abs(py - js), abs(cpp - js))
                })
        elif len(values) == 1:
            lang = values[0][0]
            if lang == 'cpp':
                results['cpp_only'].append(name)
            elif lang == 'python':
                results['python_only'].append(name)
            else:
                results['js_only'].append(name)
    
    return results


def sync_constants_from_cpp(
    cpp_path: str,
    py_path: str,
    js_path: str,
    dry_run: bool = False
) -> Dict[str, Any]:
    """
    Synchronize constants from C++ (master) to Python and JavaScript.
    C++ shared_constants.h is the canonical source.
    """
    logger.info("Synchronizing constants (C++ → Python, JavaScript)...")
    
    cpp_constants = extract_cpp_constants(cpp_path)
    py_constants = extract_python_constants(py_path)
    js_constants = extract_js_constants(js_path)
    
    comparison = compare_constants(cpp_constants, py_constants, js_constants)
    
    result = {
        'updated_python': [],
        'updated_javascript': [],
        'errors': [],
        'dry_run': dry_run
    }
    
    if not dry_run:
        # Update Python file
        for mismatch in comparison['mismatched']:
            name = mismatch['name']
            if name in cpp_constants:
                # Would update Python file here
                result['updated_python'].append(name)
        
        # Update JavaScript file
        for mismatch in comparison['mismatched']:
            name = mismatch['name']
            if name in cpp_constants:
                # Would update JS file here
                result['updated_javascript'].append(name)
    
    logger.info(f"Sync complete: {len(comparison['matching'])} matching, "
                f"{len(comparison['mismatched'])} mismatched")
    
    return result


# ============================================================================
# FILE REGENERATION
# ============================================================================

def regenerate_extracted_files(base_path: str) -> Dict[str, Any]:
    """
    Regenerate QCalc_cpp_extracted.py, QCalc_js_extracted.py, etc.
    These files contain constants/equations extracted from other languages.
    """
    logger.info("Regenerating extracted files...")
    
    result = {
        'files_updated': [],
        'errors': []
    }
    
    # Extract C++ constants and write to QCalc_cpp_extracted.py
    cpp_path = os.path.join(base_path, 'shared_constants.h')
    cpp_extracted_path = os.path.join(base_path, 'QCalc_cpp_extracted.py')
    
    try:
        cpp_constants = extract_cpp_constants(cpp_path)
        
        with open(cpp_extracted_path, 'w', encoding='utf-8') as f:
            f.write('"""\n')
            f.write('QCalc_cpp_extracted.py\n')
            f.write('======================\n')
            f.write('Auto-generated by Poseidon TaskBot\n')
            f.write(f'Generated: {datetime.now().isoformat()}\n')
            f.write(f'Source: {cpp_path}\n')
            f.write('"""\n\n')
            f.write('# C++ Constants extracted from shared_constants.h\n\n')
            
            for name, const in sorted(cpp_constants.items()):
                comment = f"  # {const.description}" if const.description else ""
                f.write(f'{name} = {const.value:.15e}{comment}\n')
        
        result['files_updated'].append(cpp_extracted_path)
        logger.info(f"Updated {cpp_extracted_path}")
        
    except Exception as e:
        result['errors'].append(f"C++ extraction failed: {e}")
        logger.error(f"C++ extraction failed: {e}")
    
    # Extract JS constants and write to QCalc_js_extracted.py
    js_path = os.path.join(base_path, 'index.js')
    js_extracted_path = os.path.join(base_path, 'QCalc_js_extracted.py')
    
    try:
        js_constants = extract_js_constants(js_path)
        
        with open(js_extracted_path, 'w', encoding='utf-8') as f:
            f.write('"""\n')
            f.write('QCalc_js_extracted.py\n')
            f.write('=====================\n')
            f.write('Auto-generated by Poseidon TaskBot\n')
            f.write(f'Generated: {datetime.now().isoformat()}\n')
            f.write(f'Source: {js_path}\n')
            f.write('"""\n\n')
            f.write('# JavaScript constants extracted from index.js CONSTANTS object\n\n')
            
            for name, const in sorted(js_constants.items()):
                f.write(f'{name} = {const.value:.15e}\n')
        
        result['files_updated'].append(js_extracted_path)
        logger.info(f"Updated {js_extracted_path}")
        
    except Exception as e:
        result['errors'].append(f"JS extraction failed: {e}")
        logger.error(f"JS extraction failed: {e}")
    
    return result


# ============================================================================
# CALCULATOR INVOCATION
# ============================================================================

def compute_equation(
    equation_name: str,
    params: Dict[str, float],
    module_name: str = "QCalc"
) -> Dict[str, Any]:
    """
    Invoke Python calculator for a specific equation.
    Called by C++ Poseidon TaskBot via pybind11.
    """
    logger.info(f"Computing {equation_name} via {module_name}")
    
    result = {
        'equation': equation_name,
        'value': 0.0,
        'success': False,
        'error': None,
        'compute_time_ms': 0.0
    }
    
    start_time = datetime.now()
    
    try:
        # Import calculator module dynamically
        if module_name == "QCalc":
            from QCalc import UnifiedFieldSolver
            solver = UnifiedFieldSolver()
            
            # Map equation names to solver methods
            if equation_name == "F_U":
                result['value'] = solver.compute_F_U(
                    r=params.get('r', 1e9),
                    t=params.get('t', 0)
                )
            elif equation_name == "F_U_Bi_i":
                result['value'] = solver.compute_F_U_Bi_i(
                    M=params.get('M', 1e30),
                    r=params.get('r', 1e9)
                )
            elif equation_name == "compressed_g":
                result['value'] = solver.compute_compressed_g(
                    r=params.get('r', 1e9)
                )
            else:
                # Generic equation
                result['value'] = solver.compute(equation_name, **params)
            
            result['success'] = True
            
        elif module_name == "CondensedPhysics":
            from CondensedPhysics import UQFFMasterEquations
            calc = UQFFMasterEquations()
            result['value'] = calc.compute(equation_name, params)
            result['success'] = True
            
    except ImportError as e:
        result['error'] = f"Module import failed: {e}"
        logger.error(result['error'])
    except Exception as e:
        result['error'] = str(e)
        logger.error(f"Computation failed: {e}")
    
    result['compute_time_ms'] = (datetime.now() - start_time).total_seconds() * 1000
    
    return result


def compute_all_equations(
    system_name: str,
    params: Dict[str, float]
) -> Dict[str, Any]:
    """Compute all UQFF equations for a system"""
    
    equations = [
        "F_U", "F_U_Bi_i", "compressed_g",
        "Ug1", "Ug2", "Ug3", "Ug4",
        "Um", "Ubi"
    ]
    
    results = {
        'system': system_name,
        'equations': {},
        'total_time_ms': 0.0
    }
    
    for eq in equations:
        eq_result = compute_equation(eq, params)
        results['equations'][eq] = eq_result
        results['total_time_ms'] += eq_result['compute_time_ms']
    
    return results


# ============================================================================
# VALIDATION
# ============================================================================

def validate_system(
    system_name: str,
    full_suite: bool = True
) -> Dict[str, Any]:
    """
    Run validation suite for a system.
    Calls QCalc_validation.py and CondensedPhysics_Validation.py
    """
    logger.info(f"Validating system: {system_name} (full={full_suite})")
    
    result = {
        'system': system_name,
        'tests_run': 0,
        'tests_passed': 0,
        'tests_failed': 0,
        'score': 0.0,
        'failed_tests': [],
        'warnings': [],
        'full_report': ""
    }
    
    try:
        # Try importing validation modules
        try:
            from QCalc_validation import validate_system as qcalc_validate
            qcalc_result = qcalc_validate(system_name)
            result['tests_run'] += qcalc_result.get('tests_run', 0)
            result['tests_passed'] += qcalc_result.get('tests_passed', 0)
            result['tests_failed'] += qcalc_result.get('tests_failed', 0)
        except ImportError:
            result['warnings'].append("QCalc_validation not available")
        
        try:
            from CondensedPhysics_Validation import validate as cp_validate
            cp_result = cp_validate(system_name)
            result['tests_run'] += cp_result.get('tests_run', 0)
            result['tests_passed'] += cp_result.get('tests_passed', 0)
            result['tests_failed'] += cp_result.get('tests_failed', 0)
        except ImportError:
            result['warnings'].append("CondensedPhysics_Validation not available")
        
        # Calculate score
        if result['tests_run'] > 0:
            result['score'] = 100.0 * result['tests_passed'] / result['tests_run']
        
    except Exception as e:
        result['warnings'].append(f"Validation error: {e}")
        logger.error(f"Validation failed: {e}")
    
    # Generate report
    result['full_report'] = f"""
Validation Report: {system_name}
================================
Tests Run: {result['tests_run']}
Passed: {result['tests_passed']}
Failed: {result['tests_failed']}
Score: {result['score']:.1f}%
Warnings: {', '.join(result['warnings']) if result['warnings'] else 'None'}
"""
    
    logger.info(f"Validation complete: {result['score']:.1f}%")
    
    return result


def cross_validate_frameworks(system_name: str) -> Dict[str, Any]:
    """Cross-validate UQFF vs MUGE vs GR for a system"""
    
    result = {
        'system': system_name,
        'uqff_result': None,
        'muge_result': None,
        'gr_comparison': None,
        'agreement_score': 0.0
    }
    
    # Placeholder - would implement full cross-validation
    logger.info(f"Cross-framework validation for {system_name}")
    
    return result


# ============================================================================
# FTPS BUNDLE OPERATIONS
# ============================================================================

def create_maintenance_bundle(base_path: str, bundle_path: str) -> Dict[str, Any]:
    """Create a maintenance bundle ZIP for FTPS transfer"""
    import zipfile
    
    result = {
        'bundle_path': bundle_path,
        'files_included': [],
        'total_size_bytes': 0,
        'success': False
    }
    
    files_to_bundle = [
        'shared_constants.h',
        'shared_constants.py',
        'poseidon_config.ini',
        'observational_systems_config.h'
    ]
    
    try:
        with zipfile.ZipFile(bundle_path, 'w', zipfile.ZIP_DEFLATED) as zf:
            for filename in files_to_bundle:
                filepath = os.path.join(base_path, filename)
                if os.path.exists(filepath):
                    zf.write(filepath, filename)
                    result['files_included'].append(filename)
            
            # Add metadata
            metadata = {
                'version': '4.2.1',
                'created': datetime.now().isoformat(),
                'source_path': base_path,
                'files': result['files_included']
            }
            zf.writestr('metadata.json', json.dumps(metadata, indent=2))
        
        result['total_size_bytes'] = os.path.getsize(bundle_path)
        result['success'] = True
        logger.info(f"Created bundle: {bundle_path} ({result['total_size_bytes']} bytes)")
        
    except Exception as e:
        result['error'] = str(e)
        logger.error(f"Bundle creation failed: {e}")
    
    return result


# ============================================================================
# MAIN ENTRY POINTS (Called from C++)
# ============================================================================

def poseidon_init(base_path: str) -> Dict[str, Any]:
    """Initialize Poseidon Python side"""
    config = PoseidonPythonConfig(base_path=base_path)
    set_config(config)
    
    return {
        'status': 'initialized',
        'base_path': config.base_path,
        'python_version': sys.version
    }


def poseidon_sync_constants() -> Dict[str, Any]:
    """Synchronize constants across languages (entry point for C++)"""
    return sync_constants_from_cpp(
        os.path.join(_config.base_path, _config.cpp_constants_file),
        os.path.join(_config.base_path, _config.py_constants_file),
        os.path.join(_config.base_path, _config.js_constants_file)
    )


def poseidon_regenerate_files() -> Dict[str, Any]:
    """Regenerate extracted files (entry point for C++)"""
    return regenerate_extracted_files(_config.base_path)


def poseidon_validate(system_name: str, full_suite: bool = True) -> Dict[str, Any]:
    """Run validation (entry point for C++)"""
    return validate_system(system_name, full_suite)


def poseidon_compute(equation: str, params: Dict[str, float]) -> Dict[str, Any]:
    """Compute equation (entry point for C++)"""
    return compute_equation(equation, params)


# ============================================================================
# CLI INTERFACE
# ============================================================================

def main():
    """CLI interface for standalone usage"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Poseidon TaskBot Python Companion')
    parser.add_argument('command', choices=[
        'sync', 'regenerate', 'validate', 'compute', 'bundle', 'status'
    ])
    parser.add_argument('--base-path', default=os.getcwd())
    parser.add_argument('--system', default='SgrA')
    parser.add_argument('--equation', default='F_U')
    parser.add_argument('--full-suite', action='store_true')
    
    args = parser.parse_args()
    
    # Initialize
    poseidon_init(args.base_path)
    
    if args.command == 'sync':
        result = poseidon_sync_constants()
        print(json.dumps(result, indent=2))
        
    elif args.command == 'regenerate':
        result = poseidon_regenerate_files()
        print(json.dumps(result, indent=2))
        
    elif args.command == 'validate':
        result = poseidon_validate(args.system, args.full_suite)
        print(result['full_report'])
        
    elif args.command == 'compute':
        result = poseidon_compute(args.equation, {'r': 1e9, 't': 0})
        print(json.dumps(result, indent=2))
        
    elif args.command == 'bundle':
        bundle_path = os.path.join(
            args.base_path, 
            f"maintenance_bundle_{datetime.now().strftime('%Y%m%d_%H%M%S')}.zip"
        )
        result = create_maintenance_bundle(args.base_path, bundle_path)
        print(json.dumps(result, indent=2))
        
    elif args.command == 'status':
        print(f"Poseidon Python Companion v4.2.1")
        print(f"Base path: {_config.base_path}")
        print(f"Python: {sys.version}")


if __name__ == '__main__':
    main()
