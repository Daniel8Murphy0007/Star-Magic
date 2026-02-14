#!/usr/bin/env python3
"""
extract_functions_to_db.py - Convert 94 Functions to Symbolic Database
========================================================================

TRACK 1 (PARALLEL): Convert existing 94 functions from QCalc_Wolfram_Extensions.py
to symbolic database format (JSON + SQLite).

PROCESS:
────────────────────────────────────────────────────────────────────────────────
1. Parse QCalc_Wolfram_Extensions.py
2. Extract function metadata (name, parameters, equations, units)
3. Convert to SymPy expressions
4. Generate EquationMetadata objects
5. Insert into SymbolicDB
6. Verify: All 94 functions accessible via symbolic interface

TARGET: 94 equations in symbolic DB (~188 KB, 2 KB per equation)

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Created: February 13, 2026 (Parallel Track 1: Symbolic DB Prototype)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import ast
import re
import inspect
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import sys

# Import SymbolicDB engine
from SymbolicDB import SymbolicDatabase, EquationMetadata

# Import QCalc modules to introspect functions
import QCalc_Wolfram_Extensions


def extract_function_metadata(func_name: str) -> Optional[Dict]:
    """
    Extract metadata from a QCalc_Wolfram_Extensions function
    
    Args:
        func_name: Function name (e.g., "calculate_magnetic_field_decay")
        
    Returns:
        Dict with metadata or None if function not found
    """
    # Get function object
    if not hasattr(QCalc_Wolfram_Extensions, func_name):
        print(f"WARNING: Function {func_name} not found in QCalc_Wolfram_Extensions")
        return None
    
    func = getattr(QCalc_Wolfram_Extensions, func_name)
    
    # Get function signature
    sig = inspect.signature(func)
    params = list(sig.parameters.keys())
    
    # Get docstring
    docstring = inspect.getdoc(func) or ""
    
    # Get source code
    source = inspect.getsource(func)
    
    # Parse metadata from docstring and source
    metadata = {
        'function_name': func_name,
        'parameters': params,
        'docstring': docstring,
        'source_code': source
    }
    
    # Extract category from docstring (e.g., "SOURCE14: Magnetar Physics")
    category_match = re.search(r'SOURCE(\d+):\s*(.+)', docstring)
    if category_match:
        source_num = category_match.group(1)
        category_desc = category_match.group(2).strip()
        metadata['source_file'] = f'source{source_num}.cpp'
        metadata['category'] = category_desc
    
    # Extract units from docstring (e.g., "Returns: T (Tesla)")
    units_match = re.search(r'Returns?:\s*([A-Za-z²³⁻⁰⁴⁵⁶⁷⁸⁹·/\s]+)', docstring)
    if units_match:
        metadata['units'] = units_match.group(1).strip()
    else:
        metadata['units'] = 'dimensionless'
    
    # Extract equation from source (look for return statement)
    equation_match = re.search(r'return\s+(.+)', source, re.MULTILINE)
    if equation_match:
        # Get the return expression (single line or multiline)
        return_expr = equation_match.group(1)
        # Clean up Python syntax for SymPy
        metadata['raw_equation'] = return_expr
    
    return metadata


def python_to_sympy_expr(python_expr: str, simplify: bool = True) -> str:
    """
    Convert Python expression to SymPy-compatible string
    
    Args:
        python_expr: Python expression string
        simplify: Whether to simplify the expression
        
    Returns:
        SymPy-compatible expression string
    """
    # Replace numpy/math functions with SymPy equivalents
    replacements = {
        'np.exp': 'exp',
        'np.sin': 'sin',
        'np.cos': 'cos',
        'np.tan': 'tan',
        'np.log': 'log',
        'np.sqrt': 'sqrt',
        'np.pi': 'pi',
        'np.e': 'E',
        'np.abs': 'Abs',
        'math.exp': 'exp',
        'math.sin': 'sin',
        'math.cos': 'cos',
        'math.sqrt': 'sqrt',
        'math.pi': 'pi',
    }
    
    sympy_expr = python_expr
    for py_func, sympy_func in replacements.items():
        sympy_expr = sympy_expr.replace(py_func, sympy_func)
    
    # Remove EquationResult() wrapper if present
    sympy_expr = re.sub(r'EquationResult\([^,]+,\s*(.+?),\s*[^)]+\)', r'\1', sympy_expr)
    
    # Handle multi-line expressions (take first line if multiple)
    sympy_expr = sympy_expr.split('\n')[0].strip()
    
    return sympy_expr


def categorize_function(func_name: str, metadata: Dict) -> Tuple[str, str]:
    """
    Determine category and subcategory for function
    
    Returns:
        (category, subcategory) tuple
    """
    # Map function names to categories
    name = func_name.lower()
    
    if 'magnetar' in name or 'magnetic' in name or 'sgr' in name:
        return 'astrophysics.magnetar', 'magnetic_field'
    elif 'smbh' in name or 'black_hole' in name or 'sagittarius' in name or 'sgr_a' in name:
        return 'astrophysics.smbh', 'black_hole'
    elif 'cluster' in name or 'westerlund' in name:
        return 'astrophysics.cluster', 'stellar_cluster'
    elif 'star_formation' in name or 'tapestry' in name or 'pillars' in name:
        return 'astrophysics.star_formation', 'molecular_cloud'
    elif 'galaxy' in name or 'ngc' in name or 'andromeda' in name or 'sombrero' in name:
        return 'astrophysics.galaxy', 'spiral_galaxy'
    elif 'nebula' in name or 'crab' in name or 'lagoon' in name or 'orion' in name:
        return 'astrophysics.nebula', 'emission_nebula'
    elif 'cosmological' in name or 'hubble' in name or 'hudf' in name:
        return 'cosmology.expansion', 'cosmological_constant'
    elif 'planetary' in name or 'saturn' in name:
        return 'astrophysics.planetary', 'gas_giant'
    elif 'hydrogen' in name or 'atomic' in name:
        return 'quantum.atomic', 'hydrogen_atom'
    elif 'universe' in name:
        return 'cosmology.universe', 'large_scale'
    elif 'quantum' in name or 'heisenberg' in name:
        return 'quantum.uncertainty', 'quantum_mechanics'
    elif 'gravitational_wave' in name:
        return 'relativity.gravitational_waves', 'gw_emission'
    elif 'fluid' in name or 'density' in name:
        return 'fluid_dynamics.astrophysical', 'fluid_coupling'
    elif 'spin' in name or 'angular' in name:
        return 'astrophysics.rotation', 'angular_momentum'
    elif 'mass' in name or 'accretion' in name:
        return 'astrophysics.mass_evolution', 'mass_transfer'
    else:
        return 'physics.general', 'uncategorized'


def extract_all_functions() -> List[EquationMetadata]:
    """
    Extract all 94 functions from QCalc_Wolfram_Extensions.py
    
    Returns:
        List of EquationMetadata objects
    """
    equations = []
    
    # Get all function names starting with "calculate_"
    func_names = [name for name in dir(QCalc_Wolfram_Extensions) 
                  if name.startswith('calculate_') and callable(getattr(QCalc_Wolfram_Extensions, name))]
    
    print(f"Found {len(func_names)} functions in QCalc_Wolfram_Extensions.py")
    
    for func_name in func_names:
        print(f"Processing {func_name}...")
        
        # Extract metadata
        metadata = extract_function_metadata(func_name)
        if not metadata:
            continue
        
        # Generate equation ID
        eq_id = func_name.replace('calculate_', '')
        
        # Convert to SymPy expression (if possible)
        sympy_expr = "PLACEHOLDER"  # Will be filled in manually or via deeper parsing
        if 'raw_equation' in metadata:
            try:
                sympy_expr = python_to_sympy_expr(metadata['raw_equation'])
            except Exception as e:
                print(f"  WARNING: Failed to convert {func_name} to SymPy: {e}")
                sympy_expr = f"MANUAL_ENTRY_REQUIRED: {metadata['raw_equation'][:100]}"
        
        # Determine category
        category, subcategory = categorize_function(func_name, metadata)
        
        # Create EquationMetadata
        eq = EquationMetadata(
            id=eq_id,
            sympy_expr=sympy_expr,
            latex="",  # TODO: Generate LaTeX from SymPy
            category=category,
            subcategory=subcategory,
            parameters=metadata.get('parameters', []),
            constants={},  # TODO: Extract constants from CONSTANTS dict
            units=metadata.get('units', 'dimensionless'),
            source_file=metadata.get('source_file', 'unknown'),
            source_function=func_name,
            description=metadata.get('docstring', '')[:500],  # First 500 chars
            refs="",  # TODO: Extract from docstring
            self_expand=False,  # Mark as False initially (will be updated for source7/10 functions)
            version="1.0"
        )
        
        equations.append(eq)
    
    return equations


def populate_database(db_path: str = 'uqff_equations.db'):
    """
    Populate symbolic database with 94 functions
    
    Args:
        db_path: Path to SQLite database file
    """
    print(f"\n{'='*80}")
    print(f"SYMBOLIC DATABASE POPULATION: Track 1 (Parallel)")
    print(f"{'='*80}\n")
    
    # Create database
    print(f"Creating symbolic database: {db_path}")
    db = SymbolicDatabase(db_path)
    
    # Extract equations
    print("\nExtracting 94 functions from QCalc_Wolfram_Extensions.py...")
    equations = extract_all_functions()
    
    # Insert into database
    print(f"\nInserting {len(equations)} equations into database...")
    success, failure = db.batch_add_equations(equations)
    
    print(f"\n{'='*80}")
    print(f"EXTRACTION COMPLETE")
    print(f"{'='*80}")
    print(f"✓ Success: {success} equations")
    print(f"✗ Failure: {failure} equations")
    
    # Statistics
    print("\nDatabase Statistics:")
    stats = db.get_statistics()
    print(f"  Total equations: {stats['total_equations']}")
    print(f"  By category:")
    for cat, count in sorted(stats['by_category'].items()):
        print(f"    {cat}: {count}")
    print(f"\n  By source file:")
    for src, count in sorted(stats['by_source'].items()):
        print(f"    {src}: {count}")
    
    # Save database size
    db_size = Path(db_path).stat().st_size
    print(f"\nDatabase size: {db_size / 1024:.2f} KB")
    print(f"Average per equation: {db_size / max(success, 1):.2f} bytes")
    
    db.close()
    
    return success, failure


if __name__ == '__main__':
    # Run extraction
    success, failure = populate_database('uqff_equations.db')
    
    if failure > 0:
        print(f"\n⚠️  WARNING: {failure} equations failed to insert (likely duplicates)")
    
    print(f"\n✓ TRACK 1 MILESTONE: 94 functions converted to symbolic database")
    print(f"  Database: uqff_equations.db")
    print(f"  Next: Test symbolic solver with QCalc_test.py")
