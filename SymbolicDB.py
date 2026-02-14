#!/usr/bin/env python3
"""
SymbolicDB.py - Symbolic Database + SymPy Solver Engine
========================================================

ARCHITECTURE: 4-Layer Stack (Option A - SQLite + SymPy, SIMPLEST)
─────────────────────────────────────────────────────────────────────────────
LAYER 1: EQUATION DATABASE (SQLite)
    - Store 5,000-8,000 equations as JSON metadata + SymPy expressions
    - Size: ~2 KB per equation = 10-15 MB total
    - Query: SQL-style searches by category, source, system type

LAYER 2: SYMPY SOLVER (Simple eval())
    - sympify(expr_str) → SymPy expression
    - expr.subs(params) → Substitute parameters
    - float(expr.evalf()) → Numerical solution
    - Performance: ~500 calc/sec (per user requirement: "just built and working")

LAYER 3: TEMPLATE FAMILIES
    - 1 template → N variants (system-specific instances)
    - Example: Mass evolution template × 100 systems = 100 equations

LAYER 4: EXECUTION ENGINE
    - Query interface (by category, source, system)
    - Result caching (optional, future optimization)
    - Self-expanding integration (port from C++ source7/source10)

DESIGN PHILOSOPHY:
─────────────────────────────────────────────────────────────────────────────
✓ COMPLETENESS OVER SPEED: "Doesn't have to be fast, just built and working"
✓ SCALABLE: 5,000-8,000 equations in 10-15 MB (vs 300,000 lines raw code)
✓ MAINTAINABLE: Equations as DATA, not code
✓ EXTENSIBLE: Easy to add new equations via JSON import
✓ QUERYABLE: SQL searches impossible with raw functions

DATA MODEL:
─────────────────────────────────────────────────────────────────────────────
equations:
    id              TEXT PRIMARY KEY    -- "magnetar_B_decay", "cluster_mass_evolution"
    sympy_expr      TEXT NOT NULL       -- "B0 * exp(-t / tau_B)"
    latex           TEXT                -- LaTeX rendering
    category        TEXT                -- "astrophysics.magnetar", "cosmology.expansion"
    subcategory     TEXT                -- "magnetic_field", "mass_evolution"
    parameters      TEXT                -- JSON: ["B0", "t", "tau_B"]
    constants       TEXT                -- JSON: {"c": 2.998e8, "G": 6.674e-11}
    units           TEXT                -- "T" (Tesla), "kg", "m/s²"
    source_file     TEXT                -- "source14.cpp", "source26.cpp"
    source_function TEXT                -- "calculate_magnetic_field_decay"
    description     TEXT                -- Human-readable physics explanation
    references      TEXT                -- Papers, observations, theoretical basis
    self_expand     INTEGER             -- 1 if from self-expanding framework, 0 otherwise
    created_date    TEXT                -- ISO timestamp
    version         TEXT                -- "1.0", "2.0-Enhanced"

templates:
    id              TEXT PRIMARY KEY    -- "mass_evolution_cluster"
    base_expr       TEXT NOT NULL       -- Template expression with placeholders
    parameters      TEXT                -- JSON: required parameters
    category        TEXT                -- Category for generated variants
    variant_count   INTEGER             -- Number of system-specific instances

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Created: February 13, 2026 (Parallel Track 1: Symbolic DB Prototype)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import sqlite3
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field
from datetime import datetime
import numpy as np

# Lazy import SymPy (only when needed for solve operations)
sympy = None  # Will be imported on first use

# ═══════════════════════════════════════════════════════════════════════════════
# SYMBOLIC DATABASE ENGINE
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class EquationMetadata:
    """Metadata for a symbolic equation in the database"""
    id: str
    sympy_expr: str
    latex: str = ""
    category: str = ""
    subcategory: str = ""
    parameters: List[str] = field(default_factory=list)
    constants: Dict[str, float] = field(default_factory=dict)
    units: str = ""
    source_file: str = ""
    source_function: str = ""
    description: str = ""
    refs: str = ""
    self_expand: bool = False
    created_date: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = "1.0"

@dataclass
class SymbolicResult:
    """Result from symbolic equation evaluation"""
    equation_id: str
    result: float
    long_form: str  # Equation with substituted values
    parameters_used: Dict[str, float]
    units: str
    category: str
    execution_time_ms: float

class SymbolicDatabase:
    """
    Symbolic Database Engine - Store and solve 5,000-8,000 equations
    
    Performance Target: 500 calculations/second (Option A - SQLite + SymPy)
    Storage Target: 10-15 MB for 5,000-8,000 equations
    """
    
    def __init__(self, db_path: str = 'equations.db'):
        """Initialize symbolic database engine"""
        self.db_path = db_path
        self.conn = sqlite3.connect(db_path)
        self.conn.row_factory = sqlite3.Row  # Access columns by name
        self._create_schema()
        
        # Self-expanding features (ported from C++ source7.cpp, source10.cpp)
        self.learning_rates: Dict[str, float] = {}  # Adaptive learning per equation
        self.dynamic_params: Dict[str, Dict[str, float]] = {}  # Runtime parameter injection
        self.metadata_cache: Dict[str, EquationMetadata] = {}  # Performance cache
        
    def _create_schema(self):
        """Create database schema (equations + templates tables)"""
        cursor = self.conn.cursor()
        
        # Main equations table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS equations (
                id              TEXT PRIMARY KEY,
                sympy_expr      TEXT NOT NULL,
                latex           TEXT,
                category        TEXT,
                subcategory     TEXT,
                parameters      TEXT,
                constants       TEXT,
                units           TEXT,
                source_file     TEXT,
                source_function TEXT,
                description     TEXT,
                refs            TEXT,
                self_expand     INTEGER DEFAULT 0,
                created_date    TEXT,
                version         TEXT DEFAULT '1.0'
            )
        ''')
        
        # Templates table (for generating equation families)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS templates (
                id              TEXT PRIMARY KEY,
                base_expr       TEXT NOT NULL,
                parameters      TEXT,
                category        TEXT,
                variant_count   INTEGER DEFAULT 0,
                description     TEXT
            )
        ''')
        
        # Create indexes for fast queries
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_category ON equations(category)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_source_file ON equations(source_file)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_self_expand ON equations(self_expand)')
        
        self.conn.commit()
    
    # ═══════════════════════════════════════════════════════════════════════════
    # EQUATION INSERTION AND MANAGEMENT
    # ═══════════════════════════════════════════════════════════════════════════
    
    def add_equation(self, eq: EquationMetadata) -> bool:
        """
        Add equation to database
        
        Args:
            eq: EquationMetadata with equation definition
            
        Returns:
            True if successful, False if equation already exists
        """
        cursor = self.conn.cursor()
        try:
            cursor.execute('''
                INSERT INTO equations VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (
                eq.id, eq.sympy_expr, eq.latex, eq.category, eq.subcategory,
                json.dumps(eq.parameters), json.dumps(eq.constants), eq.units,
                eq.source_file, eq.source_function, eq.description, eq.refs,
                1 if eq.self_expand else 0, eq.created_date, eq.version
            ))
            self.conn.commit()
            return True
        except sqlite3.IntegrityError:
            return False  # Equation ID already exists
    
    def batch_add_equations(self, equations: List[EquationMetadata]) -> Tuple[int, int]:
        """
        Add multiple equations in a single transaction
        
        Returns:
            (success_count, failure_count)
        """
        success = 0
        failure = 0
        for eq in equations:
            if self.add_equation(eq):
                success += 1
            else:
                failure += 1
        return success, failure
    
    def get_equation_metadata(self, eq_id: str) -> Optional[EquationMetadata]:
        """Retrieve equation metadata by ID"""
        # Check cache first
        if eq_id in self.metadata_cache:
            return self.metadata_cache[eq_id]
        
        cursor = self.conn.cursor()
        cursor.execute('SELECT * FROM equations WHERE id = ?', (eq_id,))
        row = cursor.fetchone()
        
        if not row:
            return None
        
        eq = EquationMetadata(
            id=row['id'],
            sympy_expr=row['sympy_expr'],
            latex=row['latex'] or "",
            category=row['category'] or "",
            subcategory=row['subcategory'] or "",
            parameters=json.loads(row['parameters']) if row['parameters'] else [],
            constants=json.loads(row['constants']) if row['constants'] else {},
            units=row['units'] or "",
            source_file=row['source_file'] or "",
            source_function=row['source_function'] or "",
            description=row['description'] or "",
            refs=row['refs'] or "",
            self_expand=bool(row['self_expand']),
            created_date=row['created_date'] or "",
            version=row['version'] or "1.0"
        )
        
        # Cache for future use
        self.metadata_cache[eq_id] = eq
        return eq
    
    # ═══════════════════════════════════════════════════════════════════════════
    # SYMPY SOLVER (LAYER 2)
    # ═══════════════════════════════════════════════════════════════════════════
    
    def solve(self, eq_id: str, param_values: Dict[str, float]) -> Optional[SymbolicResult]:
        """
        Solve equation with given parameters
        
        Args:
            eq_id: Equation ID
            param_values: Dict of parameter values
            
        Returns:
            SymbolicResult or None if equation not found
        """
        import time
        start_time = time.time()
        
        # Lazy import SymPy
        global sympy
        if sympy is None:
            import sympy as sp
            sympy = sp
        
        # Get equation metadata
        eq_meta = self.get_equation_metadata(eq_id)
        if not eq_meta:
            return None
        
        # Merge constants with provided parameters
        all_params = {**eq_meta.constants, **param_values}
        
        # Check for dynamic parameters (self-expanding feature)
        if eq_id in self.dynamic_params:
            all_params.update(self.dynamic_params[eq_id])
        
        # Parse SymPy expression
        try:
            expr = sympy.sympify(eq_meta.sympy_expr)
        except Exception as e:
            print(f"ERROR: Failed to parse expression for {eq_id}: {e}")
            return None
        
        # Substitute parameters
        try:
            result_expr = expr.subs(all_params)
            result_val = float(result_expr.evalf())
        except Exception as e:
            print(f"ERROR: Failed to evaluate {eq_id} with params {param_values}: {e}")
            return None
        
        # Build long-form equation string
        long_form = f"{eq_meta.sympy_expr}"
        for param, val in all_params.items():
            long_form += f"\n  {param} = {val}"
        long_form += f"\n  Result = {result_val} {eq_meta.units}"
        
        execution_time = (time.time() - start_time) * 1000  # Convert to ms
        
        return SymbolicResult(
            equation_id=eq_id,
            result=result_val,
            long_form=long_form,
            parameters_used=all_params,
            units=eq_meta.units,
            category=eq_meta.category,
            execution_time_ms=execution_time
        )
    
    # ═══════════════════════════════════════════════════════════════════════════
    # QUERY INTERFACE (LAYER 4)
    # ═══════════════════════════════════════════════════════════════════════════
    
    def query_by_category(self, category: str) -> List[str]:
        """Get all equation IDs in a category"""
        cursor = self.conn.cursor()
        cursor.execute('SELECT id FROM equations WHERE category LIKE ?', (f'%{category}%',))
        return [row['id'] for row in cursor.fetchall()]
    
    def query_by_source(self, source_file: str) -> List[str]:
        """Get all equation IDs from a source file"""
        cursor = self.conn.cursor()
        cursor.execute('SELECT id FROM equations WHERE source_file = ?', (source_file,))
        return [row['id'] for row in cursor.fetchall()]
    
    def query_self_expanding(self) -> List[str]:
        """Get all equations with self-expanding capabilities"""
        cursor = self.conn.cursor()
        cursor.execute('SELECT id FROM equations WHERE self_expand = 1')
        return [row['id'] for row in cursor.fetchall()]
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get database statistics"""
        cursor = self.conn.cursor()
        
        # Total equations
        cursor.execute('SELECT COUNT(*) FROM equations')
        total = cursor.fetchone()[0]
        
        # By category
        cursor.execute('SELECT category, COUNT(*) FROM equations GROUP BY category')
        by_category = {row[0]: row[1] for row in cursor.fetchall()}
        
        # By source file
        cursor.execute('SELECT source_file, COUNT(*) FROM equations GROUP BY source_file')
        by_source = {row[0]: row[1] for row in cursor.fetchall()}
        
        # Self-expanding count
        cursor.execute('SELECT COUNT(*) FROM equations WHERE self_expand = 1')
        self_expand_count = cursor.fetchone()[0]
        
        return {
            'total_equations': total,
            'by_category': by_category,
            'by_source': by_source,
            'self_expanding_count': self_expand_count
        }
    
    # ═══════════════════════════════════════════════════════════════════════════
    # SELF-EXPANDING FEATURES (Ported from C++ source7.cpp, source10.cpp)
    # ═══════════════════════════════════════════════════════════════════════════
    
    def set_dynamic_parameter(self, eq_id: str, param_name: str, value: float):
        """
        Set runtime parameter for equation (no recompilation needed)
        
        Ported from: source10.cpp line 401 (setDynamicParameter)
        """
        if eq_id not in self.dynamic_params:
            self.dynamic_params[eq_id] = {}
        self.dynamic_params[eq_id][param_name] = value
    
    def get_dynamic_parameter(self, eq_id: str, param_name: str, default: float = 0.0) -> float:
        """
        Get runtime parameter for equation
        
        Ported from: source10.cpp line 408 (getDynamicParameter)
        """
        if eq_id in self.dynamic_params and param_name in self.dynamic_params[eq_id]:
            return self.dynamic_params[eq_id][param_name]
        return default
    
    def set_learning_rate(self, eq_id: str, rate: float):
        """
        Set adaptive learning rate for equation optimization
        
        Ported from: source10.cpp line 118 (learning_rate = 0.01)
        """
        self.learning_rates[eq_id] = rate
    
    def export_state(self, filename: str):
        """
        Export dynamic state (learning rates, parameters) to JSON
        
        Ported from: source10.cpp line 350 (exportState)
        """
        state = {
            'version': '2.0-Enhanced',
            'timestamp': datetime.now().isoformat(),
            'learning_rates': self.learning_rates,
            'dynamic_params': self.dynamic_params
        }
        with open(filename, 'w') as f:
            json.dump(state, f, indent=2)
    
    def import_state(self, filename: str):
        """
        Import dynamic state from JSON
        
        Ported from: source10.cpp line 376 (importState)
        """
        with open(filename, 'r') as f:
            state = json.load(f)
        
        self.learning_rates = state.get('learning_rates', {})
        self.dynamic_params = state.get('dynamic_params', {})
    
    def close(self):
        """Close database connection"""
        self.conn.close()


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def create_symbolic_db(db_path: str = 'equations.db') -> SymbolicDatabase:
    """Create and initialize symbolic database"""
    return SymbolicDatabase(db_path)


if __name__ == '__main__':
    # Quick test
    db = create_symbolic_db('test_equations.db')
    
    # Add test equation
    eq = EquationMetadata(
        id='test_magnetic_decay',
        sympy_expr='B0 * exp(-t / tau_B)',
        latex='B(t) = B_0 e^{-t/\\tau_B}',
        category='astrophysics.magnetar',
        subcategory='magnetic_field',
        parameters=['B0', 't', 'tau_B'],
        units='T',
        source_file='source14.cpp',
        source_function='calculate_magnetic_field_decay',
        description='Exponential magnetic field decay for magnetars',
        refs='',
        self_expand=False
    )
    db.add_equation(eq)
    
    # Solve equation
    result = db.solve('test_magnetic_decay', {'B0': 1e10, 't': 1e11, 'tau_B': 1.26e14})
    if result:
        print(f"Result: {result.result:.3e} {result.units}")
        print(f"Execution time: {result.execution_time_ms:.2f} ms")
    
    # Statistics
    stats = db.get_statistics()
    print(f"\nDatabase statistics: {stats}")
    
    db.close()
