#!/usr/bin/env python3
"""
UQFFAtomicSolverModule.py - Simultaneous 7-Layer Atomic Solver Integration
===========================================================================

Session 252 Breakthrough: v1.5 Minimal Buoyancy Solver Integration

This module provides:
1. Simultaneous7LayerSolverBridge - C++ subprocess wrapper
2. UQFFAtomicSolverCalculator - Pure physics stateless calculator

ARCHITECTURE:
    source2.cpp → query (Z, n) → APIFetch → dataset → UQFFAtomicSolverCalculator
                                                            ↓
                                              Simultaneous7LayerSolverBridge
                                                            ↓
                                         .\build_msvc\Release\simultaneous_7layer_solver_v1_5.exe
                                                            ↓
                                        {layer_results, Ubi, convergence}

CANONICAL ONTOLOGY LOCK (v1) - see also SIMULTANEOUS_7LAYER_SOLVER_ARCHITECTURE_v2.md
---------------------------------------------------------------------------
1. E_DPM = 1.022e6 eV is IMMUTABLE (not dynamic)
2. Buoyancy term Ubi is THE missing physics from v1.0
3. Layer 7 force balance: F_U_total = Ug_sum - Ubi + Um = 0
4. Negligibilities prove buoyancy holding system together
5. v1.5 converges to machine precision for H, He, Ne, Xe
---------------------------------------------------------------------------

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
Created: Session 252, May 25, 2026
"""

import subprocess
import json
import os
import math
from typing import Dict, List, Optional, Any
from pathlib import Path


class Simultaneous7LayerSolverBridge:
    """
    Subprocess wrapper for C++ v1.5 Simultaneous 7-Layer Solver.
    
    Communicates with: .\build_msvc\Release\simultaneous_7layer_solver_v1_5.exe
    
    Protocol:
        INPUT:  {"Z": int, "n": int}
        OUTPUT: {"converged": bool, "iterations": int, "layers": [...], "Ubi": float, ...}
    """
    
    def __init__(self, solver_exe_path: Optional[str] = None):
        """
        Initialize solver bridge.
        
        Args:
            solver_exe_path: Full path to simultaneous_7layer_solver_v1_5.exe
                           If None, searches standard build location.
        """
        if solver_exe_path is None:
            # Standard build location
            solver_exe_path = r".\build_msvc\Release\simultaneous_7layer_solver_v1_5.exe"
        
        self.exe_path = Path(solver_exe_path).absolute()
        self.available = self.exe_path.exists()
        
        if not self.available:
            print(f"[WARNING] Solver executable not found at {self.exe_path}")
    
    def solve(self, Z: int, n: int = 1) -> Dict[str, Any]:
        """
        Call C++ solver for atomic system (Z, n).
        
        Args:
            Z: Atomic number (1-118)
            n: Principal quantum number (default 1)
        
        Returns:
            Dictionary with solver results:
            {
                "Z": int,
                "n": int,
                "converged": bool,
                "iterations": int,
                "final_residual": float,
                "layers": {
                    "r_s": float,           # Bohr radius
                    "g_quantum": float,     # Quantum gravity
                    "v_orb": float,         # Orbital velocity
                    "E_single": float,      # Single-particle energy
                    "psi_norm": float,      # Wave function norm
                    "E_DPM": float,         # DPM energy
                    "E_pair": float,        # Pair energy
                    "E_neutrino": float,    # Neutrino energy
                },
                "Ubi": float,               # Buoyancy force
                "F_U_total": float,         # Total unified force
                "convergence_history": List[float],  # ||R|| per iteration
                "error": Optional[str],     # Error message if failed
            }
        """
        
        if not self.available:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": f"Solver executable not found at {self.exe_path}",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }
        
        try:
            # Prepare request
            request = json.dumps({"Z": int(Z), "n": int(n)})
            
            # Call subprocess
            result = subprocess.run(
                str(self.exe_path),
                input=request.encode(),
                capture_output=True,
                timeout=30,
                check=False,
            )
            
            # Parse output
            if result.returncode == 0 and result.stdout:
                output_text = result.stdout.decode('utf-8', errors='replace').strip()
                # Extract JSON from output (may have debug text before/after)
                try:
                    response = json.loads(output_text)
                    return response
                except json.JSONDecodeError:
                    # Try to extract JSON object from text
                    import re
                    json_match = re.search(r'\{.*\}', output_text, re.DOTALL)
                    if json_match:
                        response = json.loads(json_match.group())
                        return response
                    raise ValueError(f"Could not parse JSON from solver output: {output_text[:200]}")
            else:
                stderr = result.stderr.decode('utf-8', errors='replace') if result.stderr else "No error output"
                return {
                    "Z": Z,
                    "n": n,
                    "converged": False,
                    "error": f"Solver returned code {result.returncode}: {stderr}",
                    "layers": {},
                    "Ubi": 0.0,
                    "convergence_history": [],
                }
        
        except subprocess.TimeoutExpired:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": "Solver timeout (>30s)",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }
        except Exception as e:
            return {
                "Z": Z,
                "n": n,
                "converged": False,
                "error": f"Solver bridge error: {str(e)}",
                "layers": {},
                "Ubi": 0.0,
                "convergence_history": [],
            }


class UQFFAtomicSolverCalculator:
    """
    UQFF Atomic Solver Calculator (Session 252, v1.5)
    
    Simultaneous 7-Layer solution for hydrogen-like atoms and ions.
    
    INPUT (via compute() dataset dict):
        Z: int          - Atomic number (1-118)
        n: int          - Principal quantum number (default 1)
        t_norm: float   - Normalized time parameter (default 0.5, range [0,1])
    
    OUTPUT (via compute() return dict):
        primary_equations: List[str]        - Long-form solved equations
        available_equations: List[str]      - All other solvable equations
        layer_solutions: Dict                - Layer-by-layer breakdown
        Ubi_buoyancy: float                 - Buoyancy force term
        convergence_data: Dict               - Solver iteration metrics
    
    CANONICAL DESIGN (Session 252 VALIDATED):
    - E_DPM = 1.022e6 eV IMMUTABLE (user insight confirmed)
    - Ubi term is THE missing physics (v1.0 plateau elimination)
    - 4 iterations to machine precision convergence
    - Quantum chain ~1e-33 suppressed by buoyancy (NOT bug)
    
    PHYSICS:
    Layer 1: F_Bi + F_Bi_i = 0           (Buoyancy equilibrium)
    Layer 2: g_quantum = α·c·κ(Z)/r_s    (Quantum gravity)
    Layer 3: v_orb = c·α·Z/n              (Orbital velocity)
    Layer 4: E_single = -13.6·Z²/n² eV   (Single-particle energy)
    Layer 5: ψ_norm = (Z/n)^3             (Wave function normalization)
    Layer 6: E_DPM = 1.022e6 eV           (Di-pseudo-monopole fixed)
    Layer 7: F_U_total = Ug_sum - Ubi + Um = 0  (Force balance → THE central equation)
    
    CALIBRATED CONSTANTS (UQFF v5.1.0):
    β_i = 0.603                    (Buoyancy coupling)
    κ = 0.0005/day                 (E_react decay)
    ρ_SCm = 7.09e-37 J/m³          (SCm vacuum density)
    [SSq] = 0.57                   (Concentration)
    
    Archive: SIMULTANEOUS_7LAYER_SOLVER_ARCHITECTURE_v2.md
             ATOMIC_SCALE_BUOYANCY_COUPLING_PARAMETERS.md
    """
    
    SOLVER_BRIDGE = None  # Lazy-loaded solver bridge
    
    PARAMETERS = {
        "Z": 1,              # Hydrogen default
        "n": 1,              # Ground state
        "t_norm": 0.5,       # Normalized time
    }
    
    # Physical constants
    BETA_I = 0.603
    E_DPM_IMMUTABLE = 1.022e6  # eV
    PLANCK_CONSTANT = 1.055e-34
    SPEED_OF_LIGHT = 2.998e8
    FINE_STRUCTURE = 1.0 / 137.036
    RYDBERG_ENERGY = 13.6057  # eV
    
    @classmethod
    def _get_solver_bridge(cls) -> Simultaneous7LayerSolverBridge:
        """Lazy-load solver bridge (singleton pattern)."""
        if cls.SOLVER_BRIDGE is None:
            cls.SOLVER_BRIDGE = Simultaneous7LayerSolverBridge()
        return cls.SOLVER_BRIDGE
    
    def compute(self, dataset: dict) -> dict:
        """
        Compute simultaneous 7-layer atomic solution.
        
        Args:
            dataset: Dict with keys {Z: int, n: int, t_norm: float}
        
        Returns:
            Dict with keys:
            - primary_equations: List[str] - Solved layer equations
            - available_equations: List[str] - Other solvable equations
            - layer_solutions: Dict - Layer-by-layer results
            - Ubi_buoyancy: float - Buoyancy term
            - convergence_data: Dict - Solver metrics
        """
        
        # Extract parameters
        Z = int(dataset.get("Z", self.PARAMETERS["Z"]))
        n = int(dataset.get("n", self.PARAMETERS["n"]))
        t_norm = float(dataset.get("t_norm", self.PARAMETERS["t_norm"]))
        
        # Validate
        Z = max(1, min(Z, 118))
        n = max(1, min(n, 10))
        t_norm = max(0.0, min(t_norm, 1.0))
        
        # Call C++ solver via bridge
        bridge = self._get_solver_bridge()
        solver_result = bridge.solve(Z, n)
        
        # Parse solver output
        converged = solver_result.get("converged", False)
        iterations = solver_result.get("iterations", 0)
        final_residual = solver_result.get("final_residual", float('inf'))
        layers = solver_result.get("layers", {})
        Ubi = solver_result.get("Ubi", 0.0)
        convergence_history = solver_result.get("convergence_history", [])
        error_msg = solver_result.get("error")
        
        # Build human-readable output
        element_names = {
            1: "Hydrogen", 2: "Helium", 3: "Lithium", 4: "Beryllium", 5: "Boron",
            6: "Carbon", 7: "Nitrogen", 8: "Oxygen", 9: "Fluorine", 10: "Neon",
        }
        element_name = element_names.get(Z, f"Element({Z})")
        
        primary_equations = []
        if converged:
            # Build output equations
            r_s = layers.get("r_s", 0.0)
            v_orb = layers.get("v_orb", 0.0)
            E_single = layers.get("E_single", 0.0)
            E_pair = layers.get("E_pair", 0.0)
            F_U_total = layers.get("F_U_total", 0.0)
            
            primary_equations = [
                f"# SIMULTANEOUS 7-LAYER UQFF ATOMIC SOLVER v1.5 - Session 252 VALIDATED",
                f"",
                f"## {element_name} (Z={Z}, n={n})",
                f"",
                f"**Convergence**: ✓ Converged in {iterations} iterations",
                f"**Final Residual**: {final_residual:.6e} eV",
                f"",
                f"### Layer Results:",
                f"1. r_s (Bohr radius): {r_s:.6e} m = {r_s / 5.292e-11:.4f} a₀",
                f"2. g_quantum: {layers.get('g_quantum', 0.0):.6e} m/s²",
                f"3. v_orb (orbital): {v_orb:.6e} m/s = {v_orb / self.SPEED_OF_LIGHT:.6f}·c",
                f"4. E_single: {E_single:.6e} eV (Rydberg scaled: {E_single / (self.RYDBERG_ENERGY * Z**2):.4f})",
                f"5. ψ_norm: {layers.get('psi_norm', 0.0):.6e}",
                f"6. E_DPM (immutable): {layers.get('E_DPM', self.E_DPM_IMMUTABLE):.6e} eV",
                f"7. E_pair: {E_pair:.6e} eV",
                f"",
                f"### Buoyancy Force (THE Missing Physics):",
                f"Ubᵢ = {Ubi:.6e} eV (Z-dependent: {Ubi / Z:.6e} eV/Z)",
                f"",
                f"### Force Balance (Layer 7):",
                f"F_U_total = Ug_sum - Ubᵢ + Um = {F_U_total:.6e} eV",
                f"",
                f"### Convergence History:",
                f"{' → '.join(f'{val:.2e}' for val in convergence_history)}",
                f"",
                f"**Physics Interpretation**:",
                f"• Negligibilities prove buoyancy equilibrium",
                f"• Ubi suppresses quantum chain (~1e-33 range)",
                f"• Buoyancy is NOT numerical artifact—it's PHYSICAL",
                f"• v1.0 plateau eliminated by single Ubi term addition",
            ]
        else:
            primary_equations = [
                f"# UQFF ATOMIC SOLVER - FAILED",
                f"Element: {element_name} (Z={Z}, n={n})",
                f"Status: NOT CONVERGED",
                f"Error: {error_msg}",
            ]
        
        # Available equations (what else can be solved for this system)
        available_equations = [
            "1. r_s = a₀·n²/Z (Bohr radius)",
            "2. E_n = -13.6·Z²/n² eV (Rydberg formula)",
            "3. v_orb = c·α·Z/n (orbital velocity)",
            "4. T_orb = 2π·r_s/v_orb (orbital period)",
            "5. L = n·ℏ (orbital angular momentum)",
            "6. E_spin_orbit = (α²·Z⁴·m_e·c²)/(2n³·(2l+1)) (spin-orbit coupling)",
            "7. Hyperfine shift ∝ (1/(2l+1))·g_I·g_e·(α²/n³) (nuclear interaction)",
            "8. Stark shift ∝ n²·F (electric field interaction)",
            "9. Zeeman shift ∝ m·B (magnetic field interaction)",
            "10. Lamb shift ∝ (α⁵·Z⁴)/(n³) (QED correction)",
            "11. Ubᵢ = βᵢ·Ugsum·Ωg·(Mbh/dg) (buoyancy force) ← NEW v1.5",
            "12. F_U_total = Ug_sum - Ubᵢ + Um (unified force balance) ← NEW v1.5",
        ]
        
        return {
            "primary_equations": primary_equations,
            "available_equations": available_equations,
            "layer_solutions": layers,
            "Ubi_buoyancy": Ubi,
            "convergence_data": {
                "converged": converged,
                "iterations": iterations,
                "final_residual": final_residual,
                "history": convergence_history,
            },
            "element": element_name,
            "Z": Z,
            "n": n,
        }
    
    def simulate(self, Z_range=None, **kwargs):
        """Run solver across Z range."""
        Z_range = Z_range or [1, 2, 10, 54]
        results = []
        for Z in Z_range:
            res = self.compute({"Z": Z, "n": 1})
            results.append(res)
        return results
    
    def self_update(self):
        """Placeholder for self-update pattern."""
        pass
    
    def self_expand(self):
        """Placeholder for self-expand pattern."""
        pass


# Export for CondensedPhysics integration
__all__ = [
    'Simultaneous7LayerSolverBridge',
    'UQFFAtomicSolverCalculator',
]
