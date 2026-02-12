#!/usr/bin/env python3
"""
OutputData.py - Output Storage for CondensedPhysics.py
=======================================================

Stores computed equation solutions for user recall.

DATA FLOW:
    InputData.py → CondensedPhysics.solve(params) → OutputData.py (this file)
                                                            ↓
                                            source2.cpp recall system

PURPOSE:
    - Store computed solutions from CondensedPhysics.py
    - Organize by query/timestamp for recall
    - Track equation sets used and computation metadata
    - Enable result comparisons across different methods

ARCHITECTURE:
    - JSON-based storage (human-readable, structured)
    - Organized by query name and timestamp
    - All results in SI units with metadata

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Created: February 11, 2026
"""

import json
import os
from datetime import datetime
from typing import Dict, List, Optional, Any
from dataclasses import dataclass, field, asdict


# ═══════════════════════════════════════════════════════════════════════════════
# COMPUTATION RESULT SCHEMA
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ComputationResult:
    """
    Stores a single computation result from CondensedPhysics.py.
    
    Attributes:
        query_name: Name of system queried (e.g., "Sagittarius A*")
        query_type: Type of calculation (e.g., "galaxy_rotation", "magnetar_field")
        timestamp: When computation was performed
        input_parameters: Dictionary of input parameters used
        long_form_equations: Step-by-step equation derivation
        solutions: Numerical results
        available_equations: List of other equations that could be solved
        simulation_set: Parameters for multi-equation simulation
        computation_time_ms: Time taken for computation
        method: Method used (e.g., "Newtonian", "NFW", "UQFF_Compressed")
    """
    query_name: str
    query_type: str
    timestamp: str
    input_parameters: Dict[str, float]
    long_form_equations: List[str]
    solutions: Dict[str, float]
    available_equations: List[str]
    simulation_set: Dict[str, Any]
    computation_time_ms: float
    method: str
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for JSON storage."""
        return asdict(self)


# ═══════════════════════════════════════════════════════════════════════════════
# OUTPUT DATA MANAGER
# ═══════════════════════════════════════════════════════════════════════════════

class OutputDataManager:
    """
    Manages output storage for CondensedPhysics.py computations.
    
    Responsibilities:
        1. Store computed solutions with metadata
        2. Organize by query name and timestamp
        3. Enable result recall for users
        4. Track computation history
    """
    
    def __init__(self, output_dir: str = "output_data"):
        """
        Initialize output data manager.
        
        Args:
            output_dir: Directory for storing output JSON files
        """
        self.output_dir = output_dir
        os.makedirs(output_dir, exist_ok=True)
        self.index_file = os.path.join(output_dir, "index.json")
        self._load_index()
    
    def _load_index(self):
        """Load index of all stored computations."""
        if os.path.exists(self.index_file):
            with open(self.index_file, 'r') as f:
                self.index = json.load(f)
        else:
            self.index = {
                'queries': {},  # {query_name: [list of computation IDs]}
                'computations': {}  # {computation_id: metadata}
            }
    
    def _save_index(self):
        """Save index to disk."""
        with open(self.index_file, 'w') as f:
            json.dump(self.index, f, indent=2)
    
    def store_result(self, result: ComputationResult) -> str:
        """
        Store computation result.
        
        Args:
            result: ComputationResult object with solutions
        
        Returns:
            Computation ID (unique identifier)
        """
        # Generate unique computation ID
        comp_id = f"{result.query_name.replace(' ', '_')}_{result.timestamp.replace(':', '').replace('-', '')}"
        
        # Store full result as JSON file
        result_file = os.path.join(self.output_dir, f"{comp_id}.json")
        with open(result_file, 'w') as f:
            json.dump(result.to_dict(), f, indent=2)
        
        # Update index
        if result.query_name not in self.index['queries']:
            self.index['queries'][result.query_name] = []
        
        self.index['queries'][result.query_name].append(comp_id)
        
        self.index['computations'][comp_id] = {
            'query_name': result.query_name,
            'query_type': result.query_type,
            'timestamp': result.timestamp,
            'method': result.method,
            'file': result_file
        }
        
        self._save_index()
        
        return comp_id
    
    def recall_by_query_name(self, query_name: str) -> List[ComputationResult]:
        """
        Recall all computations for a specific query name.
        
        Args:
            query_name: Name of system (e.g., "Sagittarius A*")
        
        Returns:
            List of ComputationResult objects (newest first)
        """
        if query_name not in self.index['queries']:
            return []
        
        results = []
        comp_ids = self.index['queries'][query_name]
        
        for comp_id in reversed(comp_ids):  # Newest first
            result_file = self.index['computations'][comp_id]['file']
            with open(result_file, 'r') as f:
                data = json.load(f)
                results.append(ComputationResult(**data))
        
        return results
    
    def recall_by_id(self, comp_id: str) -> Optional[ComputationResult]:
        """
        Recall specific computation by ID.
        
        Args:
            comp_id: Computation ID
        
        Returns:
            ComputationResult object or None
        """
        if comp_id not in self.index['computations']:
            return None
        
        result_file = self.index['computations'][comp_id]['file']
        with open(result_file, 'r') as f:
            data = json.load(f)
            return ComputationResult(**data)
    
    def list_all_queries(self) -> List[str]:
        """
        List all query names in storage.
        
        Returns:
            List of query names
        """
        return list(self.index['queries'].keys())
    
    def get_query_history(self, query_name: str) -> List[Dict[str, Any]]:
        """
        Get computation history for a query (metadata only).
        
        Args:
            query_name: Name of system
        
        Returns:
            List of metadata dictionaries (timestamp, method, type)
        """
        if query_name not in self.index['queries']:
            return []
        
        comp_ids = self.index['queries'][query_name]
        return [self.index['computations'][cid] for cid in reversed(comp_ids)]
    
    def compare_methods(self, query_name: str) -> Dict[str, Any]:
        """
        Compare different methods for the same query.
        
        Args:
            query_name: Name of system
        
        Returns:
            Dictionary with method comparisons
        """
        results = self.recall_by_query_name(query_name)
        
        if not results:
            return {'error': f'No computations found for {query_name}'}
        
        comparison = {
            'query_name': query_name,
            'num_computations': len(results),
            'methods_used': {},
            'solution_comparison': {}
        }
        
        for result in results:
            method = result.method
            
            if method not in comparison['methods_used']:
                comparison['methods_used'][method] = {
                    'count': 0,
                    'avg_computation_time_ms': 0,
                    'solutions': []
                }
            
            comparison['methods_used'][method]['count'] += 1
            comparison['methods_used'][method]['solutions'].append(result.solutions)
        
        # Calculate averages
        for method_data in comparison['methods_used'].values():
            if method_data['solutions']:
                method_data['avg_computation_time_ms'] = sum(
                    r.computation_time_ms for r in results if r.method == method
                ) / method_data['count']
        
        return comparison
    
    def export_query_results(self, query_name: str, output_file: str):
        """
        Export all results for a query to single JSON file.
        
        Args:
            query_name: Name of system
            output_file: Path to output JSON file
        """
        results = self.recall_by_query_name(query_name)
        
        export_data = {
            'query_name': query_name,
            'num_results': len(results),
            'export_timestamp': datetime.now().isoformat(),
            'results': [r.to_dict() for r in results]
        }
        
        with open(output_file, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        print(f"Exported {len(results)} results to {output_file}")


# ═══════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("OutputData.py - Output Storage Manager")
    print("=" * 70)
    
    # Create manager
    manager = OutputDataManager()
    
    # Example: Store a computation result
    result = ComputationResult(
        query_name="Sagittarius A*",
        query_type="galaxy_rotation",
        timestamp=datetime.now().isoformat(),
        input_parameters={
            'M': 4.1e6 * 1.989e30,
            'r': 26000 * 3.086e16
        },
        long_form_equations=[
            "v(r) = sqrt(G * M / r)",
            "v(r) = sqrt(6.674e-11 * 8.156e36 / 8.02e20)",
            "v(r) = 2.31e5 m/s"
        ],
        solutions={
            'v_rotation': 2.31e5,
            'v_rotation_km_s': 231.0
        },
        available_equations=[
            "NFW dark matter profile",
            "Tully-Fisher relation",
            "MOND dynamics"
        ],
        simulation_set={
            'equations': ['v(r) = sqrt(GM/r)'],
            'parameters': {'M': 8.156e36, 'r_range': (1e19, 1e21)}
        },
        computation_time_ms=12.5,
        method="Newtonian"
    )
    
    comp_id = manager.store_result(result)
    print(f"Stored computation: {comp_id}")
    
    # Recall by query name
    recalled = manager.recall_by_query_name("Sagittarius A*")
    print(f"\nRecalled {len(recalled)} result(s) for Sagittarius A*")
    
    if recalled:
        print(f"Latest computation:")
        print(f"  Method: {recalled[0].method}")
        print(f"  Solutions: {recalled[0].solutions}")
    
    # List all queries
    queries = manager.list_all_queries()
    print(f"\nAll queries in storage: {queries}")
    
    # Get history
    history = manager.get_query_history("Sagittarius A*")
    print(f"\nComputation history:")
    for h in history:
        print(f"  {h['timestamp']} - {h['method']} ({h['query_type']})")
