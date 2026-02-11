#!/usr/bin/env python3
"""
OPData.py - UQFF Output Data Storage
====================================

Stores computed solutions from QCalc.py for recall.

ARCHITECTURE:
    APIFetch.py → QCalc.py → OPData.py (this file)
                              ↓
                         User recall via query

STORAGE STRUCTURE:
    QUERY_RESULTS = {
        'query_id': {
            'timestamp': str,
            'input_params': dict,
            'long_form_equations': list,
            'solutions': dict,
            'available_equations': list
        }
    }

Author: Daniel T. Murphy (daniel.murphy00@gmail.com)
Framework: UQFF 99.9% Solvability (Star-Magic)
Copyright: © 2025-2026 Daniel T. Murphy - All Rights Reserved
"""

import json
import os
from datetime import datetime
from typing import Dict, List, Any, Optional

# ═══════════════════════════════════════════════════════════════════════════════
# IN-MEMORY QUERY RESULTS STORAGE
# ═══════════════════════════════════════════════════════════════════════════════

QUERY_RESULTS: Dict[str, Dict[str, Any]] = {}


# ═══════════════════════════════════════════════════════════════════════════════
# STORAGE CLASS
# ═══════════════════════════════════════════════════════════════════════════════

class OutputDataStore:
    """
    Manages storage and retrieval of computed UQFF solutions.
    
    Features:
    - In-memory storage for current session
    - JSON file persistence for long-term storage
    - Query-based recall by ID or search
    """
    
    def __init__(self, storage_file: str = "uqff_results.json"):
        """
        Initialize output data store.
        
        Args:
            storage_file: Path to JSON file for persistent storage
        """
        self.storage_file = storage_file
        self.results: Dict[str, Dict[str, Any]] = {}
        self._load_from_file()
    
    def store(self, result: Dict[str, Any]) -> str:
        """
        Store a computation result.
        
        Args:
            result: Complete result dict from QCalc.solve()
            
        Returns:
            query_id for later retrieval
        """
        query_id = result.get('query_id', f"query_{datetime.now().timestamp()}")
        
        self.results[query_id] = {
            'timestamp': result.get('timestamp', datetime.now().isoformat()),
            'input_params': result.get('input_params', {}),
            'long_form_equations': result.get('long_form_equations', []),
            'solutions': result.get('solutions', {}),
            'available_equations': result.get('available_equations', []),
            'simulation_set': result.get('simulation_set', {})
        }
        
        # Update global storage
        QUERY_RESULTS[query_id] = self.results[query_id]
        
        # Persist to file
        self._save_to_file()
        
        return query_id
    
    def recall(self, query_id: str) -> Optional[Dict[str, Any]]:
        """
        Recall a stored result by query_id.
        
        Args:
            query_id: The ID returned when result was stored
            
        Returns:
            Complete result dict or None if not found
        """
        return self.results.get(query_id)
    
    def search(self, search_term: str) -> List[Dict[str, Any]]:
        """
        Search stored results by name or parameter.
        
        Args:
            search_term: String to search for in query names/params
            
        Returns:
            List of matching results
        """
        matches = []
        search_lower = search_term.lower()
        
        for query_id, result in self.results.items():
            # Search in query_id
            if search_lower in query_id.lower():
                matches.append({'query_id': query_id, **result})
                continue
            
            # Search in input_params
            params = result.get('input_params', {})
            if any(search_lower in str(v).lower() for v in params.values()):
                matches.append({'query_id': query_id, **result})
        
        return matches
    
    def list_all(self) -> List[str]:
        """List all stored query IDs."""
        return list(self.results.keys())
    
    def get_latest(self, n: int = 10) -> List[Dict[str, Any]]:
        """
        Get the n most recent results.
        
        Args:
            n: Number of results to return
            
        Returns:
            List of recent results sorted by timestamp
        """
        sorted_results = sorted(
            [{'query_id': k, **v} for k, v in self.results.items()],
            key=lambda x: x.get('timestamp', ''),
            reverse=True
        )
        return sorted_results[:n]
    
    def delete(self, query_id: str) -> bool:
        """
        Delete a stored result.
        
        Args:
            query_id: The ID of the result to delete
            
        Returns:
            True if deleted, False if not found
        """
        if query_id in self.results:
            del self.results[query_id]
            if query_id in QUERY_RESULTS:
                del QUERY_RESULTS[query_id]
            self._save_to_file()
            return True
        return False
    
    def clear_all(self) -> int:
        """
        Clear all stored results.
        
        Returns:
            Number of results cleared
        """
        count = len(self.results)
        self.results.clear()
        QUERY_RESULTS.clear()
        self._save_to_file()
        return count
    
    def _save_to_file(self):
        """Persist results to JSON file."""
        try:
            with open(self.storage_file, 'w') as f:
                json.dump(self.results, f, indent=2, default=str)
        except Exception as e:
            print(f"Warning: Could not save to {self.storage_file}: {e}")
    
    def _load_from_file(self):
        """Load results from JSON file if exists."""
        if os.path.exists(self.storage_file):
            try:
                with open(self.storage_file, 'r') as f:
                    self.results = json.load(f)
                    # Update global storage
                    QUERY_RESULTS.update(self.results)
            except Exception as e:
                print(f"Warning: Could not load from {self.storage_file}: {e}")
    
    def export_to_csv(self, output_file: str = "uqff_solutions.csv"):
        """
        Export solutions to CSV format.
        
        Args:
            output_file: Path to output CSV file
        """
        import csv
        
        rows = []
        for query_id, result in self.results.items():
            base_row = {
                'query_id': query_id,
                'timestamp': result.get('timestamp', ''),
            }
            # Add input params
            for k, v in result.get('input_params', {}).items():
                base_row[f'input_{k}'] = v
            # Add solutions
            for k, v in result.get('solutions', {}).items():
                base_row[f'solution_{k}'] = v
            rows.append(base_row)
        
        if rows:
            fieldnames = set()
            for row in rows:
                fieldnames.update(row.keys())
            
            with open(output_file, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=sorted(fieldnames))
                writer.writeheader()
                writer.writerows(rows)
            
            print(f"Exported {len(rows)} results to {output_file}")
    
    def get_statistics(self) -> Dict[str, Any]:
        """Get statistics about stored results."""
        if not self.results:
            return {'count': 0, 'message': 'No results stored'}
        
        timestamps = [r.get('timestamp', '') for r in self.results.values()]
        timestamps = [t for t in timestamps if t]
        
        return {
            'count': len(self.results),
            'earliest': min(timestamps) if timestamps else None,
            'latest': max(timestamps) if timestamps else None,
            'query_ids': list(self.results.keys())[:10],  # First 10
        }


# ═══════════════════════════════════════════════════════════════════════════════
# GLOBAL STORE INSTANCE
# ═══════════════════════════════════════════════════════════════════════════════

STORE = OutputDataStore()


# ═══════════════════════════════════════════════════════════════════════════════
# CONVENIENCE FUNCTIONS
# ═══════════════════════════════════════════════════════════════════════════════

def store(result: Dict[str, Any]) -> str:
    """Store a result. Returns query_id."""
    return STORE.store(result)


def recall(query_id: str) -> Optional[Dict[str, Any]]:
    """Recall a result by query_id."""
    return STORE.recall(query_id)


def search(term: str) -> List[Dict[str, Any]]:
    """Search results by term."""
    return STORE.search(term)


def list_queries() -> List[str]:
    """List all query IDs."""
    return STORE.list_all()


def get_latest(n: int = 10) -> List[Dict[str, Any]]:
    """Get n most recent results."""
    return STORE.get_latest(n)


# ═══════════════════════════════════════════════════════════════════════════════
# MODULE TEST
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("=" * 80)
    print("OPData.py - Output Data Storage Test")
    print("=" * 80)
    
    # Test storage
    test_result = {
        'query_id': 'test_query_001',
        'timestamp': datetime.now().isoformat(),
        'input_params': {'M': 1.989e30, 'r': 1.496e11, 'T': 5778},
        'long_form_equations': [
            {'name': 'Ug1', 'result': 9.8, 'unit': 'm/s²'}
        ],
        'solutions': {'Ug1': 9.8, 'Ub': -0.5, 'F_U': 9.3},
        'available_equations': ['compute_escape_velocity', 'compute_luminosity']
    }
    
    # Store
    query_id = store(test_result)
    print(f"Stored result with ID: {query_id}")
    
    # Recall
    recalled = recall(query_id)
    print(f"Recalled: {recalled is not None}")
    
    # List
    all_ids = list_queries()
    print(f"All query IDs: {all_ids}")
    
    # Statistics
    stats = STORE.get_statistics()
    print(f"Statistics: {stats}")
