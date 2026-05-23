#!/usr/bin/env python3
"""
_merge_trinity_derivations.py — Merge extracted derivations into Trinity

Reads extracted derivations from master_closures.csv and intelligently merges them
into _uqff_derivation_equations.py registry.

Strategy:
  1. Extract all 41 derivations from CSV
  2. Create DerivationEquation objects
  3. Group by domain
  4. Generate Python code for registry insertion
  5. Show merge preview
"""

import csv
import re
from typing import Dict, List, Tuple
from pathlib import Path


def load_extracted_derivations() -> List[Dict]:
    """Load derivations already extracted by _populate_trinity_bulk.py"""
    derivations = []
    csv_path = "master_closures.csv"
    
    print("Loading derivations from master_closures.csv...")
    
    with open(csv_path, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            constant_name = row.get('closure', '').strip()
            equation_latex = row.get('raw_output', '').strip()
            status = row.get('status', 'UNKNOWN').strip()
            error_pct = row.get('error_pct', '0')
            category = row.get('category', 'UNKNOWN').strip()
            label = row.get('label', '').strip()
            predicted = row.get('predicted', '0')
            
            # Skip empty rows or paper indices
            if not constant_name or constant_name.startswith('PAPER_'):
                continue
            
            # Try to get constant value
            try:
                const_value = float(predicted) if predicted and predicted != '' else 0.0
            except (ValueError, TypeError):
                const_value = 0.0
            
            try:
                error_val = float(error_pct) if error_pct and error_pct != '' else 0.0
            except (ValueError, TypeError):
                error_val = 0.0
            
            # Skip if no constant name or equation
            if not equation_latex:
                equation_latex = f"{{\\text{{{constant_name}}}}}"
            
            session_num = extract_session_heuristic(constant_name, label)
            domain = map_category_to_domain(category)
            
            derivation = {
                'constant_name': constant_name,
                'constant_value': const_value,
                'session_number': session_num,
                'domain': domain,
                'equation_latex': equation_latex,
                'description': label if label else f"Calibrated constant for {constant_name}",
                'status': status,
                'error_pct': error_val,
                'category': category,
            }
            
            derivations.append(derivation)
    
    print(f"[OK] Loaded {len(derivations)} derivations")
    return derivations


def extract_session_heuristic(constant_name: str, label: str) -> int:
    """Extract session number via heuristic."""
    for text in [constant_name, label]:
        match = re.search(r'[Ss]ession\s*(\d+)|S(\d{2,3})', text)
        if match:
            sess = match.group(1) or match.group(2)
            if sess:
                return int(sess)
    
    # Default ranges
    if constant_name.startswith(('F_', 'PHI_', 'SSQ', 'N_')):
        return 237
    elif constant_name.startswith(('RHO_VAC', 'BETA_', 'KAPPA')):
        return 250
    elif constant_name.startswith(('OMEGA_', 'TCMB', 'SIGMA')):
        return 371
    elif constant_name.startswith(('M_', 'ALPHA_')):
        return 286
    else:
        return 400


def map_category_to_domain(category: str) -> str:
    """Map CSV category to human-readable domain."""
    mapping = {
        'DERIVATION_FIRST_PRINCIPLES': 'First Principles',
        'CALIBRATION_FROM_PARAMETERS': 'Calibration & Tuning',
        'ORPHAN_DUPLICATE': 'Legacy & Cross-Validation',
        'RESEARCH_TRACE': 'Research Notes',
        'OPEN_PREDICTION': 'Open Predictions',
    }
    return mapping.get(category, 'Miscellaneous')


def group_by_domain(derivations: List[Dict]) -> Dict[str, List[Dict]]:
    """Group derivations by domain for organized insertion."""
    grouped = {}
    for deriv in derivations:
        domain = deriv['domain']
        if domain not in grouped:
            grouped[domain] = []
        grouped[domain].append(deriv)
    return grouped


def generate_derivation_code(deriv: Dict) -> str:
    """Generate Python code for a single derivation equation."""
    const_name = deriv['constant_name']
    const_value = deriv['constant_value']
    session = deriv['session_number']
    domain = deriv['domain']
    equation = deriv['equation_latex'].replace('"', r'\"')
    description = deriv['description'].replace('"', r'\"')
    status = deriv['status']
    error_pct = deriv['error_pct']
    
    code = f'''
    @staticmethod
    def {const_name.lower()}_derivation() -> DerivationEquation:
        """Constant: {const_name} (Session {session})."""
        return DerivationEquation(
            constant_name="{const_name}",
            constant_value={const_value},
            session_number={session},
            domain="{domain}",
            equation_latex=r"{equation}",
            description="{description}",
            derivation_steps=[
                "Extracted from master_closures.csv (calibrated constant)"
            ],
            source_paper=f"Session {{session}} data",
            status="{status}",
            error_pct={error_pct},
        )
'''
    return code


def generate_registry_entries(derivations: List[Dict]) -> str:
    """Generate registry dictionary entries for DerivationRegistry."""
    registry_code = ""
    
    # Group by session range
    ranges = {
        'SESSIONS_203_340': [],
        'SESSIONS_341_439': [],
        'SESSIONS_440_539': [],
        'SESSIONS_540_665': [],
        'SESSIONS_725_785': [],
    }
    
    for deriv in derivations:
        session = deriv['session_number']
        const_name = deriv['constant_name']
        
        if 203 <= session <= 340:
            ranges['SESSIONS_203_340'].append(const_name)
        elif 341 <= session <= 439:
            ranges['SESSIONS_341_439'].append(const_name)
        elif 440 <= session <= 539:
            ranges['SESSIONS_440_539'].append(const_name)
        elif 540 <= session <= 665:
            ranges['SESSIONS_540_665'].append(const_name)
        elif 725 <= session <= 785:
            ranges['SESSIONS_725_785'].append(const_name)
    
    for range_name, consts in ranges.items():
        if consts:
            registry_code += f"\n# Add to {range_name}:\n"
            for const in consts:
                registry_code += f'    "{const}": {const.lower()}_derivation(),\n'
    
    return registry_code


def print_merge_preview(grouped: Dict[str, List[Dict]]) -> None:
    """Print summary preview of merge operations."""
    print("\n" + "="*80)
    print("DERIVATION MERGE PREVIEW")
    print("="*80 + "\n")
    
    for domain, derivs in sorted(grouped.items()):
        print(f"\n{domain} ({len(derivs)} entries):")
        for deriv in derivs[:5]:  # Show first 5
            print(f"  • {deriv['constant_name']:<30} S{deriv['session_number']:03d}  "
                  f"Status: {deriv['status']:<15} Error: {deriv['error_pct']:.4f}%")
        if len(derivs) > 5:
            print(f"  ... and {len(derivs)-5} more")
    
    total = sum(len(v) for v in grouped.values())
    print(f"\n{'='*80}")
    print(f"TOTAL DERIVATIONS TO MERGE: {total}")
    print(f"{'='*80}\n")


def main():
    print("\n" + "="*80)
    print("TRINITY DERIVATION MERGE SCRIPT")
    print("="*80 + "\n")
    
    # Load derivations
    derivations = load_extracted_derivations()
    
    # Group by domain
    grouped = group_by_domain(derivations)
    
    # Show preview
    print_merge_preview(grouped)
    
    # Generate code snippets
    print("\nGenerating code snippets for domain classes...")
    print("\nSample Code (First Derivation in Each Domain):\n")
    
    for domain in sorted(grouped.keys()):
        deriv = grouped[domain][0]
        print(f"\n--- {domain} ---")
        print(generate_derivation_code(deriv))
    
    # Generate registry entries
    print("\n\nRegistry Insertion Code:\n")
    print(generate_registry_entries(derivations))
    
    print("\n" + "="*80)
    print("MERGE COMPLETE — Ready to insert into _uqff_derivation_equations.py")
    print("="*80 + "\n")


if __name__ == '__main__':
    main()
