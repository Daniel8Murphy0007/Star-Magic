#!/usr/bin/env python3
"""
Wolfram Full Extraction - Phase 1: Source File Deep Scan
Scans all source1-178.cpp files to extract physics terms, equations, systems, and self-expanding methods
"""

import re
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Physics patterns to extract
PATTERNS = {
    'constants': [
        r'const\s+double\s+([A-Z_][A-Z0-9_]*)\s*=\s*([\d\.eE\+\-]+)',  # const double G = 6.674e-11
        r'#define\s+([A-Z_][A-Z0-9_]*)\s+([\d\.eE\+\-]+)',              # #define PI 3.14159
        r'static\s+constexpr\s+double\s+([A-Z_a-z][A-Za-z0-9_]*)\s*=\s*([\d\.eE\+\-]+)',
    ],
    'equations': [
        r'double\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\([^)]*\)\s*{[^}]*return\s+([^;]+);',  # Functions returning calculations
        r'([A-Z_][A-Z0-9_]*)\s*=\s*([^;]+);',  # Variable assignments with complex RHS
    ],
    'systems': [
        r'SystemParams\s+([a-zA-Z_][a-zA-Z0-9_]*)\s*\{',  # SystemParams declarations
        r'systems\["([^"]+)"\]\s*=\s*SystemParams',        # Map insertions
        r'struct\s+([A-Z][A-Za-z0-9_]*System)\s*\{',       # System structures
    ],
    'self_update': [
        r'setDynamicParameter\s*\(\s*"([^"]+)"\s*,\s*([^)]+)\)',
        r'getDynamicParameter\s*\(\s*"([^"]+)"\s*\)',
        r'registerDynamicTerm\s*\(\s*([^)]+)\)',
    ],
    'self_expansion': [
        r'class\s+([A-Z][A-Za-z0-9_]*)\s*:\s*public\s+Module2_0Enhanced',
        r'void\s+addEnhancement\s*\([^)]*\)',
        r'\.registerDynamicTerm\s*\(',
    ],
    'self_simulation': [
        r'void\s+simulate\s*\([^)]*\)',
        r'void\s+exportState\s*\([^)]*\)',
        r'for\s*\([^;]*time[^;]*;[^)]*\)',  # Time evolution loops
    ],
}

def extract_from_file(file_path):
    """Extract all physics entities from a single source file"""
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
    except Exception as e:
        return {'error': str(e)}
    
    results = defaultdict(list)
    
    # Extract constants
    for pattern in PATTERNS['constants']:
        matches = re.findall(pattern, content)
        for name, value in matches:
            results['constants'].append({'name': name, 'value': value, 'line': content[:content.find(name)].count('\n') + 1})
    
    # Extract equations (function names and expressions)
    for pattern in PATTERNS['equations']:
        matches = re.findall(pattern, content, re.DOTALL)
        for name, expr in matches:
            if len(expr) > 10 and len(expr) < 500:  # Filter reasonable expressions
                results['equations'].append({'name': name, 'expression': expr.strip()[:200]})
    
    # Extract systems
    for pattern in PATTERNS['systems']:
        matches = re.findall(pattern, content)
        for match in matches:
            system_name = match if isinstance(match, str) else match[0]
            results['systems'].append({'name': system_name})
    
    # Extract self-update methods
    for pattern in PATTERNS['self_update']:
        matches = re.findall(pattern, content)
        for match in matches:
            param_name = match if isinstance(match, str) else match[0]
            results['self_update'].append({'method': 'setDynamicParameter', 'param': param_name})
    
    # Extract self-expansion indicators
    for pattern in PATTERNS['self_expansion']:
        if re.search(pattern, content):
            results['self_expansion'].append({'type': 'Module2_0Enhanced'})
    
    # Extract self-simulation methods
    for pattern in PATTERNS['self_simulation']:
        if re.search(pattern, content):
            results['self_simulation'].append({'type': 'time_evolution'})
    
    # Count unique entities
    summary = {
        'constants_count': len(results['constants']),
        'equations_count': len(results['equations']),
        'systems_count': len(results['systems']),
        'self_update_count': len(results['self_update']),
        'self_expansion': bool(results['self_expansion']),
        'self_simulation': bool(results['self_simulation']),
        'total_entities': len(results['constants']) + len(results['equations']) + len(results['systems'])
    }
    
    return {
        'summary': summary,
        'entities': dict(results)
    }

def scan_all_source_files():
    """Scan all source1-178.cpp files"""
    base_path = Path('.')
    results = {}
    total_stats = defaultdict(int)
    
    print("=" * 80)
    print("PHASE 1: COMPREHENSIVE SOURCE FILE ANALYSIS")
    print("=" * 80)
    print()
    
    # Scan source1.cpp through source178.cpp
    files_found = 0
    files_missing = []
    
    for i in range(1, 179):
        file_name = f'source{i}.cpp'
        file_path = base_path / file_name
        
        if file_path.exists():
            files_found += 1
            print(f"Scanning {file_name}...", end=' ')
            
            data = extract_from_file(file_path)
            results[file_name] = data
            
            if 'summary' in data:
                summary = data['summary']
                print(f"✓ {summary['total_entities']} entities ({summary['constants_count']}c, {summary['equations_count']}e, {summary['systems_count']}s)")
                
                total_stats['constants'] += summary['constants_count']
                total_stats['equations'] += summary['equations_count']
                total_stats['systems'] += summary['systems_count']
                total_stats['total'] += summary['total_entities']
            else:
                print(f"✗ Error: {data.get('error', 'Unknown')}")
        else:
            files_missing.append(file_name)
    
    print()
    print("=" * 80)
    print("SCAN SUMMARY")
    print("=" * 80)
    print(f"Files found: {files_found}/178")
    print(f"Files missing: {len(files_missing)}")
    print(f"Total constants: {total_stats['constants']}")
    print(f"Total equations: {total_stats['equations']}")
    print(f"Total systems: {total_stats['systems']}")
    print(f"Total entities: {total_stats['total']}")
    print()
    
    # Save results to JSON
    output_file = 'wolfram_extraction/source_inventory_complete.json'
    Path('wolfram_extraction').mkdir(exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        json.dump({
            'scan_date': datetime.now().isoformat(),
            'files_scanned': files_found,
            'files_missing': files_missing,
            'total_statistics': dict(total_stats),
            'file_results': results
        }, f, indent=2)
    
    print(f"✓ Results saved to: {output_file}")
    print()
    
    return results, total_stats

if __name__ == '__main__':
    results, stats = scan_all_source_files()
    
    print("=" * 80)
    print("NEXT STEPS")
    print("=" * 80)
    print("1. Review source_inventory_complete.json for extracted entities")
    print("2. Run Phase 2: Wolfram Knowledge Validation")
    print("3. Generate companion sourceXXX_wolfram.cpp files")
    print()
