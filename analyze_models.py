#!/usr/bin/env python3
"""Comprehensive analysis of CondensedPhysics.py models and resources."""

import re
from collections import defaultdict

def analyze():
    with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
        lines = f.readlines()

    print('='*100)
    print('6. MODEL UPGRADE STATUS (compute_clean_equation)')
    print('='*100)

    # Find all models
    models = []
    for i, line in enumerate(lines):
        if re.match(r'^class \w+Model', line):
            match = re.search(r'class (\w+Model)', line)
            if match:
                models.append((i+1, match.group(1)))

    # Check upgrade status for each model
    upgraded = []
    needs_upgrade = []

    for start_line, name in models:
        # Find end of model
        end_line = len(lines)
        for j in range(start_line, len(lines)):
            if j > start_line and lines[j].startswith('class '):
                end_line = j
                break
        
        model_lines = end_line - start_line
        
        # Check for upgrade indicators
        has_clean_eq = any('def compute_clean_equation' in lines[i] for i in range(start_line-1, end_line))
        has_resonance = any('def compute_resonance' in lines[i] for i in range(start_line-1, end_line))
        has_validate = any('def validate_model' in lines[i] for i in range(start_line-1, end_line))
        has_run_tests = any('def run_tests' in lines[i] for i in range(start_line-1, end_line))
        
        record = {
            'name': name,
            'line': start_line,
            'size': model_lines,
            'clean_eq': has_clean_eq,
            'resonance': has_resonance,
            'validate': has_validate,
            'run_tests': has_run_tests
        }
        
        if has_clean_eq:
            upgraded.append(record)
        else:
            needs_upgrade.append(record)

    print(f'UPGRADED MODELS (have compute_clean_equation): {len(upgraded)}')
    print('-'*100)
    for m in upgraded:
        print(f"  {m['line']:5} | {m['name']:40} | {m['size']:4} lines | CleanEq:True  | Res:{str(m['resonance']):5} | Val:{str(m['validate']):5} | Test:{str(m['run_tests']):5}")

    print()
    print(f'MODELS NEEDING UPGRADE: {len(needs_upgrade)}')
    print('-'*100)
    for m in needs_upgrade:
        print(f"  {m['line']:5} | {m['name']:40} | {m['size']:4} lines | CleanEq:False | Res:{str(m['resonance']):5} | Val:{str(m['validate']):5} | Test:{str(m['run_tests']):5}")

    print()
    print('='*100)
    print('7. FUNCTIONS/METHODS COUNT PER MODEL')
    print('='*100)
    
    # Count methods per model
    method_counts = []
    for start_line, name in models:
        end_line = len(lines)
        for j in range(start_line, len(lines)):
            if j > start_line and lines[j].startswith('class '):
                end_line = j
                break
        
        method_count = sum(1 for i in range(start_line-1, end_line) if '    def ' in lines[i])
        method_counts.append((name, method_count, start_line))
    
    # Sort by method count descending
    method_counts.sort(key=lambda x: x[1], reverse=True)
    
    print('Top 20 models by method count:')
    for name, count, line in method_counts[:20]:
        print(f"  {count:3} methods | {name:40} (line {line})")
    
    print()
    print('='*100)
    print('8. DUPLICATE FUNCTION NAMES CHECK')
    print('='*100)
    
    # Find all function definitions
    all_funcs = defaultdict(list)
    for i, line in enumerate(lines):
        if re.match(r'^def \w+', line):
            match = re.search(r'def (\w+)', line)
            if match:
                all_funcs[match.group(1)].append(i+1)
    
    # Find duplicates (more than 1 occurrence)
    duplicates = {name: locs for name, locs in all_funcs.items() if len(locs) > 1}
    
    if duplicates:
        print(f'DUPLICATE TOP-LEVEL FUNCTIONS FOUND: {len(duplicates)}')
        for name, locs in duplicates.items():
            print(f"  {name}: lines {locs}")
    else:
        print('No duplicate top-level function names found.')
    
    print()
    print('='*100)
    print('9. CONSTANTS DICTIONARY SIZE')
    print('='*100)
    
    # Find CONSTANTS dict and count entries
    in_constants = False
    constants_start = 0
    constants_end = 0
    constants_count = 0
    
    for i, line in enumerate(lines):
        if line.strip().startswith('CONSTANTS = {'):
            in_constants = True
            constants_start = i + 1
        elif in_constants:
            if line.strip() == '}':
                in_constants = False
                constants_end = i + 1
                break
            if ':' in line and not line.strip().startswith('#'):
                constants_count += 1
    
    print(f'CONSTANTS dictionary: lines {constants_start}-{constants_end}')
    print(f'Total entries: ~{constants_count} key-value pairs')
    
    print()
    print('='*100)
    print('10. MAY 2025 MODELS SPECIFIC STATUS')
    print('='*100)
    
    may_2025_models = [
        'NGC2264Model', 'UGC10214Model', 'NGC4676Model', 'RedSpiderNebulaModel',
        'NGC3372Model', 'AGCarinaeModel', 'M42Model', 'TarantulaNebulaModel',
        'NGC2841Model', 'MysticMountainModel'
    ]
    
    for m in needs_upgrade:
        if m['name'] in may_2025_models:
            print(f"  {m['line']:5} | {m['name']:25} | {m['size']:4} lines | NEEDS: compute_clean_equation")
    
    print()
    print('='*100)
    print('SUMMARY')
    print('='*100)
    print(f'Total Lines: {len(lines)}')
    print(f'Total Models: {len(models)}')
    print(f'Upgraded Models (have compute_clean_equation): {len(upgraded)}')
    print(f'Models Needing Upgrade: {len(needs_upgrade)}')
    print(f'May 2025 Models Needing compute_clean_equation: 10')

if __name__ == '__main__':
    analyze()
