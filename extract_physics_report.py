#!/usr/bin/env python3
"""
Physics Term Extraction Report Generator
Scans all source*.cpp files for classes, functions, constants, WOLFRAM_TERMs
Target: 3000+ unique physics terms for integration into MAIN_1_CoAnQi.cpp
"""

import re
import glob
import csv
import json

def main():
    files = sorted(glob.glob('source*.cpp'))
    results = []
    
    print("Scanning all source files for physics patterns...")
    
    for f in files:
        try:
            with open(f, 'r', encoding='utf-8', errors='ignore') as file:
                content = file.read()
                
                # Extract classes/structs
                classes = re.findall(r'^(?:class|struct)\s+(\w+)', content, re.MULTILINE)
                
                # Extract functions (double, void, std::vector, etc.)
                functions = re.findall(r'^\s*(?:double|void|int|std::vector|std::map|std::string|bool)\s+(\w+)\s*\(', content, re.MULTILINE)
                
                # Extract constants
                constants = re.findall(r'const(?:expr)?\s+(?:double|int)\s+(\w+)\s*=', content)
                
                # Extract WOLFRAM_TERM
                wolfram = re.findall(r'#define\s+WOLFRAM_TERM', content)
                
                results.append({
                    'file': f,
                    'classes': len(classes),
                    'functions': len(functions),
                    'constants': len(constants),
                    'wolfram_terms': len(wolfram),
                    'total': len(classes) + len(functions) + len(constants) + len(wolfram),
                    'sample_classes': ', '.join(classes[:5]) if classes else '',
                    'sample_functions': ', '.join(functions[:5]) if functions else '',
                    'sample_constants': ', '.join(constants[:5]) if constants else '',
                    'class_list': classes,
                    'function_list': functions,
                    'constant_list': constants
                })
        except Exception as e:
            print(f'Error reading {f}: {e}')
    
    # Write to CSV
    with open('physics_extraction_report.csv', 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['file', 'classes', 'functions', 'constants', 'wolfram_terms', 'total', 
                      'sample_classes', 'sample_functions', 'sample_constants']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = {k: v for k, v in r.items() if k in fieldnames}
            writer.writerow(row)
    
    # Calculate totals
    total_classes = sum(r['classes'] for r in results)
    total_functions = sum(r['functions'] for r in results)
    total_constants = sum(r['constants'] for r in results)
    total_wolfram = sum(r['wolfram_terms'] for r in results)
    grand_total = sum(r['total'] for r in results)
    
    # Generate summary report
    summary = {
        'total_files': len(results),
        'total_patterns': grand_total,
        'classes': total_classes,
        'functions': total_functions,
        'constants': total_constants,
        'wolfram_terms': total_wolfram,
        'target_met': grand_total >= 3000,
        'files_analyzed': [r['file'] for r in results]
    }
    
    with open('physics_extraction_summary.json', 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=2)
    
    # Print summary
    print(f'\n{"="*70}')
    print(f'PHYSICS EXTRACTION COMPLETE - {len(results)} FILES ANALYZED')
    print(f'{"="*70}')
    print(f'TOTAL PHYSICS PATTERNS FOUND: {grand_total:,}')
    print(f'  - Classes/Structs: {total_classes:,}')
    print(f'  - Functions: {total_functions:,}')
    print(f'  - Constants: {total_constants:,}')
    print(f'  - WOLFRAM_TERM macros: {total_wolfram:,}')
    print(f'{"="*70}')
    print(f'\nDetailed breakdown saved to: physics_extraction_report.csv')
    print(f'Summary JSON saved to: physics_extraction_summary.json')
    print(f'\nTARGET GOAL: 3000+ terms')
    print(f'STATUS: {"✅ EXCEEDED" if grand_total >= 3000 else "❌ NOT MET"} ({grand_total:,} / 3000)')
    print(f'\n{"="*70}')
    
    # Show top 10 files by term count
    print(f'\nTOP 10 FILES BY PHYSICS TERM COUNT:')
    sorted_results = sorted(results, key=lambda x: x['total'], reverse=True)
    for i, r in enumerate(sorted_results[:10], 1):
        print(f'  {i:2}. {r["file"]:25} - {r["total"]:4} terms (C:{r["classes"]:3} F:{r["functions"]:4} K:{r["constants"]:3})')
    
    return summary

if __name__ == '__main__':
    main()
