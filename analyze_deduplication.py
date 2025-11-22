#!/usr/bin/env python3
"""
Deduplication Analysis: Extracted Patterns vs Existing MAIN_1_CoAnQi.cpp
Compares 4,890 extracted patterns against current MAIN_1_CoAnQi.cpp implementation
"""

import re
import json
import glob
from pathlib import Path
from collections import defaultdict

def extract_patterns_from_main(main_file_path):
    """Extract all classes, functions, constants from MAIN_1_CoAnQi.cpp"""
    
    with open(main_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    # Extract patterns
    existing = {
        'classes': set(),
        'functions': set(),
        'constants': set()
    }
    
    # Classes and structs
    class_pattern = r'^(?:class|struct)\s+(\w+)'
    for match in re.finditer(class_pattern, content, re.MULTILINE):
        existing['classes'].add(match.group(1))
    
    # Functions (with return types)
    func_pattern = r'^\s*(?:double|void|int|std::vector|std::map|std::string|bool)\s+(\w+)\s*\('
    for match in re.finditer(func_pattern, content, re.MULTILINE):
        existing['functions'].add(match.group(1))
    
    # Constants
    const_pattern = r'const(?:expr)?\s+(?:double|int)\s+(\w+)\s*='
    for match in re.finditer(const_pattern, content, re.MULTILINE):
        existing['constants'].add(match.group(1))
    
    return existing

def extract_patterns_from_sources(source_dir="."):
    """Extract patterns directly from all source*.cpp files"""
    
    extracted = {
        'classes': defaultdict(list),
        'functions': defaultdict(list),
        'constants': defaultdict(list)
    }
    
    # Find all source files
    source_files = []
    for pattern in ['source*.cpp', 'Source*.cpp']:
        source_files.extend(glob.glob(str(Path(source_dir) / pattern)))
    
    print(f"Scanning {len(source_files)} source files...")
    
    for filepath in source_files:
        filename = Path(filepath).name
        
        try:
            with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
        except Exception as e:
            print(f"  Warning: Could not read {filename}: {e}")
            continue
        
        # Extract classes
        class_pattern = r'^(?:class|struct)\s+(\w+)'
        for match in re.finditer(class_pattern, content, re.MULTILINE):
            cls_name = match.group(1)
            extracted['classes'][cls_name].append(filename)
        
        # Extract functions
        func_pattern = r'^\s*(?:double|void|int|std::vector|std::map|std::string|bool)\s+(\w+)\s*\('
        for match in re.finditer(func_pattern, content, re.MULTILINE):
            func_name = match.group(1)
            extracted['functions'][func_name].append(filename)
        
        # Extract constants
        const_pattern = r'const(?:expr)?\s+(?:double|int)\s+(\w+)\s*='
        for match in re.finditer(const_pattern, content, re.MULTILINE):
            const_name = match.group(1)
            extracted['constants'][const_name].append(filename)
    
    return extracted

def compare_patterns(existing, extracted):
    """Compare extracted patterns against existing MAIN_1_CoAnQi.cpp"""
    
    results = {
        'classes': {
            'already_integrated': {},
            'net_new': {},
            'total_extracted': len(extracted['classes']),
            'total_existing': len(existing['classes'])
        },
        'functions': {
            'already_integrated': {},
            'net_new': {},
            'total_extracted': len(extracted['functions']),
            'total_existing': len(existing['functions'])
        },
        'constants': {
            'already_integrated': {},
            'net_new': {},
            'total_extracted': len(extracted['constants']),
            'total_existing': len(existing['constants'])
        }
    }
    
    # Analyze classes
    for cls, sources in extracted['classes'].items():
        if cls in existing['classes']:
            results['classes']['already_integrated'][cls] = sources
        else:
            results['classes']['net_new'][cls] = sources
    
    # Analyze functions
    for func, sources in extracted['functions'].items():
        if func in existing['functions']:
            results['functions']['already_integrated'][func] = sources
        else:
            results['functions']['net_new'][func] = sources
    
    # Analyze constants
    for const, sources in extracted['constants'].items():
        if const in existing['constants']:
            results['constants']['already_integrated'][const] = sources
        else:
            results['constants']['net_new'][const] = sources
    
    return results

def generate_report(results):
    """Generate comprehensive deduplication report"""
    
    print("\n" + "="*100)
    print("DEDUPLICATION ANALYSIS REPORT")
    print("Comparing 4,890 extracted patterns vs existing MAIN_1_CoAnQi.cpp")
    print("="*100 + "\n")
    
    # Summary statistics
    total_already = (len(results['classes']['already_integrated']) + 
                    len(results['functions']['already_integrated']) + 
                    len(results['constants']['already_integrated']))
    
    total_net_new = (len(results['classes']['net_new']) + 
                    len(results['functions']['net_new']) + 
                    len(results['constants']['net_new']))
    
    total_extracted = (results['classes']['total_extracted'] + 
                      results['functions']['total_extracted'] + 
                      results['constants']['total_extracted'])
    
    print(f"EXECUTIVE SUMMARY:")
    print(f"  Total Extracted Patterns:        {total_extracted:,}")
    print(f"  Already Integrated (duplicates):  {total_already:,} ({100*total_already/total_extracted:.1f}%)")
    print(f"  Net-New Unique Patterns:          {total_net_new:,} ({100*total_net_new/total_extracted:.1f}%)")
    print()
    
    # Category breakdown
    print("CATEGORY BREAKDOWN:")
    print("-" * 100)
    
    categories = ['classes', 'functions', 'constants']
    for cat in categories:
        cat_name = cat.upper()
        total_ext = results[cat]['total_extracted']
        total_exist = results[cat]['total_existing']
        already = len(results[cat]['already_integrated'])
        net_new = len(results[cat]['net_new'])
        
        print(f"\n{cat_name}:")
        print(f"  Extracted from source files:        {total_ext:,}")
        print(f"  Already in MAIN_1_CoAnQi.cpp:       {total_exist:,}")
        print(f"  Duplicates (already integrated):    {already:,} ({100*already/total_ext:.1f}%)")
        print(f"  Net-New (unique to add):            {net_new:,} ({100*net_new/total_ext:.1f}%)")
    
    print("\n" + "="*100)
    print("NET-NEW PATTERNS TO INTEGRATE (Top 20 per category)")
    print("="*100 + "\n")
    
    # Show samples of net-new patterns
    for cat in categories:
        print(f"\n{cat.upper()} (Top 20 net-new):")
        print("-" * 100)
        net_new_items = list(results[cat]['net_new'].items())[:20]
        
        if not net_new_items:
            print("  (All patterns already integrated)")
        else:
            for name, sources in net_new_items:
                source_list = ', '.join(sources[:3])
                if len(sources) > 3:
                    source_list += f" (+{len(sources)-3} more)"
                print(f"  • {name:<40} from: {source_list}")
    
    print("\n" + "="*100)
    print("ALREADY INTEGRATED SAMPLES (Top 10 per category)")
    print("="*100 + "\n")
    
    # Show samples of duplicates
    for cat in categories:
        print(f"\n{cat.upper()} (Top 10 already integrated):")
        print("-" * 100)
        already_items = list(results[cat]['already_integrated'].items())[:10]
        
        if not already_items:
            print("  (None found)")
        else:
            for name, sources in already_items:
                source_list = ', '.join(sources[:3])
                if len(sources) > 3:
                    source_list += f" (+{len(sources)-3} more)"
                print(f"  • {name:<40} from: {source_list}")
    
    # Save detailed results
    output_file = "deduplication_analysis_results.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print("\n" + "="*100)
    print(f"DETAILED RESULTS SAVED TO: {output_file}")
    print("="*100 + "\n")
    
    # Recommendation
    print("RECOMMENDATION:")
    print("-" * 100)
    if total_net_new > 0:
        print(f"✓ Proceed with Option 2 (Direct In-File Integration)")
        print(f"✓ Integrate {total_net_new:,} net-new unique patterns only")
        print(f"✓ Skip {total_already:,} duplicates already in MAIN_1_CoAnQi.cpp")
        print(f"✓ Expected file size increase: ~{total_net_new//20:,} lines")
    else:
        print("⚠ All extracted patterns already integrated!")
        print("⚠ No new physics to add - MAIN_1_CoAnQi.cpp is complete")
    print()

def main():
    # File paths
    main_file = Path("MAIN_1_CoAnQi.cpp")
    
    if not main_file.exists():
        print(f"ERROR: {main_file} not found!")
        return
    
    print("="*100)
    print("STEP 1: Analyzing existing MAIN_1_CoAnQi.cpp patterns...")
    print("="*100)
    existing = extract_patterns_from_main(main_file)
    
    print(f"\nFound in MAIN_1_CoAnQi.cpp ({main_file.stat().st_size//1024} KB):")
    print(f"  Classes:   {len(existing['classes']):,}")
    print(f"  Functions: {len(existing['functions']):,}")
    print(f"  Constants: {len(existing['constants']):,}")
    print(f"  TOTAL:     {sum(len(v) for v in existing.values()):,}")
    print()
    
    print("="*100)
    print("STEP 2: Extracting patterns from all source*.cpp files...")
    print("="*100)
    extracted = extract_patterns_from_sources()
    
    print(f"\nExtracted from all source files:")
    print(f"  Classes:   {len(extracted['classes']):,}")
    print(f"  Functions: {len(extracted['functions']):,}")
    print(f"  Constants: {len(extracted['constants']):,}")
    print(f"  TOTAL:     {sum(len(v) for v in extracted.values()):,}")
    print()
    
    print("="*100)
    print("STEP 3: Comparing patterns for duplication...")
    print("="*100)
    results = compare_patterns(existing, extracted)
    
    generate_report(results)

if __name__ == "__main__":
    main()
