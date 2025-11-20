#!/usr/bin/env python3
"""
Comprehensive Physics Class Extractor
Scans source1.cpp through source176.cpp for ALL physics modules
"""

import re
import os
import csv
from pathlib import Path

# Excluded base classes and utilities
EXCLUDED_CLASSES = {
    'PhysicsTerm', 'UQFFModule', 'ModuleInterface',
    'QWidget', 'QApplication', 'QDialog', 'QMainWindow', 'QGraphicsItem',
    'QAbstractListModel', 'QSyntaxHighlighter', 'QPushButton', 'QUndoCommand',
    'BaseErrorListener', 'MathBaseVisitor', 'antlr4',
    'SymEngineAllocator', 'Units', 'PerlinNoise', 'Shader', 'Camera', 'Bone', 'SIMPlugin',
    'MathErrorListener', 'SymEngineVisitor', 'VarCollectorVisitor', 'MathHighlighter',
    'DraggableButton', 'InsertCommand', 'MacroCommand', 'ControlPointItem',
    'EquationSuggestModel', 'ScientificCalculatorDialog', 'FluidSolver'
}

# Physics type keywords for classification
PHYSICS_KEYWORDS = {
    'nuclear': ['nuclear', 'hydrogen', 'helium', 'fusion', 'fission', 'isotope', 'atomic', 'proton', 'neutron'],
    'gravity': ['gravity', 'gravitational', 'compressed', 'compression', 'mass', 'blackhole', 'smbh'],
    'magnetic': ['magnetic', 'field', 'magnetar', 'dipole', 'lorentz'],
    'vacuum': ['vacuum', 'aether', 'quantum', 'coupling', 'dynamic'],
    'resonance': ['resonance', 'frequency', 'harmonic', 'oscillation'],
    'galactic': ['galaxy', 'ngc', 'andromeda', 'm51', 'spiral', 'elliptical'],
    'stellar': ['star', 'stellar', 'nebula', 'orion', 'carina', 'vela', 'dwarf', 'supernova'],
    'cosmological': ['universe', 'cosmological', 'bigbang', 'expansion', 'hubble'],
    'unified': ['uqff', 'unified', 'field', 'multi'],
    'other': []
}

def classify_physics_type(class_name, file_content):
    """Determine physics category from class name and context"""
    class_lower = class_name.lower()
    
    for phys_type, keywords in PHYSICS_KEYWORDS.items():
        if phys_type == 'other':
            continue
        for keyword in keywords:
            if keyword in class_lower:
                return phys_type
    
    # Check context around class definition
    context_lower = file_content.lower()
    for phys_type, keywords in PHYSICS_KEYWORDS.items():
        if phys_type == 'other':
            continue
        score = sum(1 for kw in keywords if kw in context_lower)
        if score >= 2:
            return phys_type
    
    return 'other'

def extract_method_signature(file_content, class_start_pos):
    """Extract compute/calculate/solve method signature if exists"""
    # Search for methods after class definition
    class_section = file_content[class_start_pos:class_start_pos + 5000]  # Check next 5000 chars
    
    method_patterns = [
        r'(\w+)\s+compute\w*\s*\([^)]*\)',
        r'(\w+)\s+calculate\w*\s*\([^)]*\)',
        r'(\w+)\s+solve\w*\s*\([^)]*\)',
        r'(\w+)\s+evaluate\w*\s*\([^)]*\)',
        r'(\w+)\s+process\w*\s*\([^)]*\)'
    ]
    
    for pattern in method_patterns:
        match = re.search(pattern, class_section, re.IGNORECASE)
        if match:
            return match.group(0).strip()
    
    return ""

def scan_source_file(filepath):
    """Extract all physics classes from a single source file"""
    results = []
    
    try:
        with open(filepath, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        # Find all class declarations
        # Pattern: class ClassName or class ClassName : public BaseClass
        class_pattern = r'class\s+(\w+)(?:\s*:\s*public\s+(\w+))?\s*[{:]'
        
        for match in re.finditer(class_pattern, content, re.MULTILINE):
            class_name = match.group(1)
            base_class = match.group(2) if match.group(2) else "standalone"
            
            # Skip excluded classes
            if class_name in EXCLUDED_CLASSES:
                continue
            
            # Calculate line number
            line_num = content[:match.start()].count('\n') + 1
            
            # Classify physics type
            physics_type = classify_physics_type(class_name, content)
            
            # Extract method signature
            method_sig = extract_method_signature(content, match.start())
            
            results.append({
                'SourceFile': os.path.basename(filepath),
                'ClassName': class_name,
                'LineNumber': line_num,
                'BaseClass': base_class,
                'PhysicsType': physics_type,
                'MethodSignature': method_sig
            })
    
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
    
    return results

def main():
    """Main extraction process"""
    all_classes = []
    
    # Scan source1.cpp through source176.cpp
    for i in range(1, 177):
        filepath = f"source{i}.cpp"
        
        if os.path.exists(filepath):
            print(f"Scanning {filepath}...", end=' ')
            classes = scan_source_file(filepath)
            all_classes.extend(classes)
            print(f"Found {len(classes)} classes")
        else:
            print(f"Skipping {filepath} (not found)")
    
    # Remove duplicates (same class name in same file)
    unique_classes = []
    seen = set()
    for cls in all_classes:
        key = (cls['SourceFile'], cls['ClassName'])
        if key not in seen:
            seen.add(key)
            unique_classes.append(cls)
    
    print(f"\n{'='*80}")
    print(f"TOTAL UNIQUE PHYSICS CLASSES FOUND: {len(unique_classes)}")
    print(f"{'='*80}\n")
    
    # Write to CSV
    output_file = 'COMPLETE_PHYSICS_CLASS_INVENTORY.csv'
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        fieldnames = ['SourceFile', 'ClassName', 'LineNumber', 'BaseClass', 'PhysicsType', 'MethodSignature']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        writer.writeheader()
        for cls in sorted(unique_classes, key=lambda x: (int(re.search(r'\d+', x['SourceFile']).group()), x['LineNumber'])):
            writer.writerow(cls)
    
    print(f"Complete inventory written to: {output_file}")
    
    # Print summary statistics
    print(f"\nPhysics Type Distribution:")
    type_counts = {}
    for cls in unique_classes:
        ptype = cls['PhysicsType']
        type_counts[ptype] = type_counts.get(ptype, 0) + 1
    
    for ptype, count in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
        print(f"  {ptype:15s}: {count:4d} classes")
    
    # Print first 20 examples
    print(f"\nFirst 20 Classes (sample):")
    print(f"{'File':<15s} {'Class':<40s} {'Line':<6s} {'Base':<20s} {'Type':<15s}")
    print(f"{'-'*110}")
    for cls in unique_classes[:20]:
        print(f"{cls['SourceFile']:<15s} {cls['ClassName']:<40s} {str(cls['LineNumber']):<6s} {cls['BaseClass']:<20s} {cls['PhysicsType']:<15s}")

if __name__ == '__main__':
    main()
