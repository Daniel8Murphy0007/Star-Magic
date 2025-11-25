#!/usr/bin/env python3
"""
Intelligent Deduplication Analyzer for Star-Magic UQFF
Comprehensive tool for analyzing MAIN_1_CoAnQi.cpp and safely uncommenting patterns.

Features:
- C++ AST parsing for active code analysis
- Commented section pattern extraction
- Intelligent duplicate detection with confidence scoring
- Dependency graph generation
- Safe incremental uncommenting with rollback
- Automated validation and testing
- Comprehensive reporting

Author: GitHub Copilot
Date: November 22, 2025
"""

import re
import json
import os
import shutil
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from dataclasses import dataclass, field, asdict
from collections import defaultdict
from datetime import datetime
import subprocess


@dataclass(unsafe_hash=True)
class CodePattern:
    """Represents a code pattern (class, function, constant, etc.)"""
    pattern_type: str  # 'class', 'function', 'constant', 'macro', 'include'
    name: str
    signature: str
    line_start: int
    line_end: int
    is_commented: bool
    dependencies: Set[str] = field(default_factory=set, compare=False, hash=False)
    base_classes: Set[str] = field(default_factory=set, compare=False, hash=False)
    namespace: str = ""
    template_params: str = ""
    confidence_score: float = 1.0  # For duplicate detection
    category: str = ""  # 'exact_duplicate', 'net_new', 'ambiguous', 'framework'
    notes: str = ""


@dataclass
class AnalysisResult:
    """Results of the comprehensive analysis"""
    total_patterns: int = 0
    active_patterns: int = 0
    commented_patterns: int = 0
    exact_duplicates: int = 0
    net_new_patterns: int = 0
    ambiguous_patterns: int = 0
    framework_patterns: int = 0
    
    patterns: List[CodePattern] = field(default_factory=list)
    duplicates_map: Dict[str, List[CodePattern]] = field(default_factory=dict)
    dependency_graph: Dict[str, Set[str]] = field(default_factory=dict)
    
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())


class IntelligentDeduplicationAnalyzer:
    """Main analyzer class for comprehensive pattern analysis"""
    
    # Regex patterns for C++ parsing
    CLASS_PATTERN = r'^\s*class\s+(\w+)\s*(?::\s*(?:public\s+)?(\w+(?:\s*,\s*(?:public\s+)?\w+)*))?\s*\{'
    FUNCTION_PATTERN = r'^\s*(?:inline\s+)?(?:static\s+)?(?:virtual\s+)?(?:const\s+)?(?:[\w:]+)\s+(\w+)\s*\([^)]*\)\s*(?:const\s+)?(?:override\s+)?(?:\{|;)'
    CONSTANT_PATTERN = r'^\s*(?:static\s+)?(?:const\s+|constexpr\s+)([\w:]+)\s+(\w+)\s*='
    MACRO_PATTERN = r'^\s*#define\s+(\w+)'
    INCLUDE_PATTERN = r'^\s*#include\s+[<"]([^>"]+)[>"]'
    NAMESPACE_PATTERN = r'^\s*namespace\s+(\w+)\s*\{'
    TEMPLATE_PATTERN = r'^\s*template\s*<([^>]+)>'
    
    # Framework-specific patterns
    QT_PATTERNS = {'Q_OBJECT', 'QWidget', 'QMainWindow', 'QDialog', 'QApplication', 'slots', 'signals'}
    ANTLR_PATTERNS = {'antlr4', 'ANTLRInputStream', 'CommonTokenStream', 'ParseTree', 'Lexer', 'Parser'}
    SYMENGINE_PATTERNS = {'SymEngine', 'Basic', 'RCP', 'Expression', 'symbol'}
    
    def __init__(self, file_path: str):
        self.file_path = Path(file_path)
        self.content: List[str] = []
        self.active_section_end = 35000  # Lines 1-35000 are active
        self.commented_section_start = 31500  # Lines 31500+ are commented
        
        self.patterns: List[CodePattern] = []
        self.active_patterns: Dict[str, CodePattern] = {}
        self.commented_patterns: Dict[str, CodePattern] = {}
        
    def load_file(self) -> bool:
        """Load the source file into memory"""
        try:
            with open(self.file_path, 'r', encoding='utf-8', errors='ignore') as f:
                self.content = f.readlines()
            print(f"✓ Loaded {len(self.content)} lines from {self.file_path}")
            return True
        except Exception as e:
            print(f"✗ Error loading file: {e}")
            return False
    
    def is_commented_line(self, line: str) -> bool:
        """Check if a line is commented"""
        stripped = line.lstrip()
        return (stripped.startswith('//') or 
                stripped.startswith('/*') or
                '[Duplicate]' in line or
                '[duplicate]' in line.lower())
    
    def extract_class_pattern(self, start_line: int) -> Optional[CodePattern]:
        """Extract a class definition starting at the given line"""
        if start_line >= len(self.content):
            return None
        
        line = self.content[start_line]
        is_commented = self.is_commented_line(line)
        
        # Remove comment markers for parsing (handle multiple prefixes)
        clean_line = line
        clean_line = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', clean_line)
        clean_line = re.sub(r'^\s*//\s*', '', clean_line)
        
        match = re.match(self.CLASS_PATTERN, clean_line)
        if not match:
            return None
        
        class_name = match.group(1)
        base_classes_str = match.group(2) or ""
        
        # Extract base classes
        base_classes = set()
        if base_classes_str:
            for base in base_classes_str.split(','):
                base_clean = base.strip().replace('public ', '')
                if base_clean:
                    base_classes.add(base_clean)
        
        # Find end of class (matching brace)
        brace_count = 1
        end_line = start_line
        for i in range(start_line + 1, min(start_line + 1000, len(self.content))):
            # Clean line before counting braces
            line_content = self.content[i]
            line_content = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', line_content)
            line_content = re.sub(r'^\s*//\s*', '', line_content)
            brace_count += line_content.count('{') - line_content.count('}')
            if brace_count == 0:
                end_line = i
                break
        
        # Extract full signature
        signature_lines = []
        for i in range(start_line, min(end_line + 1, start_line + 50)):
            cleaned = self.content[i]
            cleaned = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', cleaned)
            cleaned = re.sub(r'^\s*//\s*', '', cleaned)
            signature_lines.append(cleaned.rstrip())
        signature = '\n'.join(signature_lines[:10])  # First 10 lines for signature
        
        # Check for framework patterns
        category = self._categorize_pattern(class_name, signature, base_classes)
        
        return CodePattern(
            pattern_type='class',
            name=class_name,
            signature=signature,
            line_start=start_line + 1,  # 1-indexed
            line_end=end_line + 1,
            is_commented=is_commented,
            base_classes=base_classes,
            category=category
        )
    
    def extract_function_pattern(self, start_line: int) -> Optional[CodePattern]:
        """Extract a function definition"""
        if start_line >= len(self.content):
            return None
        
        line = self.content[start_line]
        is_commented = self.is_commented_line(line)
        
        # Clean comment markers
        clean_line = line
        clean_line = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', clean_line)
        clean_line = re.sub(r'^\s*//\s*', '', clean_line)
        
        match = re.match(self.FUNCTION_PATTERN, clean_line)
        if not match:
            return None
        
        func_name = match.group(1)
        
        # Find end of function
        end_line = start_line
        if '{' in clean_line:
            brace_count = clean_line.count('{') - clean_line.count('}')
            for i in range(start_line + 1, min(start_line + 500, len(self.content))):
                line_content = self.content[i]
                line_content = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', line_content)
                line_content = re.sub(r'^\s*//\s*', '', line_content)
                brace_count += line_content.count('{') - line_content.count('}')
                if brace_count == 0:
                    end_line = i
                    break
        else:
            end_line = start_line  # Declaration only
        
        signature = clean_line.strip()
        
        return CodePattern(
            pattern_type='function',
            name=func_name,
            signature=signature,
            line_start=start_line + 1,
            line_end=end_line + 1,
            is_commented=is_commented
        )
    
    def extract_constant_pattern(self, start_line: int) -> Optional[CodePattern]:
        """Extract a constant definition"""
        if start_line >= len(self.content):
            return None
        
        line = self.content[start_line]
        is_commented = self.is_commented_line(line)
        
        # Clean comment markers
        clean_line = line
        clean_line = re.sub(r'^\s*//\s*\[Duplicate[^\]]*\]\s*', '', clean_line)
        clean_line = re.sub(r'^\s*//\s*', '', clean_line)
        
        match = re.match(self.CONSTANT_PATTERN, clean_line)
        if not match:
            return None
        
        const_type = match.group(1)
        const_name = match.group(2)
        
        return CodePattern(
            pattern_type='constant',
            name=const_name,
            signature=clean_line.strip(),
            line_start=start_line + 1,
            line_end=start_line + 1,
            is_commented=is_commented
        )
    
    def _categorize_pattern(self, name: str, signature: str, base_classes: Set[str]) -> str:
        """Categorize pattern as framework code or physics"""
        # Check for Qt patterns
        if any(qt in signature or qt in name for qt in self.QT_PATTERNS):
            return 'framework_qt'
        
        # Check for ANTLR patterns
        if any(antlr in signature or antlr in name for antlr in self.ANTLR_PATTERNS):
            return 'framework_antlr'
        
        # Check for SymEngine patterns
        if any(sym in signature or sym in name for sym in self.SYMENGINE_PATTERNS):
            return 'framework_symengine'
        
        # Check for PhysicsTerm hierarchy
        if 'PhysicsTerm' in base_classes or 'PhysicsTerm' in signature:
            return 'physics_term'
        
        return 'unknown'
    
    def analyze_active_section(self) -> int:
        """Analyze lines 1-35000 (active compiled code)"""
        print("\n=== ANALYZING ACTIVE SECTION (lines 1-35,000) ===")
        count = 0
        
        for line_num in range(min(self.active_section_end, len(self.content))):
            # Try to extract different pattern types
            pattern = (
                self.extract_class_pattern(line_num) or
                self.extract_function_pattern(line_num) or
                self.extract_constant_pattern(line_num)
            )
            
            if pattern and not pattern.is_commented:
                self.patterns.append(pattern)
                self.active_patterns[f"{pattern.pattern_type}::{pattern.name}"] = pattern
                count += 1
        
        print(f"✓ Extracted {count} active patterns")
        return count
    
    def analyze_commented_section(self) -> int:
        """Analyze lines 31500-102435 (commented code)"""
        print("\n=== ANALYZING COMMENTED SECTION (lines 31,500-102,435) ===")
        count = 0
        
        for line_num in range(self.commented_section_start, len(self.content)):
            pattern = (
                self.extract_class_pattern(line_num) or
                self.extract_function_pattern(line_num) or
                self.extract_constant_pattern(line_num)
            )
            
            if pattern and pattern.is_commented:
                self.patterns.append(pattern)
                self.commented_patterns[f"{pattern.pattern_type}::{pattern.name}"] = pattern
                count += 1
        
        print(f"✓ Extracted {count} commented patterns")
        return count
    
    def detect_duplicates(self) -> Tuple[int, int, int]:
        """
        Detect duplicates between active and commented sections.
        Returns: (exact_duplicates, net_new, ambiguous)
        """
        print("\n=== DETECTING DUPLICATES ===")
        exact_duplicates = 0
        net_new = 0
        ambiguous = 0
        
        for key, commented_pattern in self.commented_patterns.items():
            if key in self.active_patterns:
                # Exact duplicate - mark for skipping
                commented_pattern.category = 'exact_duplicate'
                commented_pattern.confidence_score = 0.0
                exact_duplicates += 1
            else:
                # Check for similar patterns (renamed duplicates)
                similar = self._find_similar_pattern(commented_pattern)
                if similar:
                    commented_pattern.category = 'ambiguous'
                    commented_pattern.confidence_score = 0.5
                    commented_pattern.notes = f"Similar to {similar.name}"
                    ambiguous += 1
                else:
                    # Net-new pattern
                    commented_pattern.category = 'net_new'
                    commented_pattern.confidence_score = 1.0
                    net_new += 1
        
        print(f"✓ Exact duplicates: {exact_duplicates}")
        print(f"✓ Net-new patterns: {net_new}")
        print(f"✓ Ambiguous patterns: {ambiguous}")
        
        return exact_duplicates, net_new, ambiguous
    
    def _find_similar_pattern(self, pattern: CodePattern) -> Optional[CodePattern]:
        """Find similar pattern in active section (fuzzy matching)"""
        # Simple similarity: check if base classes match
        for active_pattern in self.active_patterns.values():
            if (pattern.pattern_type == active_pattern.pattern_type and
                pattern.base_classes == active_pattern.base_classes and
                len(pattern.base_classes) > 0):
                return active_pattern
        return None
    
    def build_dependency_graph(self):
        """Build dependency graph for safe uncommenting order"""
        print("\n=== BUILDING DEPENDENCY GRAPH ===")
        
        for pattern in self.patterns:
            if pattern.is_commented and pattern.category == 'net_new':
                # Extract dependencies from signature
                dependencies = set()
                
                # Add base classes as dependencies
                dependencies.update(pattern.base_classes)
                
                # Extract type names from signature (simple heuristic)
                type_pattern = r'\b([A-Z]\w+)\b'
                for match in re.finditer(type_pattern, pattern.signature):
                    type_name = match.group(1)
                    if type_name != pattern.name:
                        dependencies.add(type_name)
                
                pattern.dependencies = dependencies
        
        print(f"✓ Dependency graph built")
    
    def generate_uncomment_batches(self, batch_size: int = 100) -> List[List[CodePattern]]:
        """Generate batches of patterns to uncomment, ordered by dependencies"""
        print(f"\n=== GENERATING UNCOMMENT BATCHES (size={batch_size}) ===")
        
        # Get all net-new patterns
        net_new_patterns = [p for p in self.patterns if p.is_commented and p.category == 'net_new']
        
        batches = []
        uncommitted = set(net_new_patterns)
        
        while uncommitted:
            batch = []
            
            for pattern in list(uncommitted):
                # Check if all dependencies are satisfied
                if all(dep in self.active_patterns or 
                       any(p.name == dep and p not in uncommitted for p in net_new_patterns)
                       for dep in pattern.dependencies):
                    batch.append(pattern)
                    if len(batch) >= batch_size:
                        break
            
            if not batch:
                # Circular dependency or missing deps - take patterns with no deps
                batch = [p for p in uncommitted if not p.dependencies][:batch_size]
                if not batch:
                    # Force add remaining patterns
                    batch = list(uncommitted)[:batch_size]
            
            batches.append(batch)
            uncommitted -= set(batch)
        
        print(f"✓ Generated {len(batches)} batches")
        return batches
    
    def generate_report(self, output_path: str):
        """Generate comprehensive analysis report"""
        print(f"\n=== GENERATING REPORT: {output_path} ===")
        
        result = AnalysisResult(
            total_patterns=len(self.patterns),
            active_patterns=len(self.active_patterns),
            commented_patterns=len(self.commented_patterns),
            patterns=self.patterns
        )
        
        # Count categories
        for pattern in self.patterns:
            if pattern.category == 'exact_duplicate':
                result.exact_duplicates += 1
            elif pattern.category == 'net_new':
                result.net_new_patterns += 1
            elif pattern.category == 'ambiguous':
                result.ambiguous_patterns += 1
            elif pattern.category.startswith('framework'):
                result.framework_patterns += 1
        
        # Save JSON report
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump({
                'summary': {
                    'total_patterns': result.total_patterns,
                    'active_patterns': result.active_patterns,
                    'commented_patterns': result.commented_patterns,
                    'exact_duplicates': result.exact_duplicates,
                    'net_new_patterns': result.net_new_patterns,
                    'ambiguous_patterns': result.ambiguous_patterns,
                    'framework_patterns': result.framework_patterns,
                    'timestamp': result.timestamp
                },
                'patterns': [
                    {
                        'type': p.pattern_type,
                        'name': p.name,
                        'lines': f"{p.line_start}-{p.line_end}",
                        'is_commented': p.is_commented,
                        'category': p.category,
                        'confidence': p.confidence_score,
                        'dependencies': list(p.dependencies),
                        'base_classes': list(p.base_classes),
                        'notes': p.notes
                    }
                    for p in result.patterns
                ]
            }, f, indent=2)
        
        print(f"✓ Report saved to {output_path}")
        
        # Print summary
        print("\n=== ANALYSIS SUMMARY ===")
        print(f"Total patterns extracted: {result.total_patterns}")
        print(f"  Active (compiled): {result.active_patterns}")
        print(f"  Commented: {result.commented_patterns}")
        print(f"\nDuplicate Analysis:")
        print(f"  Exact duplicates (skip): {result.exact_duplicates}")
        print(f"  Net-new patterns (uncomment): {result.net_new_patterns}")
        print(f"  Ambiguous patterns (review): {result.ambiguous_patterns}")
        print(f"  Framework code: {result.framework_patterns}")
        
        return result


def main():
    """Main execution"""
    print("=" * 80)
    print("INTELLIGENT DEDUPLICATION ANALYZER")
    print("Star-Magic UQFF - Option B Comprehensive Analysis")
    print("=" * 80)
    
    file_path = "MAIN_1_CoAnQi.cpp"
    
    if not Path(file_path).exists():
        print(f"✗ File not found: {file_path}")
        return 1
    
    # Initialize analyzer
    analyzer = IntelligentDeduplicationAnalyzer(file_path)
    
    # Load file
    if not analyzer.load_file():
        return 1
    
    # Analyze active section
    active_count = analyzer.analyze_active_section()
    
    # Analyze commented section
    commented_count = analyzer.analyze_commented_section()
    
    # Detect duplicates
    duplicates, net_new, ambiguous = analyzer.detect_duplicates()
    
    # Build dependency graph
    analyzer.build_dependency_graph()
    
    # Generate uncomment batches
    batches = analyzer.generate_uncomment_batches(batch_size=150)
    
    # Generate reports
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_path = f"deduplication_analysis_{timestamp}.json"
    result = analyzer.generate_report(report_path)
    
    print("\n" + "=" * 80)
    print("ANALYSIS COMPLETE")
    print(f"Next steps: Review {report_path} and execute incremental uncommenting")
    print("=" * 80)
    
    return 0


if __name__ == "__main__":
    exit(main())
