#!/usr/bin/env python3
"""
Phase 1: Extract closure metadata from papers 001-1211 for ledger integration.

This script scans all PAPER_*.md files in whitepapers/ directory and attempts to:
1. Extract measurement/prediction data from paper content
2. Infer closure labels from paper titles and abstracts
3. Build ledger-ready CSV rows with {closure, predicted, observed, error_pct, ...}

Output: _phase1_paper_closures.json and _phase1_extraction_report.txt
"""

import json
import re
from pathlib import Path
from collections import defaultdict
from typing import Optional, Dict, List, Tuple
import yaml

# Configuration
WHITEPAPERS_DIR = Path('whitepapers')
OUTPUT_JSON = Path('_phase1_paper_closures.json')
OUTPUT_REPORT = Path('_phase1_extraction_report.txt')

# Regex patterns for measurement extraction
PATTERNS = {
    'percent': r'([\d\.eE\-\+]+)\s*%',
    'ratio': r'ratio[:\s]+([\d\.eE\-\+]+)',
    'error': r'error[:\s]+([\d\.eE\-\+%]+)',
    'reduction': r'([\d\.]+)\s*%\s*reduction',
    'tension': r'tension[:\s]+([\d\.eE\-\+]+)',
    'observed': r'observed[:\s]+([\'"]?[\d\.eE\-\+]+[\'"]?)',
    'predicted': r'predicted[:\s]+([\'"]?[\d\.eE\-\+]+[\'"]?)',
    'measured': r'measured[:\s]+([\'"]?[\d\.eE\-\+]+[\'"]?)',
    'value': r'(?:value|result)[:\s]+([\'"]?[\d\.eE\-\+]+[\'"]?)',
}

def extract_yaml(content: str) -> dict:
    """Extract YAML frontmatter from paper."""
    match = re.match(r'^---\n(.*?)\n---', content, re.DOTALL)
    if match:
        try:
            return yaml.safe_load(match.group(1)) or {}
        except:
            return {}
    return {}

def extract_closure_from_content(paper_id: str, title: str, content: str) -> Optional[Dict]:
    """
    Attempt to extract closure data from paper content.
    
    Returns dict with: {label, predicted, observed, error_pct, inferred_from}
    """
    
    # Extract first 2000 chars of abstract for pattern matching
    abstract_match = re.search(r'## Abstract\n(.*?)(?:\n##|\n---|\Z)', content, re.DOTALL)
    abstract = abstract_match.group(1)[:2000] if abstract_match else content[:2000]
    
    closure_data = {
        'paper_id': paper_id,
        'label': None,
        'predicted': None,
        'observed': None,
        'error_pct': None,
        'inferred_from': []
    }
    
    # Try to build label from title + paper_id
    # Pattern: PAPER_NNN_DescriptiveTitle
    label = paper_id + '_' + re.sub(r'[^A-Za-z0-9_]', '_', title.lower())[:50]
    label = re.sub(r'_+', '_', label).strip('_')
    closure_data['label'] = label
    
    # Try to extract numeric values from abstract
    predicted_vals = []
    observed_vals = []
    errors = []
    
    # Pattern: "X% reduction" or similar
    for pattern_name, pattern in PATTERNS.items():
        matches = re.findall(pattern, abstract, re.IGNORECASE)
        for match in matches[:3]:  # Limit to first 3 matches per pattern
            try:
                val = float(match.strip('%').strip('"\''))
                if pattern_name in ['reduction', 'percent', 'error', 'tension']:
                    errors.append(val)
                elif pattern_name in ['predicted']:
                    predicted_vals.append(val)
                elif pattern_name in ['observed', 'measured']:
                    observed_vals.append(val)
                elif pattern_name in ['ratio', 'value']:
                    # Could be either; prefer observed
                    observed_vals.append(val)
            except ValueError:
                pass
    
    # Assign extracted values
    if predicted_vals:
        closure_data['predicted'] = predicted_vals[0]
        closure_data['inferred_from'].append('predicted_pattern')
    
    if observed_vals:
        closure_data['observed'] = observed_vals[0]
        closure_data['inferred_from'].append('observed_pattern')
    
    if errors:
        closure_data['error_pct'] = errors[0]
        closure_data['inferred_from'].append('error_pattern')
    
    # Try to infer from title
    if 'UQFF' in title and ('Analysis' in title or 'Calculation' in title):
        closure_data['inferred_from'].append('title_analysis')
    if 'Damping' in title:
        closure_data['inferred_from'].append('damping_detected')
    if 'Merger' in title:
        closure_data['inferred_from'].append('merger_detected')
    
    return closure_data

def build_ledger_row(paper_id: str, closure_data: Dict, next_id: int) -> Dict:
    """Build a ledger CSV row from extracted closure data."""
    
    status = 'OK'
    if not closure_data['predicted'] and not closure_data['observed']:
        status = 'RESEARCH_TRACE'  # Paper exists, closure incomplete
    elif closure_data['predicted'] and closure_data['observed']:
        status = 'OK'
    else:
        status = 'RESEARCH_TRACE'
    
    # Calculate error % if both predicted and observed exist
    error_pct = closure_data['error_pct']
    if not error_pct and closure_data['predicted'] and closure_data['observed']:
        if closure_data['predicted'] != 0:
            error_pct = abs((closure_data['observed'] - closure_data['predicted']) / closure_data['predicted']) * 100
    
    row = {
        'ID': next_id,
        'closure': closure_data['label'],
        'predicted': closure_data['predicted'] or '',
        'observed': closure_data['observed'] or '',
        'error_pct': round(error_pct, 3) if error_pct else '',
        'status': status,
        'cvw_stamp': 'v2.0.0',
        'sm_anchor': 'CVW v2.0.0 — G6 SM Anchor Gate compliant',
        'label': paper_id,
        'raw_output': json.dumps(closure_data['inferred_from']),
        'category': 'paper_index',
        'name': paper_id,
        'script': f'whitepapers/{paper_id}.md'
    }
    
    return row

def main():
    """Execute Phase 1 closure extraction."""
    
    # Scan all papers
    papers = sorted(WHITEPAPERS_DIR.glob('PAPER_*.md'))
    print(f'[Phase 1] Found {len(papers)} papers in whitepapers/')
    
    extracted = []
    ledger_rows = []
    next_id = 641  # Start after existing 640 rows
    
    report_lines = [
        f'Phase 1 Closure Extraction Report',
        f'=' * 60,
        f'Total papers scanned: {len(papers)}',
        f'',
        f'Extraction Results:',
        f'-' * 60,
    ]
    
    stats = defaultdict(int)
    
    for i, paper_path in enumerate(papers):
        if i % 200 == 0:
            print(f'  Processing {i}/{len(papers)}...')
        
        try:
            content = paper_path.read_text(encoding='utf-8', errors='ignore')
            yaml_data = extract_yaml(content)
            
            paper_id = yaml_data.get('paper_id', paper_path.stem)
            title = yaml_data.get('title', paper_path.stem)
            
            # Extract closure
            closure_data = extract_closure_from_content(paper_id, title, content)
            extracted.append(closure_data)
            
            # Build ledger row
            row = build_ledger_row(paper_id, closure_data, next_id)
            ledger_rows.append(row)
            
            # Track statistics
            stats[row['status']] += 1
            if closure_data['predicted']:
                stats['has_predicted'] += 1
            if closure_data['observed']:
                stats['has_observed'] += 1
            
            next_id += 1
            
        except Exception as e:
            report_lines.append(f'ERROR processing {paper_path.stem}: {str(e)[:60]}')
            stats['errors'] += 1
    
    # Write outputs
    print(f'[Phase 1] Writing outputs...')
    
    # JSON output (full extraction data)
    with open(OUTPUT_JSON, 'w', encoding='utf-8') as f:
        json.dump({
            'total_papers': len(papers),
            'extracted_closures': extracted,
            'ledger_rows': ledger_rows,
            'statistics': dict(stats)
        }, f, indent=2)
    
    # Report
    report_lines.extend([
        f'',
        f'Statistics:',
        f'  Total: {len(papers)}',
        f'  Status=OK: {stats["OK"]}',
        f'  Status=RESEARCH_TRACE: {stats["RESEARCH_TRACE"]}',
        f'  With predicted value: {stats["has_predicted"]}',
        f'  With observed value: {stats["has_observed"]}',
        f'  Errors: {stats["errors"]}',
        f'',
        f'Next IDs assigned: {641}-{next_id-1}',
        f'Total new ledger rows: {len(ledger_rows)}',
        f'',
        f'Next Phase: Merge extracted_closures into master_closures.csv',
    ])
    
    with open(OUTPUT_REPORT, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f'[Phase 1] Complete!')
    print(f'  Outputs:')
    print(f'    {OUTPUT_JSON}')
    print(f'    {OUTPUT_REPORT}')
    print(f'')
    print(f'Statistics:')
    for key, val in sorted(stats.items()):
        print(f'  {key}: {val}')

if __name__ == '__main__':
    main()
