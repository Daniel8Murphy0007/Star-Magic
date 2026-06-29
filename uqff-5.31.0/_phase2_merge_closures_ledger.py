#!/usr/bin/env python3
"""
Phase 2: Merge extracted paper closures into master_closures.csv ledger.

This script:
1. Reads extracted closures from _phase1_paper_closures.json
2. Reads current master_closures.csv
3. Appends new rows to ledger
4. Validates schema consistency
5. Outputs updated master_closures.csv + backup
"""

import json
import csv
from pathlib import Path
from datetime import datetime

INPUT_EXTRACTED = Path('_phase1_paper_closures.json')
LEDGER_PATH = Path('master_closures.csv')
LEDGER_BACKUP = Path('master_closures.csv.backup_phase2')
OUTPUT_LEDGER = Path('master_closures.csv')

# Expected ledger schema
EXPECTED_COLS = [
    'ID', 'closure', 'predicted', 'observed', 'error_pct', 'status',
    'cvw_stamp', 'sm_anchor', 'label', 'raw_output', 'category', 'name', 'script'
]

def main():
    """Execute Phase 2 ledger merge."""
    
    print('[Phase 2] Starting ledger merge...')
    
    # Load extracted closures
    print('[Phase 2] Loading extracted closures...')
    with open(INPUT_EXTRACTED, 'r', encoding='utf-8') as f:
        extracted_data = json.load(f)
    
    ledger_rows = extracted_data['ledger_rows']
    print(f'  Loaded {len(ledger_rows)} extracted closure rows')
    
    # Load current ledger
    print('[Phase 2] Loading current ledger...')
    with open(LEDGER_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        current_rows = list(reader)
        current_cols = reader.fieldnames or []
    
    print(f'  Current ledger: {len(current_rows)} rows, {len(current_cols)} columns')
    
    # Backup current ledger
    LEDGER_PATH.write_text(LEDGER_PATH.read_text(encoding='utf-8'), encoding='utf-8')
    LEDGER_BACKUP.write_text(LEDGER_PATH.read_text(encoding='utf-8'), encoding='utf-8')
    print(f'  Backup saved: {LEDGER_BACKUP}')
    
    # Validate schema consistency
    print('[Phase 2] Validating schema...')
    current_col_set = set(current_cols)
    expected_col_set = set(EXPECTED_COLS)
    
    missing_cols = expected_col_set - current_col_set
    extra_cols = current_col_set - expected_col_set
    
    if missing_cols:
        print(f'  WARNING: Missing columns: {missing_cols}')
    if extra_cols:
        print(f'  NOTE: Extra columns in ledger: {extra_cols}')
    
    # Merge rows
    print('[Phase 2] Merging rows...')
    
    # Track stats
    added_count = 0
    duplicate_count = 0
    
    existing_closures = {r.get('closure', ''): r for r in current_rows}
    
    for new_row in ledger_rows:
        closure_name = new_row.get('closure', '')
        
        # Check for duplicates
        if closure_name in existing_closures:
            print(f'  DUPLICATE: {closure_name} (skipping)')
            duplicate_count += 1
            continue
        
        # Ensure all expected columns present
        for col in EXPECTED_COLS:
            if col not in new_row:
                new_row[col] = ''
        
        current_rows.append(new_row)
        added_count += 1
    
    print(f'  Added {added_count} new rows')
    print(f'  Skipped {duplicate_count} duplicates')
    print(f'  Final ledger size: {len(current_rows)} rows')
    
    # Write updated ledger
    print('[Phase 2] Writing updated ledger...')
    
    # Use original column order, adding any new columns
    output_cols = list(current_cols)
    for col in EXPECTED_COLS:
        if col not in output_cols:
            output_cols.append(col)
    
    with open(OUTPUT_LEDGER, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=output_cols)
        writer.writeheader()
        
        # Write rows in order of ID if ID column exists
        rows_to_write = current_rows
        if 'ID' in output_cols:
            try:
                rows_to_write = sorted(
                    current_rows,
                    key=lambda r: int(r.get('ID', '0') or '0')
                )
            except:
                pass  # Keep original order if sorting fails
        
        writer.writerows(rows_to_write)
    
    print(f'  Updated ledger: {OUTPUT_LEDGER}')
    print(f'  Total rows: {len(rows_to_write)}')
    
    # Validate output
    print('[Phase 2] Validating output...')
    with open(OUTPUT_LEDGER, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        verified_rows = list(reader)
    
    print(f'  Verified: {len(verified_rows)} rows read back')
    
    # Statistics
    status_counts = {}
    for row in verified_rows:
        status = row.get('status', 'UNKNOWN')
        status_counts[status] = status_counts.get(status, 0) + 1
    
    print(f'')
    print(f'Final Ledger Statistics:')
    for status, count in sorted(status_counts.items()):
        print(f'  {status}: {count}')
    
    print(f'')
    print(f'[Phase 2] Ledger merge complete!')

if __name__ == '__main__':
    main()
