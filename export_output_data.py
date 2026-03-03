#!/usr/bin/env python3
"""
export_output_data.py - Export CondensedPhysics_OutputData.py to JSON for source2.cpp

Called by SessionLogWidget (Tab 9) to retrieve query history for display.
Exports OUTPUT_STORE to query_history.json for C++ GUI consumption.

Usage:
    python export_output_data.py

Output:
    query_history.json - Complete query history in JSON format
"""

import sys
import os

try:
    from CondensedPhysics_OutputData import OUTPUT_STORE
    
    # Export to JSON
    output_file = "query_history.json"
    OUTPUT_STORE.export_to_json(output_file)
    
    print(f"✓ Exported {len(OUTPUT_STORE._query_history)} queries to {output_file}", file=sys.stderr)
    sys.exit(0)
    
except ImportError as e:
    print(f"⚠️ Error: Could not import CondensedPhysics_OutputData: {e}", file=sys.stderr)
    # Create empty output so SessionLogWidget doesn't crash
    import json
    with open("query_history.json", 'w') as f:
        json.dump({"query_history": [], "results": {}}, f)
    sys.exit(1)
    
except Exception as e:
    print(f"⚠️ Error exporting query history: {e}", file=sys.stderr)
    # Create empty output
    import json
    with open("query_history.json", 'w') as f:
        json.dump({"query_history": [], "results": {}}, f)
    sys.exit(1)
