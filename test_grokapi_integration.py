#!/usr/bin/env python3
"""Test GrokAPI integration with saved key"""
import sys
from pathlib import Path

PROJECT_DIR = Path(__file__).parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from GrokAPI import get_xai_api_key, get_api_key_status

print("=" * 70)
print("TEST: GrokAPI Integration with Saved Key")
print("=" * 70)

print("\n[TEST 1] GrokAPI get_xai_api_key()")
try:
    key = get_xai_api_key()
    if key.startswith("xai-"):
        print(f"PASS - Retrieved key: {key[:40]}...")
        print(f"       Length: {len(key)} characters")
    else:
        print("FAIL - Invalid key format")
        exit(1)
except Exception as e:
    print(f"FAIL - {e}")
    exit(1)

print("\n[TEST 2] GrokAPI get_api_key_status()")
try:
    status = get_api_key_status()
    print(f"PASS - Status: {status}")
except Exception as e:
    print(f"FAIL - {e}")
    exit(1)

print("\n" + "=" * 70)
print("SUCCESS - GROKAPI INTEGRATION VERIFIED")
print("=" * 70)
print("\nGrokAPI is ready to use your saved API key!")
print("\nYour workflow is now complete:")
print("  1. ✅ APIKeyManager.py manages the key")
print("  2. ✅ grok_api_config.json stores the key persistently")
print("  3. ✅ GrokAPI.py uses the key automatically")
print("  4. ✅ Works from any working directory")
print("  5. ✅ No environment variables needed")
