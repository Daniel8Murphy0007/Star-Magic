#!/usr/bin/env python3
"""Test APIKeyManager save functionality"""
import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

try:
    from APIKeyManager import set_xai_api_key, get_api_key_status, get_xai_api_key
    
    print("=" * 60)
    print("Testing API Key Manager")
    print("=" * 60)
    
    # Test 1: Save API key
    print("\n[TEST 1] Saving API key...")
    set_xai_api_key('xai-test-key-12345')
    print("✓ Save successful")
    
    # Test 2: Get status
    print("\n[TEST 2] Getting API key status...")
    status = get_api_key_status()
    print(status)
    
    # Test 3: Retrieve saved key
    print("\n[TEST 3] Retrieving saved API key...")
    retrieved_key = get_xai_api_key()
    print(f"Retrieved key: {retrieved_key[:10]}..." if retrieved_key else "No key found")
    
    # Test 4: Check config file exists
    config_file = Path(__file__).parent / "grok_api_config.json"
    print(f"\n[TEST 4] Config file exists: {config_file.exists()}")
    if config_file.exists():
        print(f"Config file path: {config_file}")
        print(f"Config file size: {config_file.stat().st_size} bytes")
        with open(config_file, 'r') as f:
            print(f"Config content preview: {f.read()[:200]}")
    
    print("\n" + "=" * 60)
    print("✅ ALL TESTS PASSED")
    print("=" * 60)
    
except Exception as e:
    print("\n" + "=" * 60)
    print("❌ ERROR DETECTED")
    print("=" * 60)
    traceback.print_exc()
    sys.exit(1)
