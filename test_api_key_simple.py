#!/usr/bin/env python3
"""
Comprehensive API Key Manager Tests
Tests all functionality of the Python API key management system
"""
import sys
import os
import json
import tempfile
from pathlib import Path

# Add project to path
sys.path.insert(0, str(Path(__file__).parent))

from APIKeyManager import (
    load_api_config,
    save_api_config, 
    get_xai_api_key,
    set_xai_api_key,
    get_api_key_status
)

def run_tests():
    """Run all API key manager tests"""
    print("=" * 70)
    print("API KEY MANAGER - COMPREHENSIVE TESTS")
    print("=" * 70)
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Get API Key Status
    print("\n[TEST 1] Get API Key Status")
    try:
        status = get_api_key_status()
        print(f"Status: {status}")
        if "Saved" in status or "config" in status:
            print("PASS")
            tests_passed += 1
        else:
            print("FAIL - unexpected status")
            tests_failed += 1
    except Exception as e:
        print(f"FAIL - {e}")
        tests_failed += 1
    
    # Test 2: Set API Key
    print("\n[TEST 2] Set API Key to config file")
    try:
        set_xai_api_key("xai-test-key-final-12345")
        print("API key set successfully")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL - {e}")
        tests_failed += 1
    
    # Test 3: Get API Key
    print("\n[TEST 3] Retrieve API Key from config")
    try:
        key = get_xai_api_key()
        if key and "xai-test-key-final-12345" in key:
            print(f"Retrieved key: {key[:15]}...")
            tests_passed += 1
        else:
            print(f"FAIL - key not found or mismatch")
            tests_failed += 1
    except Exception as e:
        print(f"FAIL - {e}")
        tests_failed += 1
    
    # Test 4: Config file exists
    print("\n[TEST 4] Verify config file exists")
    config_file = Path(__file__).parent / "grok_api_config.json"
    if config_file.exists():
        print(f"Config file: {config_file}")
        size = config_file.stat().st_size
        print(f"Config file size: {size} bytes")
        tests_passed += 1
    else:
        print(f"FAIL - config file not found")
        tests_failed += 1
    
    # Test 5: Check GrokAPI integration
    print("\n[TEST 5] GrokAPI Integration")
    try:
        from GrokAPI import get_xai_api_key as grok_get_key, get_api_key_status as grok_status
        status = grok_status()
        print(f"GrokAPI status: {status}")
        tests_passed += 1
    except Exception as e:
        print(f"FAIL - {e}")
        tests_failed += 1
    
    # Summary
    print("\n" + "=" * 70)
    print(f"RESULTS: {tests_passed} PASSED, {tests_failed} FAILED")
    if tests_failed == 0:
        print("SUCCESS - All tests passed!")
    else:
        print("FAILURE - Some tests failed")
    print("=" * 70)
    
    return tests_failed == 0

if __name__ == "__main__":
    success = run_tests()
    sys.exit(0 if success else 1)
