#!/usr/bin/env python3
"""Test saving the user's provided API key"""
from APIKeyManager import set_xai_api_key, get_xai_api_key, get_api_key_status
import json
from pathlib import Path

test_key = "xai-SA7qHpjU2kC6H5B7QEiB9OfEGzuFj4sB4EK4P2X3pj7HGSelVhQDgFOHKOJ99kGVUCd4o3pj3MU1q6AP"

print("=" * 70)
print("TEST: SAVING PROVIDED API KEY")
print("=" * 70)

print(f"\nTest key: {test_key[:40]}...")
print(f"Key length: {len(test_key)} characters")

# Step 1: Save
print("\n[STEP 1] Saving API key...")
try:
    set_xai_api_key(test_key)
    print("PASS - Save successful")
except Exception as e:
    print(f"FAIL - Save failed: {e}")
    import traceback
    traceback.print_exc()
    exit(1)

# Step 2: Retrieve
print("\n[STEP 2] Retrieving API key...")
try:
    retrieved = get_xai_api_key()
    if retrieved == test_key:
        print("PASS - Key matches exactly!")
    else:
        print(f"FAIL - Key mismatch")
        print(f"  Expected length: {len(test_key)}")
        print(f"  Retrieved length: {len(retrieved)}")
        exit(1)
except Exception as e:
    print(f"FAIL - Retrieval failed: {e}")
    exit(1)

# Step 3: Check status
print("\n[STEP 3] Checking API key status...")
try:
    status = get_api_key_status()
    print(f"Status: {status}")
    if "Saved" in status:
        print("PASS - Status indicates saved key")
    else:
        print("WARNING - Status unclear")
except Exception as e:
    print(f"FAIL - Status check failed: {e}")

# Step 4: Verify config file
print("\n[STEP 4] Verifying config file...")
config_file = Path(__file__).parent / "grok_api_config.json"
if config_file.exists():
    with open(config_file, 'r') as f:
        config = json.load(f)
    saved_key = config.get("api_keys", {}).get("xai_grok", "")
    
    if saved_key == test_key:
        print(f"PASS - Config file contains correct key")
        print(f"  File: {config_file}")
        print(f"  Size: {config_file.stat().st_size} bytes")
    else:
        print(f"FAIL - Config file key mismatch")
        print(f"  Expected: {test_key[:30]}...")
        print(f"  Got: {saved_key[:30]}...")
        exit(1)
else:
    print(f"FAIL - Config file not found at {config_file}")
    exit(1)

print("\n" + "=" * 70)
print("SUCCESS - API KEY SAVE AND VERIFY PASSED")
print("=" * 70)
print(f"\nYour API key has been successfully saved to grok_api_config.json")
print(f"GrokAPI will use this key automatically without environment variables.")
