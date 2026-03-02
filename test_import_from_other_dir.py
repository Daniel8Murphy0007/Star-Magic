#!/usr/bin/env python3
"""Test APIKeyManager import from different working directory"""
import sys
import os

# Change to different directory
os.chdir('c:\\Users\\tmsjd')
print(f"Working directory: {os.getcwd()}")

# Add project to path
sys.path.insert(0, r'c:\Users\tmsjd\source\repos\Daniel8Murphy0007\Star-Magic')

try:
    print("\nImporting GrokAPI...")
    from GrokAPI import get_xai_api_key, get_api_key_status
    print("✅ GrokAPI imported successfully")
    
    print("\nImporting APIKeyManager directly...")
    from APIKeyManager import get_xai_api_key as get_key, get_api_key_status as get_status
    print("✅ APIKeyManager imported successfully")
    
    print("\nAPI Status:")
    print(get_api_key_status())
    
    print("\n✅ ALL IMPORTS SUCCESSFUL")
except Exception as e:
    print(f"❌ Import failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)
