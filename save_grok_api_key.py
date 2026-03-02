#!/usr/bin/env python3
"""
Save your actual Grok API key to the config file.
This enables persistent storage without environment variables.
"""
from APIKeyManager import set_xai_api_key, get_api_key_status
import os

# Get API key from environment variable or prompt user
api_key = os.environ.get('XAI_API_KEY', '')

if not api_key:
    print("=" * 70)
    print("GROK API KEY SETUP")
    print("=" * 70)
    print("\nTo get your Grok API key:")
    print("  1. Visit: https://console.xai.com/api-keys")
    print("  2. Create a new API key")
    print("  3. Copy it below (it starts with 'xai-')")
    print()
    api_key = input("Enter your Grok API key: ").strip()

if not api_key:
    print("❌ No API key provided. Exiting.")
    exit(1)

if not api_key.startswith('xai-'):
    print("⚠️  Warning: API key doesn't start with 'xai-'")
    print("   This might not be a valid Grok API key")
    confirm = input("Continue anyway? (y/n): ").strip().lower()
    if confirm != 'y':
        exit(1)

print("\nSaving API key to config file...")
try:
    set_xai_api_key(api_key)
    print("✅ API key saved successfully!\n")
    print(get_api_key_status())
    print("\n" + "=" * 70)
    print("✅ Setup complete! API key will be used automatically.")
    print("   No app restart needed - Python layer uses it immediately.")
    print("=" * 70)
except Exception as e:
    print(f"❌ Failed to save API key: {e}")
    exit(1)
