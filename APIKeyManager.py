#!/usr/bin/env python3
"""
API Key Configuration Manager for SuperGrok4
Handles persistent storage and retrieval of API credentials
"""

import json
import os
from pathlib import Path
from typing import Dict, Optional

CONFIG_FILE = Path(__file__).parent / "grok_api_config.json"

# Default config structure
DEFAULT_CONFIG = {
    "api_keys": {
        "xai_grok": "",
        "openai": ""
    },
    "preferences": {
        "default_provider": "grok",
        "auto_save": True,
        "remember_model": True,
        "last_model": "grok-4-1-fast-reasoning"
    },
    "version": "1.0"
}

def load_api_config() -> Dict:
    """Load API configuration from file."""
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load config file: {e}")
            return DEFAULT_CONFIG.copy()
    return DEFAULT_CONFIG.copy()

def save_api_config(config: Dict) -> bool:
    """Save API configuration to file."""
    try:
        with open(CONFIG_FILE, 'w', encoding='utf-8') as f:
            json.dump(config, f, indent=4, ensure_ascii=False)
        return True
    except IOError as e:
        print(f"Error: Could not save config file: {e}")
        return False

def get_xai_api_key() -> str:
    """
    Get xAI API key from multiple sources (in priority order):
    1. Config file (saved by user in Tab 7)
    2. Environment variable (XAI_API_KEY)
    3. Return empty string if not found
    """
    # Try config file first
    config = load_api_config()
    key = config.get('api_keys', {}).get('xai_grok', '').strip()
    if key:
        return key
    
    # Fall back to environment variable
    key = os.environ.get('XAI_API_KEY', '').strip()
    if key:
        return key
    
    return ""

def set_xai_api_key(key: str) -> bool:
    """Save xAI API key to config file."""
    config = load_api_config()
    config['api_keys']['xai_grok'] = key.strip()
    config['preferences']['auto_save'] = True
    return save_api_config(config)

def get_config_file_path() -> str:
    """Return path to config file for display."""
    return str(CONFIG_FILE)

def config_exists() -> bool:
    """Check if config file exists."""
    return CONFIG_FILE.exists()

def has_saved_api_key() -> bool:
    """Check if API key is saved in config file (not just env var)."""
    config = load_api_config()
    return bool(config.get('api_keys', {}).get('xai_grok', '').strip())

def get_api_key_status() -> str:
    """Return human-readable status of API key availability."""
    if has_saved_api_key():
        return "✅ Saved in config file (grok_api_config.json)"
    elif os.environ.get('XAI_API_KEY'):
        return "✅ Set via environment variable (XAI_API_KEY)"
    else:
        return "❌ Not configured - Click 'Configure API Key' to set up"

# Test function
if __name__ == "__main__":
    print("Testing API Key Manager...")
    print(f"Config file location: {CONFIG_FILE}")
    print(f"Config file exists: {config_exists()}")
    print(f"Saved API key present: {has_saved_api_key()}")
    print(f"Status: {get_api_key_status()}")
    print(f"Current XAI API key: {get_xai_api_key()[:10]}..." if get_xai_api_key() else "No API key found")
