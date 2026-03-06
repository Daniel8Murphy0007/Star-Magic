#!/usr/bin/env python3
"""
Fetch Grok share thread 0904a12a using xAI API with live search.
The shared conversation URL is: https://x.com/i/grok/share/0904a12a5c2b4a639389ae084391b94f
"""
import sys
import json
import requests
from pathlib import Path

PROJECT_DIR = Path(__file__).parent
if str(PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(PROJECT_DIR))

from APIKeyManager import get_xai_api_key

API_KEY = get_xai_api_key()
URL = "https://api.x.ai/v1/chat/completions"
SHARE_URL = "https://x.com/i/grok/share/0904a12a5c2b4a639389ae084391b94f"

headers = {
    "Content-Type": "application/json",
    "Authorization": f"Bearer {API_KEY}"
}

# Try with live_search tool to browse to the shared URL
payload_search = {
    "model": "grok-3",
    "messages": [
        {
            "role": "user",
            "content": (
                f"Please access this Grok shared conversation URL and return its COMPLETE, VERBATIM content "
                f"including all user messages and all Grok responses, every equation, every code block, "
                f"and every technical detail: {SHARE_URL}\n\n"
                f"Return the full conversation transcript with all physics equations, UQFF terms, "
                f"mathematical derivations, and any code that was generated."
            )
        }
    ],
    "tools": [
        {
            "type": "function",
            "function": {
                "name": "web_search",
                "description": "Search the web",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "query": {"type": "string"}
                    },
                    "required": ["query"]
                }
            }
        }
    ],
    "temperature": 0.1
}

# First attempt: grok-3 with live_search
print("Attempt 1: grok-3 with search tool...")
resp = requests.post(URL, headers=headers, json=payload_search, timeout=120)
print(f"  Status: {resp.status_code}")

if resp.status_code == 200:
    data = resp.json()
    content = data.get("choices", [{}])[0].get("message", {}).get("content", "")
    print(f"  Content length: {len(content)}")
    print(f"\n{'='*80}\nRESPONSE:\n{'='*80}")
    print(content)
else:
    print(f"  Error: {resp.text[:500]}")
    
    # Second attempt: try grok-3-latest without tool, direct URL in prompt
    print("\nAttempt 2: grok-3-latest without tool...")
    payload2 = {
        "model": "grok-3-latest",
        "messages": [
            {
                "role": "user",
                "content": (
                    f"I shared a Grok conversation at: {SHARE_URL}\n\n"
                    f"This conversation is from MY OWN Grok account. Can you access and retrieve "
                    f"the content of this shared conversation? Please return the complete transcript."
                )
            }
        ],
        "temperature": 0.1
    }
    resp2 = requests.post(URL, headers=headers, json=payload2, timeout=120)
    print(f"  Status: {resp2.status_code}")
    if resp2.status_code == 200:
        data2 = resp2.json()
        content2 = data2.get("choices", [{}])[0].get("message", {}).get("content", "")
        print(f"  Content length: {len(content2)}")
        print(content2)
    else:
        print(f"  Error: {resp2.text[:500]}")

# Also check available models
print("\n\nChecking available models...")
resp_models = requests.get("https://api.x.ai/v1/models", headers=headers, timeout=30)
if resp_models.status_code == 200:
    models = resp_models.json()
    model_ids = [m["id"] for m in models.get("data", [])]
    print("Available models:", model_ids)
else:
    print(f"Models error: {resp_models.status_code}")
