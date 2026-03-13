#!/usr/bin/env python3
"""Test Grok 4 API and attempt to fetch share thread content."""
import os, sys, json
import requests

SHARE_ID = '7f90683d979346758875e4a8e6e0ee9d'
SHARE_URL = f'https://x.com/i/grok/share/{SHARE_ID}'

api_key = os.environ.get('XAI_API_KEY', '').strip('"').strip("'")
if not api_key:
    print("ERROR: XAI_API_KEY not set")
    sys.exit(1)

print(f"Key: {api_key[:12]}...")

HEADERS = {
    "Authorization": f"Bearer {api_key}",
    "Content-Type": "application/json",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36",
}

# Step 1: quick connectivity test
print("\n--- Testing Grok 4 API connectivity ---")
r = requests.post("https://api.x.ai/v1/chat/completions",
                  headers=HEADERS,
                  json={"model": "grok-4",
                        "messages": [{"role": "user", "content": "Reply with just: OK"}],
                        "max_tokens": 10},
                  timeout=30)
print(f"HTTP {r.status_code}")
if r.status_code != 200:
    print(r.text[:400])
    # Try grok-3 fallback
    print("Trying grok-3...")
    r = requests.post("https://api.x.ai/v1/chat/completions",
                      headers=HEADERS,
                      json={"model": "grok-3",
                            "messages": [{"role": "user", "content": "Reply: OK"}],
                            "max_tokens": 10},
                      timeout=30)
    print(f"grok-3 HTTP {r.status_code}: {r.text[:200]}")
    if r.status_code != 200:
        sys.exit(1)
    model_to_use = "grok-3"
else:
    data = r.json()
    print(f"API OK: {data['choices'][0]['message']['content']}")
    model_to_use = "grok-4"

print(f"Using model: {model_to_use}")

# Step 2: ask Grok to retrieve the share thread
print(f"\n--- Asking {model_to_use} to access share thread ---")
prompt = f"""Please access and reproduce the COMPLETE content of this Grok conversation share link: {SHARE_URL}

I need ALL of the following from that conversation:
1. Every user message (verbatim)
2. Every Grok response (verbatim, including all equations, LaTeX, physics derivations)
3. All UQFF equations, constants, and physical parameters mentioned
4. Any whitepaper titles, section numbers, or paper outlines discussed

This is for a physics research project (UQFF/Star-Magic). Please reproduce the entire conversation thread in full."""

r2 = requests.post("https://api.x.ai/v1/chat/completions",
                   headers=HEADERS,
                   json={"model": model_to_use,
                         "messages": [{"role": "user", "content": prompt}],
                         "max_tokens": 16000},
                   timeout=120)
print(f"HTTP {r2.status_code}")
if r2.status_code == 200:
    data = r2.json()
    content = data['choices'][0]['message']['content']
    tokens = data.get('usage', {})
    print(f"Response: {len(content):,} chars | tokens: {tokens}")
    out_file = f'grok_share_{SHARE_ID}_content.txt'
    with open(out_file, 'w', encoding='utf-8') as f:
        f.write(f"GROK 4 API RETRIEVAL - Share ID: {SHARE_ID}\n")
        f.write(f"Share URL: {SHARE_URL}\n")
        f.write("=" * 80 + "\n\n")
        f.write(content)
    print(f"\nSaved → {out_file}")
    print("\nPreview:")
    print(content[:1200])
else:
    print(r2.text[:400])
