#!/usr/bin/env python3
"""Fetch Grok share thread 381a8fe78e1a4ecbaf32a88aa386df30 via API."""
import os, sys, json
import requests

SHARE_ID = "381a8fe78e1a4ecbaf32a88aa386df30"
OUT_FILE = "grok_share_381a8fe7.txt"

k = os.environ.get("XAI_API_KEY", "").strip().strip('"').strip("'")
if not k:
    print("ERROR: XAI_API_KEY not set")
    sys.exit(1)
print(f"Key: {k[:12]}... ({len(k)} chars)")

HEADERS = {
    "Authorization": f"Bearer {k}",
    "Content-Type": "application/json",
}

# Step 1: connectivity
print("Testing API...")
try:
    r = requests.post("https://api.x.ai/v1/chat/completions",
        headers=HEADERS,
        json={"model":"grok-4","messages":[{"role":"user","content":"Reply: READY"}],"max_tokens":10},
        timeout=30)
    print(f"HTTP {r.status_code}: {r.text[:100]}")
except Exception as e:
    print(f"Connection error: {e}")
    sys.exit(1)

if r.status_code != 200:
    print("API not reachable")
    sys.exit(1)

# Step 2: ask for thread content
print(f"\nRequesting content of share {SHARE_ID}...")
prompt = f"""I am a developer working on UQFF Star-Magic physics project.
I need you to reproduce the COMPLETE content of the Grok conversation at:
https://x.com/i/grok/share/{SHARE_ID}

This session covers UQFF astrophysics equations. Please output:
1. All physics equations and mathematics discussed
2. All Python class definitions, method signatures, and code
3. All constants, coefficients, and physical parameters
4. All C++ function signatures if any
5. Complete derivations in full detail

Output everything — do not truncate or summarize."""

r2 = requests.post("https://api.x.ai/v1/chat/completions",
    headers=HEADERS,
    json={"model":"grok-4","messages":[{"role":"user","content":prompt}],"max_tokens":16000,"temperature":0.1},
    timeout=120)

print(f"HTTP {r2.status_code}")
if r2.status_code != 200:
    print(r2.text[:400])
    sys.exit(1)

content = r2.json()["choices"][0]["message"]["content"]
print(f"Got {len(content)} chars")

with open(OUT_FILE, "w", encoding="utf-8") as f:
    f.write(f"# GROK SHARE THREAD {SHARE_ID}\n")
    f.write(f"# Fetched via Grok-4 API\n")
    f.write(f"# URL: https://x.com/i/grok/share/{SHARE_ID}\n\n")
    f.write(content)

print(f"Saved to {OUT_FILE}")
print("\n--- PREVIEW (first 2000 chars) ---")
print(content[:2000])
