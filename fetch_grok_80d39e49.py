#!/usr/bin/env python3
"""
Fetch Grok share URL 80d39e496f3a49bc9106276bd5a4e400 using xAI API
"""
import json
import requests

CONV_URL = 'https://x.com/i/grok/share/80d39e496f3a49bc9106276bd5a4e400'
CONV_ID = '80d39e496f3a49bc9106276bd5a4e400'
OUTPUT = f'grok_share_{CONV_ID}_content.txt'
API_KEY = 'xai-REDACTED_KEY'

url = 'https://api.x.ai/v1/chat/completions'
headers = {'Content-Type': 'application/json', 'Authorization': f'Bearer {API_KEY}'}

prompt = f"""Access this Grok conversation URL and return the COMPLETE VERBATIM text of EVERY message:
{CONV_URL}

I need the FULL conversation — every user message and every Grok response, word for word, unabridged.
Include ALL physics equations, calculator code, UQFF formulas, astrophysical data, validation results, 
technical content, mathematical derivations, and code blocks EXACTLY as written.

DO NOT SUMMARIZE ANYTHING. Return the complete raw transcript."""

print("=" * 80)
print(f"Fetching Grok thread: {CONV_ID}")
print(f"Target URL: {CONV_URL}")
print("=" * 80)

payload = {
    'model': 'grok-4-1-fast-reasoning',
    'messages': [
        {'role': 'system', 'content': 'You are a research assistant. When given a Grok share URL, retrieve and return the complete verbatim conversation transcript with all technical content.'},
        {'role': 'user', 'content': prompt}
    ],
    'temperature': 0.0,
    'max_tokens': 16000  # Get maximum possible content
}

print(f'\n--- Calling xAI API (model: grok-4-1-fast-reasoning) ---')
r = requests.post(url, headers=headers, json=payload, timeout=180)
print(f'Status: {r.status_code}')

if r.status_code == 200:
    data = r.json()
    content = data['choices'][0]['message'].get('content', '')
    
    if content:
        print(f'Content length: {len(content)} chars')
        print(f'\nFirst 2000 chars:\n{content[:2000]}')
        print(f'\n...\n')
        print(f'Last 500 chars:\n{content[-500:]}')
        
        with open(OUTPUT, 'w', encoding='utf-8') as f:
            f.write(content)
        
        print(f'\n{"=" * 80}')
        print(f'✓ SUCCESS: Saved to {OUTPUT}')
        print(f'✓ Content size: {len(content):,} bytes')
        print("=" * 80)
    else:
        print('ERROR: Empty response from API')
else:
    print(f'ERROR: HTTP {r.status_code}')
    print(r.text[:400])
