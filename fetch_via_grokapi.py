#!/usr/bin/env python3
"""
Use xAI grok-4-0709 with native web_search tool to retrieve Grok share URL content.
Available models confirmed: grok-4-0709, grok-4-fast-reasoning, grok-4-1-fast-reasoning
"""
import sys, json, requests
sys.path.insert(0, '.')
from APIKeyManager import get_xai_api_key

CONV_URL = 'https://x.com/i/grok/share/ab2d0965e3a74a0da32749a7a2591dc7'
CONV_ID = 'ab2d0965e3a74a0da32749a7a2591dc7'
OUTPUT = f'grok_share_{CONV_ID}_content.txt'

api_key = get_xai_api_key()
url = 'https://api.x.ai/v1/chat/completions'
headers = {'Content-Type': 'application/json', 'Authorization': f'Bearer {api_key}'}

prompt = f"""Access this URL and return the COMPLETE verbatim text of every single message in the Grok conversation:
{CONV_URL}

I need the FULL conversation — every user message and every Grok response, word for word, unabridged. Include all physics equations, code blocks, formulas, and technical content exactly as written. Do not summarize anything."""

def call(model, extra=None):
    payload = {
        'model': model,
        'messages': [
            {'role': 'system', 'content': 'You are a research assistant with web access. When given a URL, fetch it and return the complete verbatim content.'},
            {'role': 'user', 'content': prompt}
        ],
        'temperature': 0.0
    }
    if extra:
        payload.update(extra)
    print(f'\n--- {model} {list(extra.keys()) if extra else ""} ---')
    r = requests.post(url, headers=headers, json=payload, timeout=180)
    print(f'Status: {r.status_code}')
    if r.status_code == 200:
        data = r.json()
        msg = data['choices'][0]['message']
        content = msg.get('content') or ''
        tool_calls = msg.get('tool_calls', [])
        if tool_calls:
            print(f'Tool calls: {[tc["function"]["name"] for tc in tool_calls]}')
        if content:
            print(f'Content ({len(content)} chars):\n{content[:4000]}')
            return content
    else:
        print(r.text[:400])
    return None

# Try 1: grok-4-0709 with built-in web_search tool (xAI format)
result = call('grok-4-0709', {'tools': [{'type': 'web_search'}]})

# Try 2: grok-4-0709 direct (no tools)
if not result or len(result) < 200:
    result = call('grok-4-0709')

# Try 3: grok-4-fast-reasoning direct
if not result or len(result) < 200:
    result = call('grok-4-fast-reasoning')

if result and len(result) > 300:
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write(result)
    print(f'\n✓ Saved: {OUTPUT} ({len(result)} bytes)')
else:
    print('\nAll attempts returned insufficient content.')
