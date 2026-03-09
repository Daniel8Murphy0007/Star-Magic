#!/usr/bin/env python3
"""
Use xAI grok-4-1-fast-reasoning with native web_search tool to retrieve Grok share URL content.
Target: https://x.com/i/grok/share/6bfa709582e34d4c9f78ecd14352d54f
"""
import sys, json, requests
sys.path.insert(0, '.')
from APIKeyManager import get_xai_api_key

CONV_URL = 'https://x.com/i/grok/share/6bfa709582e34d4c9f78ecd14352d54f'
CONV_ID = '6bfa709582e34d4c9f78ecd14352d54f'
OUTPUT = f'grok_share_{CONV_ID}_content.txt'

api_key = get_xai_api_key()
url = 'https://api.x.ai/v1/chat/completions'
headers = {'Content-Type': 'application/json', 'Authorization': f'Bearer {api_key}'}

prompt = f"""Access this URL and return the COMPLETE verbatim text of every single message in the Grok conversation:
{CONV_URL}

I need the FULL conversation — every user message and every Grok response, word for word, unabridged. Include all physics equations, code blocks, formulas, and technical content exactly as written. Do not summarize anything. Return the entire conversation transcript."""

def call(model, extra=None):
    payload = {
        'model': model,
        'messages': [
            {'role': 'system', 'content': 'You are a research assistant with web access. When given a URL, fetch it and return the complete verbatim content without any summarization.'},
            {'role': 'user', 'content': prompt}
        ],
        'temperature': 0.0
    }
    if extra:
        payload.update(extra)
    print(f'\n--- MODEL: {model} ---')
    if extra:
        print(f'Extra params: {list(extra.keys())}')
    
    try:
        r = requests.post(url, headers=headers, json=payload, timeout=300)
        print(f'Status: {r.status_code}')
        
        if r.status_code == 200:
            data = r.json()
            msg = data['choices'][0]['message']
            content = msg.get('content') or ''
            tool_calls = msg.get('tool_calls', [])
            
            if tool_calls:
                print(f'Tool calls: {[tc["function"]["name"] for tc in tool_calls]}')
            
            if content:
                print(f'Content length: {len(content)} chars')
                print(f'First 2000 chars:\n{content[:2000]}')
                print(f'...\nLast 500 chars:\n{content[-500:]}')
                return content
        else:
            print(f'Error response: {r.text[:500]}')
    except Exception as e:
        print(f'Exception: {e}')
    
    return None

print('=' * 80)
print(f'Fetching Grok thread: {CONV_ID}')
print(f'Target URL: {CONV_URL}')
print('=' * 80)

# Try 1: grok-4-1-fast-reasoning with web_search (recommended model)
result = call('grok-4-1-fast-reasoning', {'tools': [{'type': 'web_search'}]})

# Try 2: grok-4-1-fast-reasoning direct (no tools, may use built-in reasoning)
if not result or len(result) < 500:
    print('\nAttempt 1 insufficient, trying direct call without tools...')
    result = call('grok-4-1-fast-reasoning')

# Try 3: grok-4-0709 with web_search (fallback)
if not result or len(result) < 500:
    print('\nAttempt 2 insufficient, trying grok-4-0709 with web_search...')
    result = call('grok-4-0709', {'tools': [{'type': 'web_search'}]})

# Try 4: grok-4-fast-reasoning (fallback)
if not result or len(result) < 500:
    print('\nAttempt 3 insufficient, trying grok-4-fast-reasoning...')
    result = call('grok-4-fast-reasoning')

if result and len(result) > 500:
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write(result)
    print(f'\n{"=" * 80}')
    print(f'✓ SUCCESS: Saved to {OUTPUT}')
    print(f'✓ Content size: {len(result):,} bytes')
    print(f'{"=" * 80}')
else:
    print(f'\n{"=" * 80}')
    print('✗ FAILURE: All attempts returned insufficient content.')
    print('The Grok share may require authentication or may not be publicly accessible.')
    print(f'{"=" * 80}')
