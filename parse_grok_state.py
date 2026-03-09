#!/usr/bin/env python3
"""
Parse the window.__INITIAL_STATE__ blob from the Grok share page.
The blob may be truncated by regex greediness — use a JS-aware extractor.
"""
import requests
import json
import re

SHARE_ID  = '7b0e961fa19846bea8aed946b0650e93'
SHARE_URL = f'https://x.com/i/grok/share/{SHARE_ID}'

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
}

r = requests.get(SHARE_URL, timeout=15, headers=HEADERS)
html = r.text
print(f'Page size: {len(html):,} chars')

# ── Extract the raw __INITIAL_STATE__ assignment text ────────────────────────
marker = 'window.__INITIAL_STATE__='
start = html.find(marker)
if start == -1:
    print('ERROR: window.__INITIAL_STATE__ not found')
    exit(1)

# Find the opening brace
brace_start = html.index('{', start)

# Walk the string counting braces to find matching close
depth = 0
i = brace_start
in_string = False
escape = False
for i, ch in enumerate(html[brace_start:], brace_start):
    if escape:
        escape = False
        continue
    if ch == '\\' and in_string:
        escape = True
        continue
    if ch == '"' and not escape:
        in_string = not in_string
        continue
    if not in_string:
        if ch == '{':
            depth += 1
        elif ch == '}':
            depth -= 1
            if depth == 0:
                break

raw_json = html[brace_start:i+1]
print(f'Extracted JSON blob: {len(raw_json):,} chars')

try:
    state = json.loads(raw_json)
    print('JSON parse: OK')
except json.JSONDecodeError as e:
    print(f'JSON parse ERROR at char {e.pos}: {e.msg}')
    print(f'Context: ...{raw_json[max(0,e.pos-80):e.pos+80]}...')
    exit(1)

# Save full state
with open('grok_share_initial_state.json', 'w', encoding='utf-8') as f:
    json.dump(state, f, indent=2, ensure_ascii=False)
print('Saved: grok_share_initial_state.json')

# ── Navigate the state for Grok share content ────────────────────────────────
print('\n=== TOP-LEVEL KEYS ===')
print(list(state.keys()))

def find_keys(obj, target_keys, path='', depth=0, results=None):
    """Recursively find values for any of the target keys."""
    if results is None:
        results = []
    if depth > 10:
        return results
    if isinstance(obj, dict):
        for k, v in obj.items():
            current_path = f'{path}.{k}'
            if any(tk.lower() in k.lower() for tk in target_keys):
                results.append((current_path, v))
            find_keys(v, target_keys, current_path, depth+1, results)
    elif isinstance(obj, list):
        for idx, item in enumerate(obj[:20]):
            find_keys(item, target_keys, f'{path}[{idx}]', depth+1, results)
    return results

# Search for conversation / message / grokShare content
target = ['grok', 'share', 'message', 'conversation', 'turn', 'response', 'query']
hits = find_keys(state, target)

print(f'\n=== FOUND {len(hits)} MATCHING KEYS ===')
for path, val in hits[:40]:
    val_preview = json.dumps(val)[:200] if not isinstance(val, (str, int, float, bool, type(None))) else str(val)[:200]
    print(f'  {path}:')
    print(f'    {val_preview}')
    print()
