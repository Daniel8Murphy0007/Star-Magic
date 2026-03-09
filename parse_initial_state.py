#!/usr/bin/env python3
"""
Extract Grok conversation from window.__INITIAL_STATE__ embedded in page textContent.
"""
import re, json, os

SHARE_ID = '7b78ffcb915f48bb90d55034c9c50b18'
INPUT    = f'grok_share_{SHARE_ID}_content.txt'
OUTPUT   = f'grok_share_{SHARE_ID}_parsed.txt'

with open(INPUT, 'r', encoding='utf-8') as f:
    raw = f.read()

print(f'Input: {len(raw)} chars')

# Extract window.__INITIAL_STATE__ JSON
match = re.search(r'window\.__INITIAL_STATE__=(\{.+)', raw, re.DOTALL)
if not match:
    print('ERROR: __INITIAL_STATE__ not found')
    # Dump a sample
    idx = raw.find('INITIAL_STATE')
    print('Context around INITIAL_STATE:', raw[max(0, idx-50):idx+200])
    exit(1)

state_str = match.group(1)
# The JSON ends at the next </script> or window. assignment
ended = re.search(r'\s*;?\s*window\.|</script>', state_str)
if ended:
    state_str = state_str[:ended.start()]

print(f'State JSON: {len(state_str)} chars')

try:
    state = json.loads(state_str)
except json.JSONDecodeError as e:
    print(f'JSON parse error: {e}')
    # Try to find and fix truncation
    print('First 500 chars of state_str:', state_str[:500])
    exit(1)

# Walk the state to find Grok conversation content
def find_grok_content(obj, path='', depth=0):
    """Recursively search for Grok conversation messages."""
    results = []
    if depth > 20:
        return results
    if isinstance(obj, dict):
        for k, v in obj.items():
            subpath = f'{path}.{k}'
            # Keys that might contain conversation content
            if k in ('message', 'text', 'content', 'response', 'body',
                     'humanMessage', 'assistantMessage', 'conversationTitle',
                     'snippet', 'grokShare', 'conversation', 'messages',
                     'grokConversation'):
                if isinstance(v, str) and len(v) > 30:
                    results.append((subpath, v))
                elif isinstance(v, (dict, list)):
                    results.extend(find_grok_content(v, subpath, depth+1))
            else:
                results.extend(find_grok_content(v, subpath, depth+1))
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            results.extend(find_grok_content(item, f'{path}[{i}]', depth+1))
    return results

print('\nTop-level keys in state:', list(state.keys())[:20])

# Look specifically in entities, grokConversation, etc.
grok_sections = {k: v for k, v in state.items() if 'grok' in k.lower() or 'conv' in k.lower()}
print('Grok/conversation sections:', list(grok_sections.keys()))

found = find_grok_content(state)
print(f'Content items found: {len(found)}')

output_parts = []
for path, text in found:
    output_parts.append(f'[{path}]\n{text}')
    print(f'\nPATH: {path}')
    print(text[:300])
    print('---')

if output_parts:
    full = '\n\n'.join(output_parts)
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write(full)
    print(f'\nSaved {len(found)} content items -> {OUTPUT}')
else:
    print('\nNo content items found in state. Checking all string values > 100 chars...')
    def all_long_strings(obj, path='', depth=0):
        r = []
        if depth > 15: return r
        if isinstance(obj, dict):
            for k, v in obj.items():
                r.extend(all_long_strings(v, f'{path}.{k}', depth+1))
        elif isinstance(obj, list):
            for i, v in enumerate(obj):
                r.extend(all_long_strings(v, f'{path}[{i}]', depth+1))
        elif isinstance(obj, str) and len(obj) > 100:
            r.append((path, obj))
        return r
    
    longs = all_long_strings(state)
    print(f'Long strings in state: {len(longs)}')
    for path, text in longs[:10]:
        print(f'{path}: {text[:100]}')
