#!/usr/bin/env python3
"""Extract Grok conversation from the successfully fetched HTML."""
import re, json

SHARE_ID = '7b78ffcb915f48bb90d55034c9c50b18'
HTML_FILE = f'grok_direct_{SHARE_ID}_html.txt'

with open(HTML_FILE, 'r', encoding='utf-8') as f:
    html_text = f.read()

print(f'HTML: {len(html_text)} chars')

# Extract __INITIAL_STATE__
m = re.search(r'window\.__INITIAL_STATE__\s*=\s*(\{.+)', html_text, re.DOTALL)
state_str = m.group(1)
# Find the end: next window. assignment or </script>
end = re.search(r'\s*;\s*(?:window\.|</script>)', state_str)
if end:
    state_str = state_str[:end.start()]

state = json.loads(state_str)
print(f'State parsed. Top-level keys: {list(state.keys())}')

entities = state.get('entities', {})
grok_share = entities.get('grokShare', {})
print(f'grokShare keys: {list(grok_share.keys())[:10]}')

# ── Pretty-print the full grokShare entity ──────────────────────────────────
OUTPUT = f'grok_conversation_{SHARE_ID}.txt'
lines = []

def extract_messages(obj, indent=0):
    pad = '  ' * indent
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, str) and len(v) > 20:
                lines.append(f'{pad}[{k}]: {v}')
            elif isinstance(v, (dict, list)):
                lines.append(f'{pad}[{k}]:')
                extract_messages(v, indent+1)
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            lines.append(f'{pad}--- item {i} ---')
            extract_messages(item, indent+1)

extract_messages(grok_share)

# Also dump raw JSON for full fidelity
raw_json = json.dumps(grok_share, indent=2)

with open(OUTPUT, 'w', encoding='utf-8') as f:
    f.write('=== GROK SHARE CONVERSATION (Extracted) ===\n\n')
    f.write('\n'.join(lines))
    f.write('\n\n=== RAW JSON ===\n\n')
    f.write(raw_json)

print(f'Written: {OUTPUT} ({len(raw_json)} chars)')
print('\n=== PREVIEW (first 3000 chars) ===')
print('\n'.join(lines[:80]))
