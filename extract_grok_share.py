#!/usr/bin/env python3
"""
Extract Grok share conversation content from the publicly accessible HTML.
Parses __NEXT_DATA__ / window.__STATE__ JSON blobs embedded in the page.
"""
import requests
import json
import re

SHARE_ID  = '7b0e961fa19846bea8aed946b0650e93'
SHARE_URL = f'https://x.com/i/grok/share/{SHARE_ID}'

HEADERS = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.9',
}

print(f'Fetching: {SHARE_URL}')
r = requests.get(SHARE_URL, timeout=15, headers=HEADERS)
print(f'Status: {r.status_code}  |  Size: {len(r.text):,} chars')
html = r.text

# ── Strategy 1: Extract __NEXT_DATA__ (Next.js SSR state) ────────────────────
print('\n--- Strategy 1: __NEXT_DATA__ ---')
m = re.search(r'<script id="__NEXT_DATA__"[^>]*>(.*?)</script>', html, re.DOTALL)
if m:
    try:
        nd = json.loads(m.group(1))
        print(f'Found __NEXT_DATA__  ({len(m.group(1)):,} chars)')
        # Pretty-print top-level keys
        print(f'Top keys: {list(nd.keys())}')
        # Recursively search for conversation/messages
        text = json.dumps(nd, indent=2)
        # Save full JSON for inspection
        with open('grok_share_next_data.json', 'w', encoding='utf-8') as f:
            f.write(text)
        print('Saved: grok_share_next_data.json')
    except json.JSONDecodeError as e:
        print(f'JSON parse error: {e}')
else:
    print('Not found.')

# ── Strategy 2: Search for embedded JSON blobs with conversation content ──────
print('\n--- Strategy 2: Inline JSON blobs ---')
blobs = re.findall(r'<script[^>]*>\s*window\.__(?:STATE|INITIAL_STATE|DATA)__\s*=\s*(\{.*?\});\s*</script>', html, re.DOTALL)
print(f'Found {len(blobs)} window.__STATE__ blobs')
for i, blob in enumerate(blobs):
    try:
        d = json.loads(blob)
        print(f'  Blob {i}: top keys = {list(d.keys())[:8]}')
    except Exception as e:
        print(f'  Blob {i}: parse error {e}  (len={len(blob)})')

# ── Strategy 3: Extract all <script> JSON blocks and search for messages ──────
print('\n--- Strategy 3: Scan all script tags for conversation text ---')
scripts = re.findall(r'<script[^>]*>(.*?)</script>', html, re.DOTALL)
print(f'Total <script> tags: {len(scripts)}')
for i, s in enumerate(scripts):
    stripped = s.strip()
    if len(stripped) < 20:
        continue
    if any(kw in stripped for kw in ['grok_share', 'conversationId', 'messages', 'grokShare', 'shareId']):
        print(f'\n  >> Script {i} contains conversation keywords (len={len(stripped):,}):')
        print(f'     {stripped[:400]}')

# ── Strategy 4: Look for og/meta tags summarising the content ────────────────
print('\n--- Strategy 4: Open Graph / meta description ---')
og_tags = re.findall(r'<meta[^>]+(?:og:|twitter:|name="description")[^>]*/>', html, re.IGNORECASE)
for tag in og_tags[:10]:
    print(f'  {tag[:200]}')

# ── Strategy 5: Look for any text that looks like conversation turns ──────────
print('\n--- Strategy 5: Raw text extraction (visible content) ---')
# Strip all tags
clean = re.sub(r'<[^>]+>', ' ', html)
clean = re.sub(r'\s+', ' ', clean).strip()
# Look for long non-whitespace runs that are likely content
# Find substrings around the share ID
idx = clean.find(SHARE_ID[:8])
if idx >= 0:
    print(f'Share ID found in text at pos {idx}: ...{clean[max(0,idx-100):idx+300]}...')
else:
    print('Share ID not found in stripped text.')
    # Just print a window of what looks like body content
    # Find first long text run (likely React-rendered text)
    matches = re.findall(r'[A-Za-z][A-Za-z0-9 ,.\-!?()]{80,}', clean)
    print(f'Long text runs found: {len(matches)}')
    for m in matches[:5]:
        print(f'  "{m[:120]}"')

print('\nDone.')
