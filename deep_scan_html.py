#!/usr/bin/env python3
"""Deep-scan the HTML for any Grok conversation content."""
import re, json, html as html_module

SHARE_ID = '7b78ffcb915f48bb90d55034c9c50b18'
HTML_FILE = f'grok_direct_{SHARE_ID}_html.txt'

with open(HTML_FILE, 'r', encoding='utf-8') as f:
    src = f.read()

print(f'HTML: {len(src)} chars')

# ── 1. Find ALL script tags ──────────────────────────────────────────────────
scripts = re.findall(r'<script[^>]*>(.*?)</script>', src, re.DOTALL | re.IGNORECASE)
print(f'Script blocks: {len(scripts)}')

# ── 2. Look for UQFF, physics, conversation clues ─────────────────────────
keywords = ['UQFF', 'uqff', 'grok', 'conversation', 'message', 'humanMessage', 
            'assistantMessage', 'shareConversation', 'share_conversation',
            'conversationId', 'shareId', 'grokResponse', 'result']
for kw in keywords:
    hits = [(m.start(), src[max(0,m.start()-30):m.start()+80]) for m in re.finditer(kw, src)]
    if hits:
        print(f'\n[{kw}] ({len(hits)} hits):')
        for pos, ctx in hits[:3]:
            print(f'  @{pos}: {ctx!r}')

# ── 3. Look for data-* attributes with JSON ──────────────────────────────────
data_attrs = re.findall(r'data-(?:state|store|props|initial)[^=]*=["\'](.*?)["\']', src)
print(f'\ndata-state/store/props attrs: {len(data_attrs)}')
for d in data_attrs[:5]:
    print(d[:200])

# ── 4. Other JSON embeds (next.js style / react SSR) ────────────────────────
json_embeds = re.findall(r'__NEXT_DATA__|__REDUX_STATE__|__APP_STATE__|id="__NEXT_DATA"', src)
print(f'\nNext.js / Redux embeds: {len(json_embeds)}')

# ── 5. Look for the GrokShare render output ─────────────────────────────────
# Twitter/X renders grok shares in a React tree seeded from __INITIAL_STATE__
# The grokShare.entities path should have the conversation
m = re.search(r'window\.__INITIAL_STATE__\s*=\s*(\{.+)', src, re.DOTALL)
if m:
    state_raw = m.group(1)
    end = re.search(r'\s*;\s*(?:window\.|</script>)', state_raw)
    if end:
        state_raw = state_raw[:end.start()]
    state = json.loads(state_raw)
    entities = state.get('entities', {})
    
    # Dump ALL non-empty entities
    print('\n--- All entities ---')
    for k, v in entities.items():
        print(f'  {k}: {type(v).__name__} len={len(v) if isinstance(v, (dict,list)) else len(str(v))} | sample={str(v)[:80]}')
    
    # Specifically check conversations
    convs = entities.get('conversations', {})
    print(f'\nconversations: {list(convs.keys())[:5]}')
    for cid, cdata in list(convs.items())[:2]:
        print(f'  [{cid}]: {json.dumps(cdata)[:300]}')
    
    # Check entries
    entries = entities.get('entries', {})
    print(f'\nentries count: {len(entries)}')
    for eid, edata in list(entries.items())[:3]:
        print(f'  [{eid}]: {json.dumps(edata)[:300]}')

# ── 6. Any inline JSON containing the share ID ──────────────────────────────
print(f'\n--- Searching for share ID in HTML ---')
share_id_hits = [m.start() for m in re.finditer(SHARE_ID, src)]
print(f'Share ID appears {len(share_id_hits)} times')
for pos in share_id_hits[:5]:
    print(f'  @{pos}: {src[max(0,pos-60):pos+100]!r}')

# ── 7. Check the og: meta tags rendered by SSR ──────────────────────────────
og_tags = re.findall(r'<meta[^>]+(?:og:|twitter:)[^>]*>', src)
print(f'\nOG/Twitter meta tags: {len(og_tags)}')
for tag in og_tags[:10]:
    print(tag)
