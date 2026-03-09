#!/usr/bin/env python3
"""
Fetch Grok share thread: 7b78ffcb915f48bb90d55034c9c50b18
Direct GraphQL API call using Firefox X.com auth cookies — no Selenium needed.
"""
import os, sys, json, sqlite3, shutil, tempfile
import urllib.request, urllib.parse

SHARE_ID    = '7b78ffcb915f48bb90d55034c9c50b18'
OUTPUT_FILE = f'grok_share_{SHARE_ID}_content.txt'
# X.com's public app-level Bearer token (from their web JS bundle)
BEARER      = 'AAAAAAAAAAAAAAAAAAAAANRILgAAAAAAnNwIzUejRCOuH5E6I38u8zUk38Y%3DUouWPwtxNUezdzobearsAHsJWnk5GbOf91kX0yn2Ez44MUA'
# GrokShare GraphQL query ID (discovered via CDP monitoring)
QUERY_ID    = '3aCm_HRrYXX8T7sas50Zlw'

FF_PROFILE   = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
FF_COOKIES   = os.path.join(FF_PROFILE, 'cookies.sqlite')

# ── Step 1: pull auth_token + ct0 from Firefox ───────────────────────────────
print('Reading X.com auth cookies from Firefox...')
tmp_db = os.path.join(tempfile.gettempdir(), 'ff_ck_direct.db')
shutil.copy2(FF_COOKIES, tmp_db)

auth_token = ct0 = None
try:
    conn = sqlite3.connect(tmp_db)
    cur  = conn.cursor()
    cur.execute("""
        SELECT name, value FROM moz_cookies
        WHERE (host LIKE '%.x.com%' OR host = 'x.com')
          AND name IN ('auth_token', 'ct0')
    """)
    for name, value in cur.fetchall():
        if name == 'auth_token':
            auth_token = value
        elif name == 'ct0':
            ct0 = value
    conn.close()
finally:
    try: os.unlink(tmp_db)
    except: pass

if not auth_token or not ct0:
    print(f'ERROR: Missing cookies — auth_token={bool(auth_token)}, ct0={bool(ct0)}')
    sys.exit(1)

print(f'auth_token: {auth_token[:12]}... ct0: {ct0[:12]}...')

# ── Step 1b: get current Bearer from X.com JS bundle ─────────────────────────
import re as _re
print('Fetching current Bearer token from X.com JS bundle...')
bearer_found = None
try:
    main_req = urllib.request.Request(
        'https://x.com/',
        headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/134.0.0.0 Safari/537.36'}
    )
    with urllib.request.urlopen(main_req, timeout=15) as r:
        html = r.read().decode('utf-8', errors='replace')
    # Find JS bundle URLs
    bundle_urls = _re.findall(r'src="(https://abs\.twimg\.com/responsive-web/client-web/[^"]+\.js)"', html)
    if not bundle_urls:
        bundle_urls = _re.findall(r'"(https://abs\.twimg\.com/responsive-web/client-web/main\.[^"]+\.js)"', html)
    print(f'Found {len(bundle_urls)} JS bundle URLs')
    # Check first few bundles for Bearer pattern
    for burl in bundle_urls[:8]:
        try:
            br = urllib.request.Request(burl, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(br, timeout=10) as r2:
                js_chunk = r2.read(200_000).decode('utf-8', errors='replace')
            # Bearer tokens are ~100 char base64 strings starting with AAAA
            match = _re.search(r'(AAAAAAAAAAAAA[A-Za-z0-9%/+=]{50,120})', js_chunk)
            if match:
                bearer_found = match.group(1)
                print(f'Bearer found in {burl.split("/")[-1]}: {bearer_found[:20]}...')
                break
        except Exception:
            continue
except Exception as e:
    print(f'Could not fetch Bearer from JS bundle: {e}')

BEARER_USE = urllib.parse.unquote(bearer_found) if bearer_found else urllib.parse.unquote(BEARER)
print(f'Using Bearer: {BEARER_USE[:20]}...')

# ── Step 2: build GraphQL request ────────────────────────────────────────────
variables = json.dumps({"grok_share_id": SHARE_ID}, separators=(',', ':'))
features  = json.dumps({
    "premium_content_api_read_enabled": False,
    "communities_web_enable_tweet_community_results_fetch": True
}, separators=(',', ':'))

url = (f'https://x.com/i/api/graphql/{QUERY_ID}/GrokShare'
       f'?variables={urllib.parse.quote(variables)}'
       f'&features={urllib.parse.quote(features)}')

cookie_str = f'auth_token={auth_token}; ct0={ct0}'
headers = {
    'Authorization':             f'Bearer {BEARER_USE}',
    'Cookie':                    cookie_str,
    'x-csrf-token':              ct0,
    'Content-Type':              'application/json',
    'Accept':                    '*/*',
    'User-Agent':                'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/134.0.0.0 Safari/537.36',
    'Referer':                   f'https://x.com/i/grok/share/{SHARE_ID}',
    'x-twitter-active-user':     'yes',
    'x-twitter-client-language': 'en',
    'Accept-Language':           'en-US,en;q=0.9',
}

print(f'Calling GrokShare API for share ID: {SHARE_ID}')
req = urllib.request.Request(url, headers=headers)

try:
    with urllib.request.urlopen(req, timeout=30) as resp:
        status = resp.status
        raw    = resp.read().decode('utf-8')
        print(f'HTTP {status} — {len(raw)} bytes received')
except urllib.error.HTTPError as e:
    raw    = e.read().decode('utf-8')
    status = e.code
    print(f'HTTP {status} error body: {raw[:500]}')
    sys.exit(1)

# ── Step 3: parse and extract conversation ───────────────────────────────────
try:
    data = json.loads(raw)
except json.JSONDecodeError:
    print('Response is not JSON:')
    print(raw[:1000])
    sys.exit(1)

# Save raw JSON for inspection
raw_file = f'grok_share_{SHARE_ID}_raw.json'
with open(raw_file, 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=2, ensure_ascii=False)
print(f'Raw JSON saved → {raw_file}')

# Navigate the response tree
def extract_text(obj, depth=0):
    """Recursively extract all text/message content from the JSON structure."""
    lines = []
    if isinstance(obj, dict):
        # Common keys that hold conversation content
        for key in ('message', 'text', 'content', 'response', 'userMessage',
                    'humanMessage', 'assistantMessage', 'result', 'body',
                    'snippet', 'conversationTitle'):
            if key in obj and isinstance(obj[key], str) and len(obj[key]) > 10:
                lines.append(f'[{key}] {obj[key]}')
        for v in obj.values():
            lines.extend(extract_text(v, depth+1))
    elif isinstance(obj, list):
        for item in obj:
            lines.extend(extract_text(item, depth+1))
    return lines

texts = extract_text(data)

if texts:
    output = '\n\n'.join(texts)
    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(output)
    print(f'\nConversation extracted → {OUTPUT_FILE}')
    print(f'Total content blocks: {len(texts)}')
    # Preview first 2000 chars
    print('\n' + '='*60)
    print(output[:2000])
    print('='*60)
else:
    # No recognised structure — print top-level keys for inspection
    print('\nNo conversation text found. Top-level response keys:')
    if isinstance(data, dict):
        def show_keys(d, prefix=''):
            for k, v in d.items():
                t = type(v).__name__
                preview = str(v)[:80] if not isinstance(v, (dict, list)) else ''
                print(f'  {prefix}{k} ({t}) {preview}')
                if isinstance(v, dict) and len(v) < 10:
                    show_keys(v, prefix + '  ')
        show_keys(data)
    print(f'\nFull raw JSON in {raw_file}')
