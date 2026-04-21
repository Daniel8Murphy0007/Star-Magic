#!/usr/bin/env python3
"""
Direct HTTP approach to fetch Grok share page.
Reads Firefox cookies, sends proper browser headers, parses conversation content.
No Selenium - avoids all automation detection.
"""
import sqlite3, shutil, os, re, json, html, gzip
import requests

SHARE_ID = '7b78ffcb915f48bb90d55034c9c50b18'
SHARE_URL = f'https://x.com/i/grok/share/{SHARE_ID}'
FF_PROFILE = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
COOKIES_DB  = os.path.join(FF_PROFILE, 'cookies.sqlite')
TMP_DB      = os.path.join(os.environ['TEMP'], 'ff_cookies_copy.sqlite')

print('=== Direct HTTP Grok Share Fetcher ===')
print(f'Target: {SHARE_URL}')

# ── 1. Read Firefox cookies ──────────────────────────────────────────────────
shutil.copy2(COOKIES_DB, TMP_DB)
conn = sqlite3.connect(TMP_DB)
rows = conn.execute(
    "SELECT name, value, host FROM moz_cookies WHERE host LIKE '%x.com%'"
).fetchall()
conn.close()

jar = requests.cookies.RequestsCookieJar()
auth_token = ct0 = None
for name, value, host in rows:
    domain = host.lstrip('.')
    jar.set(name, value, domain=domain)
    if name == 'auth_token': auth_token = value
    if name == 'ct0':        ct0 = value

print(f'Cookies loaded: {len(rows)}')
print(f'auth_token: {"YES (" + auth_token[:12] + "...)" if auth_token else "MISSING"}')
print(f'ct0:        {"YES (" + ct0[:12] + "...)" if ct0 else "MISSING"}')

if not auth_token or not ct0:
    print('ERROR: Missing auth cookies'); exit(1)

# ── 2. Build headers that look like a real browser ───────────────────────────
headers = {
    'User-Agent':      'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:136.0) Gecko/20100101 Firefox/136.0',
    'Accept':          'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.5',
    'Accept-Encoding': 'gzip, deflate, br',
    'DNT':             '1',
    'Upgrade-Insecure-Requests': '1',
    'Sec-Fetch-Dest':  'document',
    'Sec-Fetch-Mode':  'navigate',
    'Sec-Fetch-Site':  'none',
    'Sec-Fetch-User':  '?1',
    'Connection':      'keep-alive',
    'X-Csrf-Token':    ct0,
    'Authorization':   f'Bearer [PROMPT_FOR_GROK_WEB_BEARER_TOKEN]',
}

session = requests.Session()
session.cookies = jar
session.headers.update(headers)

# ── 3. GET the share page ─────────────────────────────────────────────────────
print(f'\nGET {SHARE_URL} ...')
resp = session.get(SHARE_URL, timeout=30, allow_redirects=True)
print(f'Status: {resp.status_code}  Content-Length: {len(resp.content)}')
print(f'Content-Type: {resp.headers.get("content-type", "?")}')

html_text = resp.text
print(f'HTML text: {len(html_text)} chars')

# ── 4. Look for conversation content ─────────────────────────────────────────
def try_extract_state(src):
    m = re.search(r'window\.__INITIAL_STATE__\s*=\s*(\{.+?);\s*(?:window\.|</script>|$)', src, re.DOTALL)
    if not m:
        m = re.search(r'window\.__INITIAL_STATE__\s*=\s*(\{.*)', src, re.DOTALL)
    if not m:
        return None
    state_str = m.group(1).rstrip(';').strip()
    try:
        return json.loads(state_str)
    except json.JSONDecodeError as e:
        print(f'JSON parse error: {e}')
        return None

state = try_extract_state(html_text)
if state:
    print('Found __INITIAL_STATE__')
    print('Top-level keys:', list(state.keys())[:15])
    grok_keys = {k: v for k, v in state.items() if 'grok' in k.lower() or 'conv' in k.lower()}
    print('Grok keys:', list(grok_keys.keys()))
    
    # Look for any non-empty entities
    entities = state.get('entities', {})
    non_empty = {k: v for k, v in entities.items() if v}
    print('Non-empty entities:', list(non_empty.keys()))
else:
    print('No __INITIAL_STATE__ found')

# ── 5. Try GrokShare API directly ────────────────────────────────────────────
print('\n--- Trying GrokShare API ---')
api_headers = {
    'User-Agent':       'Mozilla/5.0 (Windows NT 10.0; Win64; x64; rv:136.0) Gecko/20100101 Firefox/136.0',
    'Accept':           '*/*',
    'Accept-Language':  'en-US,en;q=0.5',
    'Accept-Encoding':  'gzip, deflate, br',
    'Content-Type':     'application/json',
    'X-Twitter-Auth-Type': 'OAuth2Session',
    'X-Twitter-Client-Language': 'en',
    'X-Csrf-Token':     ct0,
    'Authorization':    'Bearer [PROMPT_FOR_GROK_WEB_BEARER_TOKEN]',
    'Referer':          f'https://x.com/i/grok/share/{SHARE_ID}',
    'Origin':           'https://x.com',
    'Connection':       'keep-alive',
    'Sec-Fetch-Dest':   'empty',
    'Sec-Fetch-Mode':   'cors',
    'Sec-Fetch-Site':   'same-origin',
}

params = {
    'variables': json.dumps({'shareId': SHARE_ID}),
    'features':  json.dumps({'responsive_web_graphql_exclude_directive_enabled': True,
                              'verified_phone_label_enabled': False})
}

graphql_url = 'https://x.com/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare'
api_resp = session.get(graphql_url, headers=api_headers, params=params, timeout=30)
print(f'GrokShare API status: {api_resp.status_code}')
print(f'GrokShare API body ({len(api_resp.text)} chars):')
print(api_resp.text[:2000])

# ── 6. Save everything ────────────────────────────────────────────────────────
with open(f'grok_direct_{SHARE_ID}_html.txt', 'w', encoding='utf-8') as f:
    f.write(html_text)
print(f'\nSaved HTML -> grok_direct_{SHARE_ID}_html.txt')

if api_resp.text and len(api_resp.text) > 20:
    with open(f'grok_direct_{SHARE_ID}_api.json', 'w', encoding='utf-8') as f:
        f.write(api_resp.text)
    print(f'Saved API response -> grok_direct_{SHARE_ID}_api.json')
