#!/usr/bin/env python3
"""
Fetch Grok share 7f90683d979346758875e4a8e6e0ee9d
Reads Edge cookies → direct API call (no Selenium needed).
"""
import os, sys, json, sqlite3, shutil, tempfile, glob, time
import urllib.request, urllib.error

SHARE_ID = '7f90683d979346758875e4a8e6e0ee9d'
OUT_TXT  = f'grok_share_{SHARE_ID}_content.txt'

# ── 1. Read Edge cookies ─────────────────────────────────────────────────────
edge_base = os.path.expandvars(r'%LOCALAPPDATA%\Microsoft\Edge\User Data')
profiles  = ['Default'] + [os.path.basename(p) for p in
             glob.glob(os.path.join(edge_base, 'Profile *'))]

auth_token = ct0 = None
for prof in profiles:
    db = os.path.join(edge_base, prof, 'Network', 'Cookies')
    if not os.path.exists(db):
        continue
    tmp = os.path.join(tempfile.gettempdir(), f'edge_{prof}.db')
    try:
        shutil.copy2(db, tmp)
        conn = sqlite3.connect(tmp)
        # Edge stores cookies encrypted — read name only to find relevant cookie
        c = conn.cursor()
        c.execute("SELECT host_key, name, encrypted_value FROM cookies "
                  "WHERE host_key LIKE '%x.com%' AND name IN ('auth_token','ct0')")
        rows = c.fetchall()
        conn.close()
        print(f'{prof}: found {len(rows)} relevant cookies (encrypted)')
        for r in rows:
            print(f'  {r[0]} :: {r[1]}')
    except Exception as e:
        print(f'Edge {prof} error: {e}')
    finally:
        try: os.unlink(tmp)
        except: pass

# ── 2. Try Firefox cookies (unencrypted) ─────────────────────────────────────
ff_profiles = glob.glob(os.path.expandvars(
    r'%APPDATA%\Mozilla\Firefox\Profiles\*\cookies.sqlite'))
print(f'\nFirefox cookie DBs: {ff_profiles}')

ff_cookies = {}
for ff_db in ff_profiles:
    tmp = os.path.join(tempfile.gettempdir(), 'ff_auth.db')
    try:
        shutil.copy2(ff_db, tmp)
        conn = sqlite3.connect(tmp)
        c = conn.cursor()
        c.execute("""SELECT name, value FROM moz_cookies
                     WHERE (host LIKE '%.x.com%' OR host='x.com'
                            OR host LIKE '%.twitter.com%' OR host='twitter.com')""")
        for name, value in c.fetchall():
            ff_cookies[name] = value
        print(f'  Found {len(ff_cookies)} FF cookies total')
        # Show auth-relevant ones
        for key in ('auth_token','ct0','twid','kdt','guest_id'):
            if key in ff_cookies:
                print(f'  {key} = {ff_cookies[key][:40]}...')
        conn.close()
    except Exception as e:
        print(f'FF error: {e}')
    finally:
        try: os.unlink(tmp)
        except: pass

print(f'\nFirefox X.com cookies collected: {list(ff_cookies.keys())}')

# ── 3. Build cookie string and call the GraphQL GrokShare API ────────────────
auth  = ff_cookies.get('auth_token', '')
ct0_v = ff_cookies.get('ct0', '')

if not auth:
    print('\n❌ No auth_token found in Firefox — session may be expired')
    sys.exit(1)

# X.com public bearer token (same one used in the browser JS bundle)
BEARER = 'AAAAAAAAAAAAAAAAAAAAANRILgAAAAAAnNwIzUejRCOuH5E6I%2BxMen7Lm0E%3D'

import urllib.parse, json as _json

variables = _json.dumps({"grok_share_id": SHARE_ID})
features  = _json.dumps({
    "premium_content_api_read_enabled": False,
    "communities_web_enable_tweet_community_results_fetch": True,
})

url = ('https://x.com/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare'
       '?variables=' + urllib.parse.quote(variables)
       + '&features=' + urllib.parse.quote(features))

cookie_str = '; '.join(f'{k}={v}' for k, v in ff_cookies.items())

print(f'\nCalling: {url[:100]}...')
req = urllib.request.Request(url)
req.add_header('Authorization', f'Bearer {BEARER}')
req.add_header('Cookie', cookie_str)
req.add_header('x-csrf-token', ct0_v)
req.add_header('User-Agent',
    'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 '
    '(KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36')
req.add_header('Accept', 'application/json')
req.add_header('Referer', f'https://x.com/i/grok/share/{SHARE_ID}')
req.add_header('x-twitter-auth-type', 'OAuth2Session')
req.add_header('x-twitter-client-language', 'en')
req.add_header('x-twitter-active-user', 'yes')

try:
    with urllib.request.urlopen(req, timeout=20) as resp:
        raw = resp.read().decode('utf-8')
        print(f'✅ HTTP {resp.status}: {len(raw):,} chars')
        raw_file = f'grok_share_{SHARE_ID}_raw.json'
        with open(raw_file, 'w', encoding='utf-8') as f:
            f.write(raw)
        print(f'Saved raw JSON → {raw_file}')
        print('\nPreview:', raw[:600])
except urllib.error.HTTPError as e:
    body = e.read().decode('utf-8', errors='replace')
    print(f'HTTP {e.code}: {body[:400]}')
    sys.exit(1)
except Exception as e:
    print(f'Error: {e}')
    sys.exit(1)
