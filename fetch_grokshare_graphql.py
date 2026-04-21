#!/usr/bin/env python3
"""
Call X.com's GrokShare GraphQL endpoint directly using Firefox auth cookies.
Endpoint found via CDP monitoring: /i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare
"""
import os, sys, sqlite3, shutil, tempfile, json
import requests

CONV_ID = 'ab2d0965e3a74a0da32749a7a2591dc7'
FF_DB = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\cookies.sqlite')

# Read Firefox cookies
tmp = os.path.join(tempfile.gettempdir(), 'ff_gs.db')
shutil.copy2(FF_DB, tmp)
conn = sqlite3.connect(tmp)
c = conn.cursor()
c.execute("SELECT host, name, value FROM moz_cookies WHERE host LIKE '%.x.com%' OR host='x.com'")
rows = c.fetchall()
conn.close()
os.unlink(tmp)

cookies = {name: value for host, name, value in rows}
ct0 = cookies.get('ct0', '')
print(f'Cookies loaded: {len(cookies)}, ct0: {ct0[:20]}...')

# GrokShare GraphQL endpoint (discovered via CDP)
url = 'https://x.com/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare'
params = {'variables': json.dumps({'grok_share_id': CONV_ID})}

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/145.0.0.0 Safari/537.36 Edg/145.0.0.0',
    'Accept': '*/*',
    'Accept-Language': 'en-US,en;q=0.9',
    'Authorization': 'Bearer [PROMPT_FOR_GROK_WEB_BEARER_TOKEN]',
    'x-csrf-token': ct0,
    'x-twitter-auth-type': 'OAuth2Session',
    'x-twitter-client-language': 'en',
    'x-twitter-client-version': 'Twitter-TweetDeck-blackbird-chrome/V900.2026.1',
    'Referer': f'https://x.com/i/grok/share/{CONV_ID}',
    'Origin': 'https://x.com',
}

session = requests.Session()
for host, name, value in rows:
    session.cookies.set(name, value, domain='.x.com')

print(f'Calling GrokShare GraphQL...')
resp = session.get(url, params=params, headers=headers, timeout=30)
print(f'Status: {resp.status_code}')
print(f'Size: {len(resp.text)} bytes')
print(f'First 2000 chars: {resp.text[:2000]}')

if resp.status_code == 200 and len(resp.text) > 100:
    out = f'grok_share_{CONV_ID}_graphql.json'
    with open(out, 'w', encoding='utf-8') as f:
        f.write(resp.text)
    print(f'\nSaved: {out}')

    # Parse and extract conversation text
    try:
        data = resp.json()
        print('\n=== Parsed structure (keys) ===')
        def show_keys(obj, depth=0, prefix=''):
            if isinstance(obj, dict):
                for k, v in obj.items():
                    print('  ' * depth + str(k))
                    if depth < 3:
                        show_keys(v, depth+1)
            elif isinstance(obj, list) and obj:
                print('  ' * depth + f'[list len={len(obj)}]')
                show_keys(obj[0], depth+1)
        show_keys(data)
    except Exception as e:
        print(f'Parse error: {e}')
