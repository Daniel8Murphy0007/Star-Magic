#!/usr/bin/env python3
"""
Direct X.com API approach: use Firefox auth_token to call Grok conversation API.
Bypasses Selenium rendering issues entirely.
"""
import os, sqlite3, shutil, tempfile, json, sys

# Read full cookies from Firefox
FF_COOKIES_DB = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\cookies.sqlite')
tmp_db = os.path.join(tempfile.gettempdir(), 'ff_api_inject.db')
shutil.copy2(FF_COOKIES_DB, tmp_db)

conn = sqlite3.connect(tmp_db)
c = conn.cursor()
c.execute("SELECT host, name, value FROM moz_cookies WHERE host LIKE '%.x.com%' OR host = 'x.com'")
ff_cookies = c.fetchall()
conn.close()
os.unlink(tmp_db)

cookies_dict = {name: value for host, name, value in ff_cookies}
auth_token = cookies_dict.get('auth_token', '')
ct0 = cookies_dict.get('ct0', '')
print(f'auth_token: {auth_token[:20]}...')
print(f'ct0: {ct0[:20]}...')

import urllib.request, urllib.error

CONV_ID = 'ab2d0965e3a74a0da32749a7a2591dc7'

# Build cookie header from all x.com cookies
cookie_header = '; '.join(f'{name}={value}' for host, name, value in ff_cookies if 'x.com' in host and not name.startswith('_ga') and name not in ('__cf_bm',))

headers = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36',
    'Accept': 'application/json',
    'Accept-Language': 'en-US,en;q=0.9',
    'Authorization': 'Bearer AAAAAAAAAAAAAAAAAAAAANRILgAAAAAAnNwIzUejRCOuH5E6I38u8zUk38Y%3DUouWPwtxNUezdzobearsAHsJWnk5GbOf91kX0yn2Ez44MUA',
    'Cookie': cookie_header,
    'x-csrf-token': ct0,
    'x-twitter-auth-type': 'OAuth2Session',
    'x-twitter-client-language': 'en',
    'Referer': f'https://x.com/i/grok/share/{CONV_ID}',
}

# Try multiple API endpoint patterns
endpoints = [
    f'https://api.x.com/2/grok/shared_conversation/{CONV_ID}',
    f'https://x.com/i/api/2/grok/shared_conversation/{CONV_ID}',
    f'https://x.com/i/api/graphql/shared_conversation?convId={CONV_ID}',
]

for url in endpoints:
    print(f'\nTrying: {url}')
    try:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=15) as resp:
            body = resp.read().decode('utf-8', errors='replace')
            print(f'Status: {resp.status}  Size: {len(body)} bytes')
            print('First 1000:', body[:1000])
            if len(body) > 100:
                with open(f'grok_api_response_{CONV_ID}.json', 'w') as f:
                    f.write(body)
                print(f'Saved: grok_api_response_{CONV_ID}.json')
                break
    except urllib.error.HTTPError as e:
        body = e.read().decode('utf-8', errors='replace')
        print(f'HTTP {e.code}: {body[:300]}')
    except Exception as ex:
        print(f'Error: {ex}')
