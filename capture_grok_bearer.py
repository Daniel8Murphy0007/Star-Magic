#!/usr/bin/env python3
"""
Capture GrokShare API request headers from Edge browser via CDP.
Injects Firefox X.com cookies, loads share page, intercepts the outgoing
GrokShare API call to extract the exact Bearer token in use.
"""
import os, sys, time, json, sqlite3, shutil, tempfile
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options

SHARE_ID    = '7b78ffcb915f48bb90d55034c9c50b18'
TARGET_URL  = f'https://x.com/i/grok/share/{SHARE_ID}'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
FF_PROFILE  = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
FF_COOKIES  = os.path.join(FF_PROFILE, 'cookies.sqlite')

# ── Read Firefox X.com cookies ───────────────────────────────────────────────
print('Reading Firefox cookies...')
tmp = os.path.join(tempfile.gettempdir(), 'ff_ck_cdp.db')
shutil.copy2(FF_COOKIES, tmp)
cookies = []
try:
    conn = sqlite3.connect(tmp)
    cur  = conn.cursor()
    cur.execute("""
        SELECT host, name, value, path, expiry, isSecure, isHttpOnly
        FROM moz_cookies
        WHERE host LIKE '%.x.com%' OR host = 'x.com'
    """)
    cookies = cur.fetchall()
    conn.close()
finally:
    try: os.unlink(tmp)
    except: pass
print(f'Found {len(cookies)} X.com cookies')

# ── Launch Edge ──────────────────────────────────────────────────────────────
opts = Options()
opts.add_argument('--no-sandbox')
opts.add_argument('--disable-dev-shm-usage')
opts.add_argument('--disable-gpu')

svc = Service(DRIVER_PATH)
print('Launching Edge...')
driver = webdriver.Edge(service=svc, options=opts)

captured_bearer = None
captured_response = None

try:
    # Seed cookies
    driver.get('https://x.com/robots.txt')
    time.sleep(2)
    driver.delete_all_cookies()
    for host, name, value, path, expiry, secure, http_only in cookies:
        try:
            driver.add_cookie({
                'domain': host if not host.startswith('.') else host,
                'name':   name,
                'value':  value,
                'path':   path or '/',
                'secure': bool(secure),
                'httpOnly': bool(http_only),
            })
        except Exception:
            pass

    # Navigate to share page
    print(f'Loading {TARGET_URL}')
    driver.get(TARGET_URL)
    
    print('Waiting for page to initialise...')
    time.sleep(15)
    print('Page title:', driver.title)

    # Verify cookies are active
    cookies_in_browser = driver.execute_script("return document.cookie")
    has_auth = 'auth_token' in (cookies_in_browser or '')
    has_ct0  = 'ct0' in (cookies_in_browser or '')
    print(f'auth_token in browser: {has_auth}, ct0 in browser: {has_ct0}')

    # Make async fetch and capture status + body
    print('Executing async fetch...')
    import json as _json
    variables = _json.dumps({"grok_share_id": SHARE_ID})
    features  = _json.dumps({
        "premium_content_api_read_enabled": False,
        "communities_web_enable_tweet_community_results_fetch": True,
    })

    driver.set_script_timeout(30)
    js_async = f"""
    const callback = arguments[arguments.length - 1];
    const ct0 = (document.cookie.match(/ct0=([^;]+)/) || [])[1] || '';
    fetch(
        '/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare' +
        '?variables=' + encodeURIComponent({repr(variables)}) +
        '&features=' + encodeURIComponent({repr(features)}),
        {{
            method: 'GET',
            credentials: 'include',
            headers: {{
                'Content-Type': 'application/json',
                'x-twitter-active-user': 'yes',
                'x-twitter-client-language': 'en',
                'x-csrf-token': ct0
            }}
        }}
    ).then(async r => {{
        const body = await r.text();
        callback(JSON.stringify({{status: r.status, ok: r.ok, body: body}}));
    }}).catch(e => callback(JSON.stringify({{error: e.toString()}})));
    """

    result = driver.execute_async_script(js_async)
    print(f'Raw result: {result[:600] if result else "(empty)"}')

    if result:
        try:
            data = _json.loads(result)
            print(f'Status: {data.get("status")}, OK: {data.get("ok")}')
            body = data.get('body', '')
            if body:
                print(f'Body ({len(body)} bytes): {body[:800]}')
                captured_response = body
            else:
                print(f'Error: {data.get("error")}')
        except Exception:
            captured_response = result

    if not captured_response:
        body_text = driver.find_element('tag name', 'body').text
        page_src  = driver.page_source
        # innerText via JS (catches more content than Selenium .text)
        inner_text = driver.execute_script("return document.body.innerText") or ''
        
        print(f'\n.text: {len(body_text)} chars')
        print(f'page_source: {len(page_src)} chars')
        print(f'innerText: {len(inner_text)} chars')

        best = max([body_text, inner_text, page_src], key=len)
        best_name = ['.text', 'innerText', 'page_source'][
            [body_text, inner_text, page_src].index(best)]
        print(f'Largest source: {best_name} ({len(best)} chars)')

        if len(best) > 500:
            out = f'grok_share_{SHARE_ID}_content.txt'
            with open(out, 'w', encoding='utf-8') as f:
                f.write(best)
            print(f'Saved → {out}')
            print('\n=== PREVIEW ===')
            print(best[:3000])
        else:
            print('Content too short — page did not render')
            print('.text preview:', body_text)

finally:
    driver.quit()

# ── Save results ─────────────────────────────────────────────────────────────
if captured_bearer:
    with open('grok_bearer_live.txt', 'w') as f:
        f.write(captured_bearer)
    print(f'\nFull Bearer saved → grok_bearer_live.txt')

if captured_response:
    raw_file = f'grok_share_{SHARE_ID}_raw.json'
    with open(raw_file, 'w', encoding='utf-8') as f:
        f.write(captured_response)
    print(f'Raw API response saved → {raw_file}')
    try:
        data = json.loads(captured_response)
        print('JSON keys:', list(data.keys())[:10])
    except Exception:
        pass
