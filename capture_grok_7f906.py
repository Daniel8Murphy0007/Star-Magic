#!/usr/bin/env python3
"""
Capture Grok share thread: 7f90683d979346758875e4a8e6e0ee9d
Uses Firefox X.com cookies injected into Edge via Selenium CDP.
"""
import os, sys, time, json, sqlite3, shutil, tempfile
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options

SHARE_ID    = '7f90683d979346758875e4a8e6e0ee9d'
TARGET_URL  = f'https://x.com/i/grok/share/{SHARE_ID}'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
FF_PROFILE  = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles')
OUT_TXT     = f'grok_share_{SHARE_ID}_content.txt'

# ── Find Firefox cookies.sqlite ─────────────────────────────────────────────
import glob
ff_cookie_candidates = glob.glob(os.path.join(FF_PROFILE, '*', 'cookies.sqlite'))
FF_COOKIES = ff_cookie_candidates[0] if ff_cookie_candidates else None
print(f'Firefox cookies: {FF_COOKIES}')

cookies = []
if FF_COOKIES and os.path.exists(FF_COOKIES):
    tmp = os.path.join(tempfile.gettempdir(), 'ff_ck_capture.db')
    shutil.copy2(FF_COOKIES, tmp)
    try:
        conn = sqlite3.connect(tmp)
        cur  = conn.cursor()
        cur.execute("""
            SELECT host, name, value, path, expiry, isSecure, isHttpOnly
            FROM moz_cookies
            WHERE host LIKE '%.x.com%' OR host = 'x.com'
               OR host LIKE '%.twitter.com%' OR host = 'twitter.com'
        """)
        cookies = cur.fetchall()
        conn.close()
    finally:
        try: os.unlink(tmp)
        except: pass
print(f'Found {len(cookies)} X.com/Twitter cookies from Firefox')

# ── Launch Edge (headless=False so we can see what happens) ─────────────────
opts = Options()
opts.add_argument('--no-sandbox')
opts.add_argument('--disable-dev-shm-usage')
opts.add_argument('--disable-gpu')
opts.add_argument('--disable-extensions')
opts.add_argument('--mute-audio')
# Uncomment for fully headless: opts.add_argument('--headless=new')

svc = Service(DRIVER_PATH)
print('Launching Edge...')
driver = webdriver.Edge(service=svc, options=opts)

captured_text = None

try:
    # Seed cookies
    driver.get('https://x.com/robots.txt')
    time.sleep(2)
    driver.delete_all_cookies()
    for host, name, value, path, expiry, secure, http_only in cookies:
        try:
            driver.add_cookie({
                'domain':   host if not host.startswith('.') else host,
                'name':     name,
                'value':    value,
                'path':     path or '/',
                'secure':   bool(secure),
                'httpOnly': bool(http_only),
            })
        except Exception:
            pass
    print(f'Injected {len(cookies)} cookies')

    # Navigate to Grok share page
    print(f'Loading {TARGET_URL}')
    driver.get(TARGET_URL)

    print('Waiting 20s for JS rendering...')
    time.sleep(20)
    print('Page title:', driver.title)

    # Try intercepting via fetch API
    _json = __import__('json')
    grok_api_url = f'https://x.com/i/api/grok/1.1/grok_share_conversation?share_id={SHARE_ID}'
    script = f"""
        const response = await fetch('{grok_api_url}', {{
            headers: {{
                'Authorization': 'Bearer [PROMPT_FOR_GROK_WEB_BEARER_TOKEN]',
                'x-csrf-token': document.cookie.match(/ct0=([^;]+)/)?.[1] || '',
                'Content-Type': 'application/json',
            }},
            credentials: 'include'
        }});
        const status = response.status;
        const text = await response.text();
        return JSON.stringify({{ok: response.ok, status, body: text}});
    """
    try:
        result = driver.execute_async_script("""
            const callback = arguments[arguments.length - 1];
            (async () => {
                const shareId = '""" + SHARE_ID + """';
                const grokUrl = 'https://x.com/i/api/grok/1.1/grok_share_conversation?share_id=' + shareId;
                try {
                    const r = await fetch(grokUrl, {credentials: 'include'});
                    const t = await r.text();
                    callback(JSON.stringify({ok: r.ok, status: r.status, body: t}));
                } catch(e) {
                    callback(JSON.stringify({ok: false, error: e.toString()}));
                }
            })();
        """)
        print(f'API call result: {result[:300] if result else "(none)"}')
        if result:
            data = _json.loads(result)
            if data.get('ok') and data.get('body'):
                captured_text = data['body']
                print(f'Got API response: {len(captured_text)} chars')
    except Exception as e:
        print(f'API fetch failed: {e}')

    # Fall back: grab rendered page text
    if not captured_text or len(captured_text) < 500:
        print('Trying rendered page innerText...')
        inner_text = driver.execute_script("return document.body.innerText") or ''
        body_text  = driver.find_element('tag name', 'body').text or ''
        page_src   = driver.page_source or ''
        best = max([inner_text, body_text, page_src], key=len)
        best_name = ['innerText', 'body.text', 'page_source'][
            [inner_text, body_text, page_src].index(best)]
        print(f'Best source: {best_name} — {len(best)} chars')
        if len(best) > 1000:
            captured_text = best
        else:
            print('WARNING: page did not render enough content')
            print('Preview:', best[:500])

    # ── Scroll and wait more if content is thin ──────────────────────────────
    if captured_text and len(captured_text) < 5000:
        print('Content thin, scrolling and waiting 15s more...')
        driver.execute_script("window.scrollTo(0, document.body.scrollHeight)")
        time.sleep(15)
        inner_text2 = driver.execute_script("return document.body.innerText") or ''
        if len(inner_text2) > len(captured_text):
            captured_text = inner_text2
            print(f'After scroll: {len(captured_text)} chars')

finally:
    driver.quit()

# ── Save ─────────────────────────────────────────────────────────────────────
if captured_text and len(captured_text) > 500:
    with open(OUT_TXT, 'w', encoding='utf-8') as f:
        f.write(captured_text)
    print(f'\n✅ Saved → {OUT_TXT} ({len(captured_text):,} chars)')
    print('\n=== PREVIEW (first 4000 chars) ===')
    print(captured_text[:4000])
else:
    print('\n❌ Could not capture meaningful content')
    print('Please manually paste the Grok thread content into:')
    print(f'  {OUT_TXT}')
