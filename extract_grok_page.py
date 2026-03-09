#!/usr/bin/env python3
"""
Extract Grok share content from rendered Edge page body.
Re-uses the Edge session approach — page renders fully even without auth_token
because the SHARE link is publicly accessible once rendered.
"""
import os, sys, time, sqlite3, shutil, tempfile
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

SHARE_ID    = '7b78ffcb915f48bb90d55034c9c50b18'
TARGET_URL  = f'https://x.com/i/grok/share/{SHARE_ID}'
OUTPUT_FILE = f'grok_share_{SHARE_ID}_content.txt'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
FF_PROFILE  = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
FF_COOKIES  = os.path.join(FF_PROFILE, 'cookies.sqlite')

# ── Read Firefox X.com cookies ───────────────────────────────────────────────
print('Reading Firefox cookies...')
tmp = os.path.join(tempfile.gettempdir(), 'ff_ck_extract.db')
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

try:
    # Inject cookies
    driver.get('https://x.com/robots.txt')
    time.sleep(2)
    driver.delete_all_cookies()
    for host, name, value, path, expiry, secure, http_only in cookies:
        try:
            driver.add_cookie({
                'domain':   host,
                'name':     name,
                'value':    value,
                'path':     path or '/',
                'secure':   bool(secure),
                'httpOnly': bool(http_only),
            })
        except Exception:
            pass

    # Load share page
    print(f'Loading {TARGET_URL}')
    driver.get(TARGET_URL)
    print('Waiting 15s for initial page load...')
    time.sleep(15)
    print('Page title:', driver.title)

    # Trigger React hydration by making a fetch call inside x.com's JS context
    # (same mechanism that caused capture_grok_bearer.py to get 957K chars)
    print('Triggering in-page fetch to activate React rendering...')
    import json as _json
    variables = _json.dumps({"grok_share_id": SHARE_ID})
    features  = _json.dumps({
        "premium_content_api_read_enabled": False,
        "communities_web_enable_tweet_community_results_fetch": True,
    })
    driver.set_script_timeout(30)
    driver.execute_async_script(f"""
        const callback = arguments[arguments.length - 1];
        const ct0 = (document.cookie.match(/ct0=([^;]+)/) || [])[1] || '';
        fetch('/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare'
            + '?variables=' + encodeURIComponent({repr(variables)})
            + '&features=' + encodeURIComponent({repr(features)}),
            {{ method:'GET', credentials:'include',
               headers:{{'x-csrf-token': ct0, 'x-twitter-active-user':'yes'}} }}
        ).then(r => r.text())
         .then(t => callback(t))
         .catch(e => callback('err'));
    """)

    # Now poll until full content appears (should be fast now)
    print('Polling for content after trigger...')
    for i in range(10):
        time.sleep(3)
        size = len(driver.find_element(By.TAG_NAME, 'body').text)
        print(f'  {(i+1)*3}s: {size} chars')
        if size > 10000:
            print('  Content loaded!')
            break

    body_text = driver.find_element(By.TAG_NAME, 'body').text
    print(f'Page body: {len(body_text)} chars')

    if len(body_text) > 200:
        with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
            f.write(body_text)
        print(f'Saved → {OUTPUT_FILE}')
        print('\n=== PREVIEW (first 3000 chars) ===')
        print(body_text[:3000])
    else:
        print(f'Content too short ({len(body_text)} chars) — page may not have loaded')
        print('Body:', body_text)

finally:
    driver.quit()
