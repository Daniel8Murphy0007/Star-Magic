#!/usr/bin/env python3
"""
Stealth Edge + CDP approach:
1. Disable webdriver detection BEFORE page load
2. Inject Firefox cookies
3. Navigate to Grok share - browser makes GrokShare API call automatically
4. Wait for full React hydration
5. Read grokShare Redux state from __INITIAL_STATE__ after client-side update
   OR intercept the network response via CDP
"""
import sqlite3, shutil, os, json, time, re

SHARE_ID   = '7b78ffcb915f48bb90d55034c9c50b18'
SHARE_URL  = f'https://x.com/i/grok/share/{SHARE_ID}'
FF_PROFILE = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
COOKIES_DB = os.path.join(FF_PROFILE, 'cookies.sqlite')
TMP_DB     = os.path.join(os.environ['TEMP'], 'ff_ck_stealth.sqlite')
EDGE_PATH  = r'C:\edgedriver_win64\msedgedriver.exe'
OUTPUT     = f'grok_stealth_{SHARE_ID}.txt'

from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By

# ── Read Firefox cookies ─────────────────────────────────────────────────────
shutil.copy2(COOKIES_DB, TMP_DB)
conn = sqlite3.connect(TMP_DB)
ff_cookies = conn.execute(
    "SELECT name, value, host, path, expiry, isSecure FROM moz_cookies WHERE host LIKE '%x.com'"
).fetchall()
conn.close()
print(f'Cookies: {len(ff_cookies)}')

ct0 = next((v for n,v,*_ in ff_cookies if n=='ct0'), None)
auth_t = next((v for n,v,*_ in ff_cookies if n=='auth_token'), None)
print(f'ct0={ct0[:10] if ct0 else "MISSING"}..., auth_token={"YES" if auth_t else "MISSING"}')

# ── Configure Edge with stealth ──────────────────────────────────────────────
options = Options()
options.add_argument('--disable-blink-features=AutomationControlled')
options.add_argument('--no-sandbox')
options.add_argument('--disable-dev-shm-usage')
options.add_experimental_option('excludeSwitches', ['enable-automation'])
options.add_experimental_option('useAutomationExtension', False)
options.add_argument('--user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/134.0.0.0 Safari/537.36 Edg/134.0.0.0')

service = Service(EDGE_PATH)
print('Launching Edge (stealth mode)...')
driver = webdriver.Edge(service=service, options=options)

# ── Patch navigator.webdriver BEFORE any navigation ──────────────────────────
driver.execute_cdp_cmd('Page.addScriptToEvaluateOnNewDocument', {
    'source': '''
        Object.defineProperty(navigator, 'webdriver', {get: () => undefined});
        Object.defineProperty(navigator, 'plugins', {get: () => [1,2,3,4,5]});
        window.chrome = {runtime: {}};
    '''
})

# ── Enable Network capture to sniff the GrokShare API response ───────────────
driver.execute_cdp_cmd('Network.enable', {})

# ── Navigate to x.com first to set cookies ───────────────────────────────────
print('Navigating to x.com...')
driver.get('https://x.com/robots.txt')
time.sleep(2)

# Clear any existing cookies and set Firefox ones
driver.delete_all_cookies()
for name, value, host, path, expiry, secure in ff_cookies:
    domain = host.lstrip('.')
    try:
        ck = {
            'name':   name,
            'value':  value,
            'domain': domain,
            'path':   path or '/',
            'secure': bool(secure),
        }
        if expiry: ck['expiry'] = int(expiry)
        driver.add_cookie(ck)
    except Exception as e:
        pass  # Some cookies may fail

print(f'Cookies set. Navigating to Grok share...')

# ── Navigate to the share URL ─────────────────────────────────────────────────
driver.get(SHARE_URL)
print('Page loading...')

# ── Wait up to 30s for GrokShare content to appear ───────────────────────────
# We look for: content longer than 200 chars appearing in the main content area
result_text = None
for i in range(30):
    time.sleep(1)
    try:
        body_text = driver.execute_script('return document.body.innerText')
        if body_text and len(body_text) > 500:
            result_text = body_text
            print(f'Got content: {len(body_text)} chars (after {i+1}s)')
            break
        else:
            if i % 5 == 0:
                print(f'  [{i+1}s] innerText: {len(body_text) if body_text else 0} chars')
    except Exception as e:
        print(f'  [{i+1}s] error: {e}')

# ── Try to get expanded content via fetch from within page ───────────────────
print('\n--- Attempting in-page fetch ---')
try:
    fetch_result = driver.execute_async_script('''
        var done = arguments[0];
        var url = arguments[1];
        fetch(url, {credentials: 'include'})
            .then(r => r.text())
            .then(t => done({status: 'ok', text: t, len: t.length}))
            .catch(e => done({status: 'err', msg: String(e)}));
    ''', SHARE_URL)
    print(f'Fetch result: {fetch_result.get("status")} len={fetch_result.get("len", 0)}')
    if fetch_result.get('len', 0) > 10000:
        page_html = fetch_result['text']
        with open(f'grok_stealth_fetch_{SHARE_ID}.html', 'w', encoding='utf-8') as f:
            f.write(page_html)
        print(f'Saved fetch HTML ({len(page_html)} chars)')
except Exception as e:
    print(f'Fetch error: {e}')

# ── Try to read Redux state (injected after hydration) ───────────────────────
print('\n--- Reading Redux/React state ---')
try:
    grok_data = driver.execute_script('''
        // Try multiple ways to access the store
        var state = null;
        
        // Method 1: window.__INITIAL_STATE__
        if (window.__INITIAL_STATE__ && window.__INITIAL_STATE__.entities) {
            var gs = window.__INITIAL_STATE__.entities.grokShare;
            if (gs && gs.entities && Object.keys(gs.entities).length > 0) {
                return JSON.stringify({method: "initial_state", data: gs});
            }
        }
        
        // Method 2: React root internal state (React 18 fiber)
        try {
            var root = document.getElementById("react-root");
            if (root) {
                var fiber = root._reactFiber || Object.keys(root).find(k => k.startsWith("__reactFiber"));
                if (fiber) return JSON.stringify({method: "fiber", found: true});
            }
        } catch(e) {}
        
        // Method 3: Look for rendered conversation text
        var all = Array.from(document.querySelectorAll("article, [data-testid], [class*='grok'], [class*='conv']"))
            .map(el => el.textContent)
            .filter(t => t.length > 100);
        if (all.length > 0) return JSON.stringify({method: "rendered", items: all.length, sample: all[0].substring(0,200)});
        
        return JSON.stringify({method: "none", bodyLen: document.body.innerText.length});
    ''')
    print(f'Redux state check: {grok_data}')
except Exception as e:
    print(f'Redux check error: {e}')

# ── Get full page source for analysis ─────────────────────────────────────────
page_src = driver.page_source
print(f'\nPage source: {len(page_src)} chars')

body_txt = driver.execute_script('return document.body.innerText') or ''
print(f'Body innerText: {len(body_txt)} chars')
if body_txt:
    print('Preview:')
    print(body_txt[:500])

# Save all
with open(OUTPUT, 'w', encoding='utf-8') as f:
    f.write(f'URL: {SHARE_URL}\n')
    f.write(f'Body innerText ({len(body_txt)} chars):\n')
    f.write(body_txt)
    if result_text and result_text != body_txt:
        f.write(f'\n\nEarly capture ({len(result_text)} chars):\n')
        f.write(result_text)

print(f'\nSaved -> {OUTPUT}')
driver.quit()
print('Done.')
