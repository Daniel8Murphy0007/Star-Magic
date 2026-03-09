#!/usr/bin/env python3
"""
Verify auth_token is properly injected and sent in network requests.
Tests by calling a simple X.com API and checking if 401 vs 403/200 response.
"""
import os, sys, time, sqlite3, shutil, tempfile, json
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options

DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
FF_PROFILE  = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
FF_COOKIES  = os.path.join(FF_PROFILE, 'cookies.sqlite')

# Read Firefox X.com cookies
tmp = os.path.join(tempfile.gettempdir(), 'ff_ck_verify.db')
shutil.copy2(FF_COOKIES, tmp)
all_cookies = []
try:
    conn = sqlite3.connect(tmp)
    cur  = conn.cursor()
    cur.execute("SELECT host, name, value, path, expiry, isSecure, isHttpOnly FROM moz_cookies WHERE host LIKE '%.x.com%' OR host = 'x.com'")
    all_cookies = cur.fetchall()
    conn.close()
finally:
    try: os.unlink(tmp)
    except: pass

opts = Options()
opts.add_argument('--no-sandbox')
opts.add_argument('--disable-dev-shm-usage')
opts.add_argument('--disable-gpu')

svc = Service(DRIVER_PATH)
driver = webdriver.Edge(service=svc, options=opts)

try:
    driver.get('https://x.com/robots.txt')
    time.sleep(2)
    driver.delete_all_cookies()
    
    set_ok = 0
    for host, name, value, path, expiry, secure, http_only in all_cookies:
        try:
            driver.add_cookie({
                'domain':   host,        # '.x.com' with dot = all subdomains
                'name':     name,
                'value':    value,
                'path':     path or '/',
                'secure':   bool(secure),
                'httpOnly': bool(http_only),
            })
            set_ok += 1
        except Exception as e:
            print(f'Cookie fail: {name} ({host}): {e}')
    
    print(f'Cookies set: {set_ok}/{len(all_cookies)}')
    
    # Verify auth_token IS in browser (even if httpOnly, it shows in driver.get_cookies())
    browser_cookies = {c['name']: c for c in driver.get_cookies()}
    print(f'auth_token in driver.get_cookies(): {"auth_token" in browser_cookies}')
    print(f'ct0 in driver.get_cookies(): {"ct0" in browser_cookies}')
    if 'auth_token' in browser_cookies:
        c = browser_cookies['auth_token']
        print(f'auth_token httpOnly: {c.get("httpOnly")}, domain: {c.get("domain")}, secure: {c.get("secure")}')

    # Test auth: call a simple authenticated endpoint
    # GET /1.1/account/verify_credentials.json returns 200 if logged in, 401 if not
    driver.set_script_timeout(15)
    result = driver.execute_async_script("""
        const cb = arguments[arguments.length - 1];
        const ct0 = (document.cookie.match(/ct0=([^;]+)/) || [])[1] || '';
        fetch('https://x.com/i/api/1.1/account/settings.json', {
            credentials: 'include',
            headers: {'x-csrf-token': ct0, 'x-twitter-active-user': 'yes'}
        }).then(r => cb(String(r.status))).catch(e => cb('err:' + e));
    """)
    print(f'Account settings endpoint: HTTP {result}')
    
    if result == '200':
        print('AUTH WORKING - session is valid')
    elif result == '403':
        print('Auth sent but forbidden (rate-limited or needs fresh Bearer)')
    elif result == '401':
        print('NOT AUTHENTICATED - auth_token not being sent')
    else:
        print(f'Unexpected response: {result}')

    # Now navigate to Grok share and wait longer
    print('\nNavigating to Grok share...')
    driver.get('https://x.com/i/grok/share/7b78ffcb915f48bb90d55034c9c50b18')
    time.sleep(20)

    # Try multiple text extraction methods
    text_sel  = driver.find_element('tag name', 'body').text
    text_cont = driver.execute_script("return document.body.textContent") or ''
    inner_txt = driver.execute_script("return document.body.innerText") or ''
    
    print(f'.text:        {len(text_sel)} chars')
    print(f'textContent:  {len(text_cont)} chars')
    print(f'innerText:    {len(inner_txt)} chars')
    
    # Preview whichever has most content
    best = max(text_sel, text_cont, inner_txt, key=len)
    best_name = max([('.text', text_sel), ('textContent', text_cont), ('innerText', inner_txt)], key=lambda x: len(x[1]))[0]
    print(f'\nBest: {best_name} ({len(best)} chars)')
    if len(best) > 400:
        with open('grok_share_7b78ffcb915f48bb90d55034c9c50b18_content.txt', 'w', encoding='utf-8') as f:
            f.write(best)
        print(f'Saved!')
        print(best[:3000])
    else:
        print('Still no content. Preview of textContent:')
        print(text_cont[:500])

finally:
    driver.quit()
