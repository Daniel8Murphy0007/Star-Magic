#!/usr/bin/env python3
"""
Use CDP network logging to capture what API calls X.com makes for Grok share page.
"""
import os, sys, time, sqlite3, shutil, tempfile, json

TARGET_URL = 'https://x.com/i/grok/share/ab2d0965e3a74a0da32749a7a2591dc7'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
FF_COOKIES_DB = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\cookies.sqlite')

def get_ff_cookies():
    tmp_db = os.path.join(tempfile.gettempdir(), 'ff_cdp.db')
    shutil.copy2(FF_COOKIES_DB, tmp_db)
    conn = sqlite3.connect(tmp_db)
    c = conn.cursor()
    c.execute("SELECT host, name, value, path, expiry, isSecure FROM moz_cookies WHERE host LIKE '%.x.com%' OR host='x.com'")
    rows = c.fetchall()
    conn.close()
    os.unlink(tmp_db)
    return rows

from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By

edge_options = Options()
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')
edge_options.set_capability('ms:loggingPrefs', {'performance': 'ALL'})

service = Service(executable_path=DRIVER_PATH)
driver = webdriver.Edge(service=service, options=edge_options)

try:
    # Enable CDP network events
    driver.execute_cdp_cmd('Network.enable', {})

    # Navigate to x.com and inject cookies
    print('Setting up auth...')
    driver.get('https://x.com/robots.txt')
    time.sleep(2)

    ff_cookies = get_ff_cookies()
    for (host, name, value, path, expiry, is_secure) in ff_cookies:
        cookie = {'name': name, 'value': value, 'domain': host, 'path': path or '/'}
        if expiry: cookie['expiry'] = expiry
        try:
            driver.add_cookie(cookie)
        except: pass

    print(f'Injected {len(ff_cookies)} cookies')
    print(f'Loading: {TARGET_URL}')
    driver.get(TARGET_URL)

    print('Waiting 30s while monitoring network...')
    time.sleep(30)

    # Capture performance logs
    logs = driver.get_log('performance')
    print(f'\nCaptured {len(logs)} performance log entries')

    # Filter for XHR/Fetch requests to API endpoints
    api_calls = []
    for log in logs:
        msg = json.loads(log['message'])
        method = msg.get('message', {}).get('method', '')
        params = msg.get('message', {}).get('params', {})

        if method in ('Network.requestWillBeSent', 'Network.responseReceived'):
            url = ''
            if method == 'Network.requestWillBeSent':
                url = params.get('request', {}).get('url', '')
            else:
                url = params.get('response', {}).get('url', '')

            if 'grok' in url.lower() or 'graphql' in url.lower() or 'ai' in url.lower():
                status = params.get('response', {}).get('status', '')
                print(f'[{method.split(".")[-1]}] {status} {url[:120]}')
                api_calls.append({'method': method, 'url': url, 'status': status})

    print(f'\n=== Summary: {len(api_calls)} Grok/GraphQL API calls ===')
    for call in api_calls[:20]:
        print(f'  {call["status"]} {call["url"][:100]}')

    # Save all API calls
    with open('grok_network_log.json', 'w') as f:
        json.dump(api_calls, f, indent=2)
    print('\nFull log saved: grok_network_log.json')

    # Also get current body text
    body = driver.find_element(By.TAG_NAME, 'body').text
    print(f'\nBody text ({len(body)} chars): {body[:300]}')

finally:
    driver.quit()
    print('Done.')
