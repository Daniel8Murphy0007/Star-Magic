#!/usr/bin/env python3
"""
Scrape Grok share thread ab2d0965e3a74a0da32749a7a2591dc7
Uses auth cookies from Firefox profile (X.com session) injected into Edge Selenium session.
"""
import os, sys, time, sqlite3, shutil, tempfile

TARGET_URL = 'https://x.com/i/grok/share/ab2d0965e3a74a0da32749a7a2591dc7'
OUTPUT_FILE = 'grok_share_ab2d0965e3a74a0da32749a7a2591dc7_content.txt'
HTML_FILE   = 'grok_share_ab2d0965e3a74a0da32749a7a2591dc7_source.html'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'

FF_PROFILE = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release')
FF_COOKIES_DB = os.path.join(FF_PROFILE, 'cookies.sqlite')

# --- Step 1: Read auth cookies from Firefox ---
print('Reading X.com auth cookies from Firefox...')
tmp_db = os.path.join(tempfile.gettempdir(), 'ff_ck_inject.db')
shutil.copy2(FF_COOKIES_DB, tmp_db)

try:
    conn = sqlite3.connect(tmp_db)
    c = conn.cursor()
    c.execute("""
        SELECT host, name, value, path, expiry, isSecure, isHttpOnly
        FROM moz_cookies
        WHERE host LIKE '%.x.com%' OR host = 'x.com'
        ORDER BY name
    """)
    ff_cookies = c.fetchall()
    conn.close()
finally:
    try: os.unlink(tmp_db)
    except: pass

if not ff_cookies:
    print('ERROR: No X.com cookies found in Firefox profile')
    sys.exit(1)

print(f'Found {len(ff_cookies)} X.com cookies in Firefox')
cookie_names = [r[1] for r in ff_cookies]
print(f'Cookie names: {cookie_names}')

# Check critical cookies
auth_token = next((r[2] for r in ff_cookies if r[1] == 'auth_token'), None)
ct0 = next((r[2] for r in ff_cookies if r[1] == 'ct0'), None)
if auth_token:
    print(f'auth_token: {auth_token[:20]}...')
if ct0:
    print(f'ct0: {ct0[:20]}...')

# --- Step 2: Launch Edge with fresh profile (no user-data-dir) ---
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By

edge_options = Options()
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')
edge_options.add_argument('--disable-blink-features=AutomationControlled')
edge_options.add_experimental_option('excludeSwitches', ['enable-automation'])
edge_options.add_experimental_option('useAutomationExtension', False)

service = Service(executable_path=DRIVER_PATH)

try:
    print('Starting Edge (fresh profile, no lock issues)...')
    driver = webdriver.Edge(service=service, options=edge_options)
except Exception as e:
    print(f'ERROR starting Edge: {e}')
    sys.exit(1)

try:
    # --- Step 3: Navigate to x.com to set cookie domain ---
    print('Navigating to x.com to set cookies...')
    driver.get('https://x.com/robots.txt')
    time.sleep(2)

    # --- Step 4: Inject Firefox cookies ---
    print('Injecting auth cookies...')
    for (host, name, value, path, expiry, is_secure, is_httponly) in ff_cookies:
        # Clean up domain (selenium wants domain without leading dot for some cookies)
        domain = host if not host.startswith('.') else host
        cookie = {
            'name': name,
            'value': value,
            'domain': domain,
            'path': path or '/',
            'secure': bool(is_secure),
        }
        if expiry:
            cookie['expiry'] = expiry
        try:
            driver.add_cookie(cookie)
        except Exception as ce:
            print(f'  Skipped cookie {name}: {ce}')
    print(f'Cookies injected.')

    # --- Step 5: Load the target URL ---
    print(f'Loading: {TARGET_URL}')
    driver.get(TARGET_URL)

    print('Waiting up to 90s for Grok conversation content...')
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC

    def grok_content_loaded(drv):
        txt = drv.find_element(By.TAG_NAME, 'body').text
        # Content is loaded when body has substantial text beyond just the nav
        return len(txt) > 500 and 'log in' not in txt.lower()

    try:
        WebDriverWait(driver, 90).until(grok_content_loaded)
        print('Content detected!')
    except:
        print('Timeout waiting for content — capturing what is available')

    # Scroll down to trigger any lazy loading
    driver.execute_script('window.scrollTo(0, document.body.scrollHeight)')
    time.sleep(5)
    driver.execute_script('window.scrollTo(0, 0)')
    time.sleep(3)

    body_text = driver.find_element(By.TAG_NAME, 'body').text

    print(f'Total content: {len(body_text)} chars')
    print('=== FIRST 3000 CHARS ===')
    print(body_text[:3000])

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(body_text)
    print(f'\nSaved to: {OUTPUT_FILE} ({len(body_text)} bytes)')

    html = driver.page_source
    with open(HTML_FILE, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'HTML saved to: {HTML_FILE} ({len(html)} bytes)')

finally:
    driver.quit()
    print('Browser closed.')
