#!/usr/bin/env python3
"""
Authenticated Grok share scraper for thread 0a48404556d54ab2a0c9e5b7cafd81be
Uses existing Edge profile (with X.com session) to access the share URL.

REQUIREMENT: Close Microsoft Edge COMPLETELY before running this script.
             If Edge is running, the profile is locked and this will fail.

Usage: python scrape_grok_0a484045_profile.py
"""
import os, sys, time
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.microsoft import EdgeChromiumDriverManager

TARGET_URL = 'https://x.com/i/grok/share/ab2d0965e3a74a0da32749a7a2591dc7'
OUTPUT_FILE = 'grok_share_ab2d0965e3a74a0da32749a7a2591dc7_content.txt'
HTML_FILE   = 'grok_share_ab2d0965e3a74a0da32749a7a2591dc7_source.html'
# Use webdriver-manager to auto-match driver to installed Edge version
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'  # fallback

edge_profile = os.path.expandvars(r'%LOCALAPPDATA%\Microsoft\Edge\User Data')
if not os.path.isdir(edge_profile):
    print(f'ERROR: Edge profile not found at {edge_profile}')
    sys.exit(1)

print(f'Using Edge profile: {edge_profile}')
print(f'Driver: {DRIVER_PATH}')

edge_options = Options()
edge_options.add_argument(f'--user-data-dir={edge_profile}')
edge_options.add_argument('--profile-directory=Default')
edge_options.add_argument('--disable-gpu')
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')
edge_options.add_argument('--disable-blink-features=AutomationControlled')
edge_options.add_experimental_option('excludeSwitches', ['enable-automation'])
edge_options.add_experimental_option('useAutomationExtension', False)

service = Service(executable_path=DRIVER_PATH)

try:
    # Try auto-matching driver first
    try:
        auto_driver_path = EdgeChromiumDriverManager().install()
        print(f'Auto-driver: {auto_driver_path}')
        service = Service(executable_path=auto_driver_path)
    except Exception as dm_err:
        print(f'webdriver-manager failed ({dm_err}), using fallback: {DRIVER_PATH}')
        service = Service(executable_path=DRIVER_PATH)
    driver = webdriver.Edge(service=service, options=edge_options)
except Exception as e:
    print(f'ERROR starting Edge: {e}')
    print('Make sure Edge is completely closed and try again.')
    sys.exit(1)

try:
    print(f'Loading: {TARGET_URL}')
    driver.get(TARGET_URL)

    # Wait for Grok's content to fully render (JS-heavy page)
    print('Waiting for page render (15s)...')
    time.sleep(15)

    body_text = driver.find_element(By.TAG_NAME, 'body').text
    if len(body_text) < 500:
        print(f'Content short ({len(body_text)} chars), waiting additional 15s...')
        time.sleep(15)
        body_text = driver.find_element(By.TAG_NAME, 'body').text

    print(f'Total content: {len(body_text)} chars')
    print('=== FIRST 5000 CHARS ===')
    print(body_text[:5000])

    with open(OUTPUT_FILE, 'w', encoding='utf-8') as f:
        f.write(body_text)
    print(f'\n✓ Saved to: {OUTPUT_FILE}')
    print(f'✓ Content size: {len(body_text)} bytes')

    html = driver.page_source
    with open(HTML_FILE, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'✓ HTML saved to: {HTML_FILE}')

finally:
    driver.quit()
    print('Browser closed.')
