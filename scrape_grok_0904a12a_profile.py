#!/usr/bin/env python3
"""
Scrape Grok thread 0904a12a using the user's existing Edge profile (with X login).
"""
import os, sys, time
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

url = 'https://x.com/i/grok/share/0904a12a5c2b4a639389ae084391b94f'
edge_profile = os.path.expandvars(r'%LOCALAPPDATA%\Microsoft\Edge\User Data')

edge_options = Options()
# Use existing Edge profile so X.com session cookies are available
edge_options.add_argument(f'--user-data-dir={edge_profile}')
edge_options.add_argument('--profile-directory=Default')
edge_options.add_argument('--disable-gpu')
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')
edge_options.add_argument('--disable-blink-features=AutomationControlled')
edge_options.add_experimental_option('excludeSwitches', ['enable-automation'])
edge_options.add_experimental_option('useAutomationExtension', False)

service = Service(executable_path='edge_driver\\msedgedriver.exe')
driver = webdriver.Edge(service=service, options=edge_options)

try:
    print(f'Loading: {url}')
    driver.get(url)
    time.sleep(15)  # Wait for Grok JS content to fully render

    body_text = driver.find_element(By.TAG_NAME, 'body').text

    if len(body_text) < 200:
        print('Content too short - may still be loading, waiting more...')
        time.sleep(10)
        body_text = driver.find_element(By.TAG_NAME, 'body').text

    print(f'Total content length: {len(body_text)} chars')
    print('=== FIRST 5000 CHARS ===')
    print(body_text[:5000])

    with open('grok_share_0904a12a5c2b4a639389ae084391b94f_content.txt', 'w', encoding='utf-8') as f:
        f.write(body_text)
    print('\nSaved to: grok_share_0904a12a5c2b4a639389ae084391b94f_content.txt')

    html = driver.page_source
    with open('grok_share_0904a12a5c2b4a639389ae084391b94f_source.html', 'w', encoding='utf-8') as f:
        f.write(html)
    print('HTML saved to: grok_share_0904a12a5c2b4a639389ae084391b94f_source.html')

finally:
    driver.quit()
    print('Browser closed.')
