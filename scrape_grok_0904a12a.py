#!/usr/bin/env python3
"""One-shot scraper for Grok share thread 0904a12a5c2b4a639389ae084391b94f"""
import time
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By

edge_options = Options()
edge_options.add_argument('--headless')
edge_options.add_argument('--disable-gpu')
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')
edge_options.add_argument('--disable-blink-features=AutomationControlled')
edge_options.add_experimental_option('excludeSwitches', ['enable-automation'])
edge_options.add_experimental_option('useAutomationExtension', False)

service = Service(executable_path='edge_driver\\msedgedriver.exe')
driver = webdriver.Edge(service=service, options=edge_options)
driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")

url = 'https://x.com/i/grok/share/0904a12a5c2b4a639389ae084391b94f'
print('Loading:', url)
driver.get(url)
time.sleep(12)  # Wait for full JS render

body_text = driver.find_element(By.TAG_NAME, 'body').text
print('=== PAGE TEXT (first 4000 chars) ===')
print(body_text[:4000])
print(f'\n=== TOTAL LENGTH: {len(body_text)} chars ===')

with open('grok_thread_0904a12a_raw.txt', 'w', encoding='utf-8') as f:
    f.write(body_text)
print('Full content saved to: grok_thread_0904a12a_raw.txt')
driver.quit()
