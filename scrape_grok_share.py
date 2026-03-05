#!/usr/bin/env python3
"""
Scrape Grok share URL for physics and mathematical content.
Usage: python scrape_grok_share.py <URL>
"""

import os
import sys
import time
import json
import re
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

# Get URL from command line
if len(sys.argv) < 2:
    print("Usage: python scrape_grok_share.py <URL>")
    sys.exit(1)

url = sys.argv[1]

# Extract ID from URL for filename
url_id_match = re.search(r'/share/([a-f0-9]+)', url)
url_id = url_id_match.group(1) if url_id_match else 'unknown'

# Initialize Edge WebDriver
edge_driver_path = os.path.join(os.getcwd(), 'edge_driver', 'msedgedriver.exe')

edge_options = Options()
edge_options.add_argument('--disable-gpu')
edge_options.add_argument('--no-sandbox')
edge_options.add_argument('--disable-dev-shm-usage')
edge_options.add_argument('--window-size=1920,1080')

service = Service(executable_path=edge_driver_path)
driver = webdriver.Edge(service=service, options=edge_options)

try:
    print(f"Loading: {url}")
    driver.get(url)
    
    # Wait for page to load
    time.sleep(10)
    
    # Get page source
    page_source = driver.page_source
    
    # Try to find content elements
    try:
        # Wait for main content
        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.TAG_NAME, "body"))
        )
        
        # Extract all text content
        body = driver.find_element(By.TAG_NAME, "body")
        full_text = body.text
        
        print(f"\n{'='*80}")
        print("SCRAPED CONTENT:")
        print('='*80)
        print(full_text)
        print('='*80)
        
        # Save to file with unique ID
        content_file = f'grok_share_{url_id}_content.txt'
        with open(content_file, 'w', encoding='utf-8') as f:
            f.write(full_text)
        
        print(f"\n✓ Content saved to: {content_file}")
        print(f"✓ Content size: {len(full_text)} bytes")
        
        # Save HTML source
        html_file = f'grok_share_{url_id}_source.html'
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(page_source)
        
        print(f"✓ HTML source saved to: {html_file}")
        
    except Exception as e:
        print(f"Error extracting content: {e}")
        # Save what we have
        html_file = f'grok_share_{url_id}_source.html'
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(page_source)
        print(f"✓ HTML source saved (extraction failed): {html_file}")
    
finally:
    driver.quit()
    print("\n✓ Browser closed")
