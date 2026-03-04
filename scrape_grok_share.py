#!/usr/bin/env python3
"""
Scrape Grok share URL for physics and mathematical content.
"""

import os
import time
import json
from selenium import webdriver
from selenium.webdriver.edge.service import Service
from selenium.webdriver.edge.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

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
    url = "https://x.com/i/grok/share/4e0ecf23920b435cb3b2f410e93699b5"
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
        
        # Save to file
        with open('grok_share_4e0ecf23_content.txt', 'w', encoding='utf-8') as f:
            f.write(full_text)
        
        print("\n✓ Content saved to: grok_share_4e0ecf23_content.txt")
        
        # Save HTML source
        with open('grok_share_4e0ecf23_source.html', 'w', encoding='utf-8') as f:
            f.write(page_source)
        
        print("✓ HTML source saved to: grok_share_4e0ecf23_source.html")
        
    except Exception as e:
        print(f"Error extracting content: {e}")
        # Save what we have
        with open('grok_share_4e0ecf23_source.html', 'w', encoding='utf-8') as f:
            f.write(page_source)
        print("✓ HTML source saved (extraction failed)")
    
finally:
    driver.quit()
    print("\n✓ Browser closed")
