#!/usr/bin/env python3
"""
Scrape Grok Thread e3a0c453 - Star Magic UQFF Physics
URL: https://x.com/i/grok/share/e3a0c4534d01419f95020d8393cf0023
Date: March 8, 2026
"""

import requests
from bs4 import BeautifulSoup
import json
import time

def scrape_grok_thread(url, thread_id):
    """Scrape Grok shared conversation thread"""
    
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/webp,*/*;q=0.8',
        'Accept-Language': 'en-US,en;q=0.5',
        'Accept-Encoding': 'gzip, deflate, br',
        'Connection': 'keep-alive',
        'Upgrade-Insecure-Requests': '1',
        'Sec-Fetch-Dest': 'document',
        'Sec-Fetch-Mode': 'navigate',
        'Sec-Fetch-Site': 'none',
        'Cache-Control': 'max-age=0'
    }
    
    print(f"Fetching URL: {url}")
    
    try:
        response = requests.get(url, headers=headers, timeout=30)
        response.raise_for_status()
        
        # Save raw HTML
        html_file = f"grok_share_{thread_id}_source.html"
        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(response.text)
        print(f"✓ Saved raw HTML: {html_file} ({len(response.text)} bytes)")
        
        # Parse with BeautifulSoup
        soup = BeautifulSoup(response.text, 'html.parser')
        
        # Extract text content
        # Method 1: Try to find conversation container
        conversation_text = []
        
        # Look for common Grok share patterns
        for elem in soup.find_all(['div', 'article', 'section', 'p', 'pre', 'code']):
            text = elem.get_text(strip=True)
            if len(text) > 50:  # Ignore very short snippets
                conversation_text.append(text)
        
        # Method 2: Extract all text
        full_text = soup.get_text(separator='\n', strip=True)
        
        # Save extracted content
        content_file = f"grok_share_{thread_id}_content.txt"
        with open(content_file, 'w', encoding='utf-8') as f:
            f.write(full_text)
        print(f"✓ Saved extracted text: {content_file} ({len(full_text)} bytes)")
        
        # Try to extract structured data
        script_tags = soup.find_all('script', type='application/json')
        if script_tags:
            json_file = f"grok_share_{thread_id}_data.json"
            with open(json_file, 'w', encoding='utf-8') as f:
                for script in script_tags:
                    try:
                        data = json.loads(script.string)
                        json.dump(data, f, indent=2)
                        f.write('\n---\n')
                    except:
                        pass
            print(f"✓ Saved JSON data: {json_file}")
        
        return {
            'success': True,
            'html_size': len(response.text),
            'text_size': len(full_text),
            'html_file': html_file,
            'content_file': content_file
        }
        
    except requests.exceptions.RequestException as e:
        print(f"✗ Error fetching URL: {e}")
        return {'success': False, 'error': str(e)}

if __name__ == "__main__":
    url = "https://x.com/i/grok/share/e3a0c4534d01419f95020d8393cf0023"
    thread_id = "e3a0c4534d01419f95020d8393cf0023"
    
    result = scrape_grok_thread(url, thread_id)
    
    if result['success']:
        print(f"\n✅ SUCCESS")
        print(f"HTML: {result['html_size']} bytes → {result['html_file']}")
        print(f"Text: {result['text_size']} bytes → {result['content_file']}")
    else:
        print(f"\n❌ FAILED: {result['error']}")
