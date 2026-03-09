#!/usr/bin/env python3
"""
Load Grok public share URL in Edge (headless), wait for React to render,
then extract the fully rendered conversation text.
No authentication required — public share link.
"""
import os
import sys
import time
import json
import re

SHARE_ID   = '7b0e961fa19846bea8aed946b0650e93'
TARGET_URL = f'https://x.com/i/grok/share/{SHARE_ID}'
DRIVER_PATH = r'C:\edgedriver_win64\msedgedriver.exe'
OUTPUT_TXT  = f'grok_share_{SHARE_ID[:8]}_conversation.txt'
OUTPUT_JSON = f'grok_share_{SHARE_ID[:8]}_state.json'

try:
    from selenium import webdriver
    from selenium.webdriver.edge.service import Service
    from selenium.webdriver.edge.options import Options
    from selenium.webdriver.common.by import By
    from selenium.webdriver.support.ui import WebDriverWait
    from selenium.webdriver.support import expected_conditions as EC
except ImportError:
    print('ERROR: selenium not installed. Run: pip install selenium')
    sys.exit(1)

opts = Options()
opts.add_argument('--headless')
opts.add_argument('--no-sandbox')
opts.add_argument('--disable-dev-shm-usage')
opts.add_argument('--window-size=1920,1080')
opts.add_argument('--disable-gpu')
opts.add_argument('--disable-blink-features=AutomationControlled')
opts.add_experimental_option('excludeSwitches', ['enable-automation'])

svc = Service(DRIVER_PATH)
print(f'Launching Edge (headless)...')
driver = webdriver.Edge(service=svc, options=opts)

conversation_data = None

try:
    print(f'Loading: {TARGET_URL}')
    driver.get(TARGET_URL)

    # Wait up to 30s for conversation content to appear
    print('Waiting for React to render conversation...')
    for attempt in range(30):
        time.sleep(1)
        # Check if grokShare state has populated
        result = driver.execute_script("""
            try {
                var state = window.__INITIAL_STATE__;
                if (!state) return null;
                var gs = state.entities && state.entities.grokShare;
                if (gs && gs.entities && Object.keys(gs.entities).length > 0) {
                    return JSON.stringify(gs);
                }
                // Also check Redux store if present
                if (window.__reduxStore) {
                    var s = window.__reduxStore.getState();
                    var gs2 = s && s.entities && s.entities.grokShare;
                    if (gs2 && Object.keys(gs2.entities || {}).length > 0) {
                        return JSON.stringify(gs2);
                    }
                }
                return 'waiting:' + JSON.stringify(Object.keys(state.entities || {}));
            } catch(e) { return 'err:' + e.toString(); }
        """)
        if result and not result.startswith('waiting') and not result.startswith('err'):
            print(f'grokShare state populated after {attempt+1}s')
            conversation_data = json.loads(result)
            break
        if attempt % 5 == 0:
            print(f'  t={attempt+1}s: {str(result)[:120]}')

    if not conversation_data:
        print('grokShare state did not populate — trying DOM extraction...')

    # ── DOM extraction: grab visible text from conversation elements ──────────
    print('\n=== DOM TEXT EXTRACTION ===')
    # Common Grok share page selectors
    selectors_to_try = [
        '[data-testid="grok-share-message"]',
        '[data-testid*="grok"]',
        '.grok-share',
        '[class*="grok"]',
        'article',
        'main',
    ]
    for sel in selectors_to_try:
        try:
            elems = driver.find_elements(By.CSS_SELECTOR, sel)
            if elems:
                print(f'Selector "{sel}": {len(elems)} elements found')
                for i, e in enumerate(elems[:5]):
                    txt = e.text.strip()
                    if txt:
                        print(f'  [{i}] {txt[:200]}')
        except Exception:
            pass

    # ── Capture any XHR/fetch responses via Performance Logs ─────────────────
    print('\n=== PAGE SOURCE SNIPPET ===')
    src = driver.page_source
    print(f'Page source length: {len(src):,} chars')

    # Try to find conversation text in page source
    patterns = [
        r'"message"\s*:\s*"([^"]{50,})"',
        r'"text"\s*:\s*"([^"]{50,})"',
        r'"content"\s*:\s*"([^"]{50,})"',
        r'"response"\s*:\s*"([^"]{50,})"',
        r'"query"\s*:\s*"([^"]{20,})"',
    ]
    found_texts = set()
    for pat in patterns:
        for m in re.finditer(pat, src):
            t = m.group(1).replace('\\n', '\n').replace('\\"', '"')
            if len(t) > 50 and t not in found_texts:
                found_texts.add(t)
                print(f'\nPattern [{pat[:30]}] match:')
                print(f'  {t[:300]}')

    # Save page source for manual inspection
    with open(f'grok_share_{SHARE_ID[:8]}_source.html', 'w', encoding='utf-8') as f:
        f.write(src)
    print(f'\nSaved full page source: grok_share_{SHARE_ID[:8]}_source.html')

    # ── Final state dump ──────────────────────────────────────────────────────
    full_state = driver.execute_script("""
        try { return JSON.stringify(window.__INITIAL_STATE__ || {}); }
        catch(e) { return '{}'; }
    """)
    if full_state and len(full_state) > 10:
        state_obj = json.loads(full_state)
        with open(OUTPUT_JSON, 'w', encoding='utf-8') as f:
            json.dump(state_obj, f, indent=2, ensure_ascii=False)
        print(f'Saved runtime state: {OUTPUT_JSON}')
        gs = state_obj.get('entities', {}).get('grokShare', {}).get('entities', {})
        if gs:
            print(f'\n=== GROK SHARE ENTITIES ({len(gs)} items) ===')
            print(json.dumps(gs, indent=2, ensure_ascii=False)[:3000])
        else:
            print('\ngrokShare.entities still empty after render.')
            print(f'entities keys: {list(state_obj.get("entities", {}).keys())}')

finally:
    driver.quit()
    print('\nEdge closed.')
