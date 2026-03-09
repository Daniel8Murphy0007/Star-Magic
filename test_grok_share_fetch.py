#!/usr/bin/env python3
"""
Test scrape method on xAI Grok share URL.
Tries 3 approaches in order:
  1. Unauthenticated page GET
  2. GraphQL GrokShare endpoint with xAI API key as Bearer
  3. xAI API - check if any /v1/ endpoint exposes shared conversations
"""
import requests
import json
import os

SHARE_ID = '7b0e961fa19846bea8aed946b0650e93'
SHARE_URL = f'https://x.com/i/grok/share/{SHARE_ID}'
XAI_KEY   = os.environ.get('XAI_API_KEY', '')  # set via environment variable

HEADERS_BROWSER = {
    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/123.0 Safari/537.36',
    'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8',
    'Accept-Language': 'en-US,en;q=0.9',
}

# ─────────────────────────────────────────────────────────────────
# TEST 1: Unauthenticated page fetch
# ─────────────────────────────────────────────────────────────────
print('=' * 60)
print('TEST 1: Unauthenticated GET of share page')
print('=' * 60)
try:
    r = requests.get(SHARE_URL, timeout=15, headers=HEADERS_BROWSER, allow_redirects=True)
    print(f'Status      : {r.status_code}')
    print(f'Final URL   : {r.url}')
    ct = r.headers.get('Content-Type', '?')
    print(f'Content-Type: {ct}')
    body = r.text
    print(f'Body length : {len(body)} chars')
    # Look for any meaningful content
    if 'grok' in body.lower() or 'conversation' in body.lower():
        print('>>> Content contains Grok/conversation markers')
    if '<title>' in body:
        import re
        title = re.search(r'<title>(.*?)</title>', body, re.IGNORECASE)
        if title:
            print(f'Page title  : {title.group(1)}')
    print(f'Body snippet: {body[:400]}')
except Exception as e:
    print(f'ERROR: {e}')

print()

# ─────────────────────────────────────────────────────────────────
# TEST 2: GrokShare GraphQL endpoint with xAI key as Bearer
# ─────────────────────────────────────────────────────────────────
print('=' * 60)
print('TEST 2: GrokShare GraphQL with xAI API key as Bearer')
print('=' * 60)
variables = json.dumps({"grok_share_id": SHARE_ID})
features  = json.dumps({
    "premium_content_api_read_enabled": False,
    "communities_web_enable_tweet_community_results_fetch": True,
})
graphql_url = (
    'https://x.com/i/api/graphql/3aCm_HRrYXX8T7sas50Zlw/GrokShare'
    f'?variables={requests.utils.quote(variables)}'
    f'&features={requests.utils.quote(features)}'
)
try:
    r2 = requests.get(graphql_url, timeout=15, headers={
        'Authorization': f'Bearer {XAI_KEY}',
        'Content-Type': 'application/json',
        'x-twitter-active-user': 'yes',
        'x-twitter-client-language': 'en',
        'User-Agent': HEADERS_BROWSER['User-Agent'],
    })
    print(f'Status      : {r2.status_code}')
    ct2 = r2.headers.get('Content-Type', '?')
    print(f'Content-Type: {ct2}')
    print(f'Body snippet: {r2.text[:500]}')
except Exception as e:
    print(f'ERROR: {e}')

print()

# ─────────────────────────────────────────────────────────────────
# TEST 3: xAI API base endpoint check
# ─────────────────────────────────────────────────────────────────
print('=' * 60)
print('TEST 3: xAI API /v1/models (verify key validity)')
print('=' * 60)
try:
    r3 = requests.get('https://api.x.ai/v1/models', timeout=15, headers={
        'Authorization': f'Bearer {XAI_KEY}',
        'Content-Type': 'application/json',
    })
    print(f'Status      : {r3.status_code}')
    data = r3.json()
    if 'data' in data:
        models = [m.get('id', '?') for m in data['data']]
        print(f'Available models: {models}')
    else:
        print(f'Response: {json.dumps(data, indent=2)[:400]}')
except Exception as e:
    print(f'ERROR: {e}')

print()
print('=' * 60)
print('SUMMARY')
print('=' * 60)
print('If TEST 1 returned 200 with content -> page is public, no auth needed')
print('If TEST 2 returned 200 with JSON    -> xAI key works on GrokShare GraphQL')
print('If TEST 2 returned 403/401          -> need Firefox cookie injection method')
print('If TEST 3 returned 200 + models     -> xAI key is valid')
