#!/usr/bin/env python3
"""
Parse Grok share HTML page source and extract conversation text.
"""
import re, sys, os
from html.parser import HTMLParser

SHARE_ID = '7b78ffcb915f48bb90d55034c9c50b18'
INPUT    = f'grok_share_{SHARE_ID}_content.txt'
OUTPUT   = f'grok_share_{SHARE_ID}_parsed.txt'

with open(INPUT, 'r', encoding='utf-8') as f:
    html = f.read()

print(f'HTML source: {len(html)} chars')

# Strategy 1: look for the JSON data blob embedded in the page (SSR / hydration)
# X.com embeds conversation data as JSON in a script tag
json_patterns = [
    r'"grok_share":\s*(\{[^}]+\})',
    r'"conversation":\s*(\[.+?\])',
    r'"messages":\s*(\[.+?\])',
    r'grokConversation["\s:]+(\{.+?\})',
]
for pat in json_patterns:
    m = re.search(pat, html, re.DOTALL)
    if m:
        print(f'Found JSON data via pattern: {pat[:40]}')
        print(m.group(0)[:300])
        print('---')

# Strategy 2: extract all visible text by stripping HTML tags
class TextExtractor(HTMLParser):
    def __init__(self):
        super().__init__()
        self.texts = []
        self.skip  = False
        self.skip_tags = {'script', 'style', 'noscript', 'head'}
        self.depth = 0
    def handle_starttag(self, tag, attrs):
        if tag in self.skip_tags:
            self.skip = True
            self.depth += 1
    def handle_endtag(self, tag):
        if tag in self.skip_tags:
            self.depth -= 1
            if self.depth <= 0:
                self.skip = False
                self.depth = 0
    def handle_data(self, data):
        if not self.skip:
            text = data.strip()
            if len(text) > 2:
                self.texts.append(text)

extractor = TextExtractor()
extractor.feed(html)
all_text = '\n'.join(extractor.texts)
print(f'\nExtracted text: {len(all_text)} chars')
print(f'Lines: {len(extractor.texts)}')

# Strategy 3: Find the conversation content — look for long paragraphs
# (messages tend to be multi-sentence paragraphs, not nav labels)
long_texts = [t for t in extractor.texts if len(t) > 80]
print(f'Long text segments (>80 chars): {len(long_texts)}')

if long_texts:
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write('\n\n'.join(long_texts))
    print(f'\nSaved long-text segments → {OUTPUT}')
    print('\n=== FIRST 5000 CHARS ===')
    print('\n\n'.join(long_texts)[:5000])
else:
    # Save all extracted text anyway
    with open(OUTPUT, 'w', encoding='utf-8') as f:
        f.write(all_text)
    print(f'Saved all text → {OUTPUT}')
    print(all_text[:3000])
