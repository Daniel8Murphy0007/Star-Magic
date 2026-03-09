import re

with open('grok_share_ab2d0965e3a74a0da32749a7a2591dc7_source.html', encoding='utf-8') as f:
    html = f.read()

idx = html.find('DanielMurp')
print('DanielMurp context:', html[max(0,idx-100):idx+400] if idx >= 0 else 'NOT FOUND')

for marker in ['conversation', 'message', 'grok', 'response', 'Human:', 'Assistant:']:
    count = html.lower().count(marker.lower())
    print(f'{marker}: {count}')

matches = re.findall(r'data-testid="[^"]{0,60}"', html[:100000])
print('Test IDs (unique):', list(set(matches))[:30])

# Find JSON data blobs
json_blobs = re.findall(r'"messages"\s*:\s*\[', html)
print(f'messages arrays found: {len(json_blobs)}')

# Look for conversation thread content
idxc = html.find('Copy Conversation')
print(f'\nCopy Conversation context:\n{html[max(0,idxc-200):idxc+500]}')
