"""Strip emoji & misc symbols safely (Python, not PowerShell)."""
import sys, pathlib, re
fp = pathlib.Path(sys.argv[1])
s = fp.read_text(encoding='utf-8')
# Remove characters in supplementary planes (most emoji)
s = ''.join(c for c in s if ord(c) < 0x10000)
# Remove BMP emoji ranges
s = re.sub(r'[\u2600-\u27BF]', '', s)
# Common single-char emoji & icons
for c in '\u2705\u274C\u2713\u2717\u26A0\u2728\u2B50\u2B07\u2B06\u2B05\u2B95':
    s = s.replace(c, '')
fp.write_text(s, encoding='utf-8')
print('ok')
