import sqlite3
import urllib.parse

db = r'C:\Users\tmsjd\AppData\Local\Microsoft\Edge\User Data\Default\Network\Cookies'

# Try: URI mode with nolock=1 (allows read of locked WAL database)
try:
    uri = 'file:' + urllib.parse.quote(db.replace('\\', '/'), safe='/:') + '?mode=ro&nolock=1'
    print(f"URI: {uri[:80]}")
    conn = sqlite3.connect(uri, uri=True, timeout=3)
    tables = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")]
    print('Tables:', tables)
    conn.close()
except Exception as e:
    print(f'URI mode error: {e}')

# Try: use immutable mode
try:
    uri2 = 'file:' + urllib.parse.quote(db.replace('\\', '/'), safe='/:') + '?immutable=1'
    conn2 = sqlite3.connect(uri2, uri=True, timeout=3)
    tables2 = [r[0] for r in conn2.execute("SELECT name FROM sqlite_master WHERE type='table'")]
    print('immutable Tables:', tables2)
    conn2.close()
except Exception as e:
    print(f'immutable mode error: {e}')
