import sqlite3, shutil, tempfile, os, time, datetime

src = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\cookies.sqlite')
tmp = os.path.join(tempfile.gettempdir(), 'ck2.db')
shutil.copy2(src, tmp)
conn = sqlite3.connect(tmp)
cur = conn.cursor()
cur.execute("SELECT name, value, expiry FROM moz_cookies WHERE (host LIKE '%.x.com%' OR host='x.com') AND name IN ('auth_token','ct0')")
now = int(time.time())
for name, val, exp in cur.fetchall():
    status = 'VALID' if exp > now else 'EXPIRED'
    try:
        exp_str = datetime.datetime.fromtimestamp(exp).isoformat()
    except Exception:
        exp_str = str(exp)
    print(f'{name}: {val[:20]}... expiry_epoch={exp} [{status}]  ({exp_str})')
conn.close()
os.unlink(tmp)
