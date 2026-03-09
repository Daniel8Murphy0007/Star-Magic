import sqlite3, shutil, tempfile, os

src = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\cookies.sqlite')
tmp = os.path.join(tempfile.gettempdir(), 'ck3.db')
shutil.copy2(src, tmp)
conn = sqlite3.connect(tmp)
cur = conn.cursor()
cur.execute("SELECT host, name, value, path, expiry, isSecure, isHttpOnly, sameSite FROM moz_cookies WHERE (host LIKE '%.x.com%' OR host='x.com') AND name IN ('auth_token','ct0')")
print(f"{'name':<15} {'host':<15} {'httpOnly':<10} {'secure':<8} {'sameSite':<10} value[:20]")
print('-'*80)
for host, name, value, path, expiry, secure, http_only, same_site in cur.fetchall():
    print(f"{name:<15} {host:<15} {str(bool(http_only)):<10} {str(bool(secure)):<8} {str(same_site):<10} {value[:20]}")
conn.close()
os.unlink(tmp)
