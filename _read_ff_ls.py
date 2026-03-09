import sqlite3, shutil, tempfile, os, json

src = os.path.expandvars(r'%APPDATA%\Mozilla\Firefox\Profiles\pawatwum.default-release\storage\default\https+++x.com\ls\data.sqlite')
tmp = os.path.join(tempfile.gettempdir(), 'ff_ls.db')
shutil.copy2(src, tmp)

conn = sqlite3.connect(tmp)
cur = conn.cursor()

# List tables
cur.execute("SELECT name FROM sqlite_master WHERE type='table'")
print('Tables:', cur.fetchall())

# Try to read all rows
try:
    rows = cur.execute("SELECT * FROM data").fetchall()
    print(f'\n{len(rows)} rows in data table')
    for row in rows:
        print('ROW:', str(row)[:200])
except Exception as e:
    print('data error:', e)
    cur.execute("SELECT * FROM sqlite_master")
    print(cur.fetchall())

conn.close()
os.unlink(tmp)
