import os, sqlite3, shutil, tempfile, glob

base = os.path.expandvars(r'%LOCALAPPDATA%\Microsoft\Edge\User Data')
profiles = ['Default'] + [os.path.basename(p) for p in glob.glob(os.path.join(base, 'Profile *'))]

for prof in profiles:
    cookies_db = os.path.join(base, prof, 'Network', 'Cookies')
    if not os.path.exists(cookies_db):
        continue
    tmp = os.path.join(tempfile.gettempdir(), f'edge_chk_{prof.replace(" ","_")}.db')
    shutil.copy2(cookies_db, tmp)
    try:
        conn = sqlite3.connect(tmp)
        c = conn.cursor()
        c.execute("SELECT host_key, name FROM cookies WHERE host_key LIKE '%x.com%' OR host_key LIKE '%twitter%'")
        rows = c.fetchall()
        conn.close()
        print(f'=== {prof}: {len(rows)} cookies ===')
        for r in rows[:10]:
            print(f'  {r}')
    except Exception as e:
        print(f'  ERROR: {e}')
    finally:
        try: os.unlink(tmp)
        except: pass
