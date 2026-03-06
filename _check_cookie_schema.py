import ctypes
import ctypes.wintypes as wt
import sqlite3
import tempfile
import os

GENERIC_READ = 0x80000000
FILE_SHARE_READ  = 0x00000001
FILE_SHARE_WRITE = 0x00000002
FILE_SHARE_DELETE = 0x00000004
OPEN_EXISTING = 3
FILE_FLAG_SEQUENTIAL_SCAN = 0x08000000

k32 = ctypes.WinDLL('kernel32', use_last_error=True)
db = r'C:\Users\tmsjd\AppData\Local\Microsoft\Edge\User Data\Default\Network\Cookies'

h = k32.CreateFileW(db, GENERIC_READ, FILE_SHARE_READ|FILE_SHARE_WRITE|FILE_SHARE_DELETE, None, OPEN_EXISTING, FILE_FLAG_SEQUENTIAL_SCAN, None)
sz = wt.LARGE_INTEGER()
k32.GetFileSizeEx(h, ctypes.byref(sz))
buf = (ctypes.c_char * sz.value)()
rd = wt.DWORD()
k32.ReadFile(h, buf, sz.value, ctypes.byref(rd), None)
k32.CloseHandle(h)

tmp = os.path.join(tempfile.gettempdir(), 'ck_tmp2.db')
with open(tmp, 'wb') as f:
    f.write(bytes(buf[:rd.value]))

conn = sqlite3.connect(tmp)
tables = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'").fetchall()]
print("Tables:", tables)
for t in tables:
    cols = [r[1] for r in conn.execute(f"PRAGMA table_info({t})").fetchall()]
    print(f"  {t}: {cols}")
conn.close()
os.remove(tmp)
