"""Try to rebuild a paper with the standard prep pipeline. Reports OK/FAIL only."""
import subprocess, sys, os
name = sys.argv[1]
md = f"whitepapers/{name}.md"
if not os.path.exists(md):
    print(f"NOTFOUND {name}"); sys.exit(2)
for script in ("_strip_emoji.py", "_fix_cmd_letter_glue.py"):
    subprocess.run([sys.executable, script, md], check=False, capture_output=True)
r = subprocess.run([sys.executable, "_regen_one.py", name], capture_output=True, text=True, encoding="utf-8", errors="replace")
out = (r.stdout or "") + (r.stderr or "")
if r.returncode == 0:
    print(f"OK   {name}")
else:
    # find first error line
    err_line = ""
    for line in out.splitlines():
        if line.startswith("!") or line.startswith("l."):
            err_line += line[:120] + " | "
            if err_line.count("|") >= 3: break
    print(f"FAIL {name}  {err_line}")
sys.exit(0)
