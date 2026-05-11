"""Re-verify all 177 papers from _reverify_list.txt. Print summary."""
import subprocess, sys
names = open("_reverify_list.txt", encoding="utf-8").read().splitlines()
ok, fail = [], []
for i, n in enumerate(names, 1):
    if not n: continue
    r = subprocess.run([sys.executable, "_try_one.py", n], capture_output=True, text=True, encoding="utf-8", errors="replace")
    line = (r.stdout or "").strip().splitlines()
    line = line[-1] if line else ""
    print(f"[{i}/{len(names)}] {line}", flush=True)
    if line.startswith("OK"):
        ok.append(n)
    else:
        fail.append(line)
print(f"\nTotal OK: {len(ok)}  FAIL: {len(fail)}")
if fail:
    open("_reverify_fails.txt", "w", encoding="utf-8").write("\n".join(fail))
    print("Fails written to _reverify_fails.txt")
