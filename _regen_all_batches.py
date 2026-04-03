"""
_regen_all_batches.py
Regenerate all PAPER_*.pdf in sequential batches of 20 (workers=2).
Runs generate_pdfs.py START END for each batch.
"""
import subprocess, sys, time, glob, re, os

BATCH_SIZE = 20

# Discover all paper numbers from source files
md_files = sorted(glob.glob("whitepapers/PAPER_*.md"))
nums = []
for f in md_files:
    m = re.search(r"PAPER_(\d+)", os.path.basename(f))
    if m:
        nums.append(int(m.group(1)))
nums = sorted(set(nums))

if not nums:
    print("No whitepapers found.")
    sys.exit(1)

total_papers = len(nums)
# Build batches
batches = []
for i in range(0, total_papers, BATCH_SIZE):
    chunk = nums[i:i + BATCH_SIZE]
    batches.append((chunk[0], chunk[-1]))

print(f"Papers: {nums[0]}–{nums[-1]}  ({total_papers} total)")
print(f"Batches of {BATCH_SIZE}: {len(batches)} batches\n")

ok_total   = 0
fail_total = 0
t_start    = time.time()

for idx, (start, end) in enumerate(batches, 1):
    batch_label = f"[{idx:3d}/{len(batches)}]  {start:04d}–{end:04d}"
    print(f"{batch_label}  starting...", flush=True)
    b_start = time.time()
    result = subprocess.run(
        [sys.executable, "generate_pdfs.py", str(start), str(end)],
        capture_output=False,   # let output stream live
    )
    b_elapsed = time.time() - b_start
    status = "OK" if result.returncode == 0 else f"EXIT={result.returncode}"
    print(f"{batch_label}  done in {b_elapsed:.0f}s  [{status}]\n", flush=True)

elapsed = time.time() - t_start
m, s = divmod(int(elapsed), 60)
print("=" * 60)
print(f"  All {len(batches)} batches complete  ({total_papers} papers)  {m}m {s}s")
print("=" * 60)
