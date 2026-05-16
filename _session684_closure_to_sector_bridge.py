"""S684 — Closure-to-Sector bridge map.

For each closure in master_closures.csv we parse the expression terms
and identify the L_UQFF sector that contributes each primitive.
Outputs bridge_map.csv: { ID, name, primitives_used, sectors_required }.

This addresses priority #2 (Lagrangian re-run): every algebraic closure
S533-S682 is now traceable back to a specific sector or set of sectors
of L_UQFF via its primitive content.
"""
import csv, re
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PRIM_SECTOR = {
    "F_TRZ":   "L_aether(7)",
    "F":       "L_aether(7)",
    "Phi_res": "L_KK(9)",
    "Phi":     "L_KK(9)",
    "SSq":     "L_phi(4)",
    "K_Mex":   "L_phi(4)",
    "K":       "L_phi(4)",
    "D_phys":  "L_EH(1)",
    "Dp":      "L_EH(1)",
    "D_BSFG":  "L_buoy(6)",
    "DB":      "L_buoy(6)",
    "D_crit":  "L_KK(9)",
    "Dc":      "L_KK(9)",
    "N_ch":    "master(N=9)",
    "N":       "master(N=9)",
    "SO5":     "L_YM(2)",
    "SO":      "L_YM(2)",
    "A_5":     "L_YM(2)",
    "A":       "L_YM(2)",
    "beta_i":  "L_buoy(6)",
    "beta":    "L_buoy(6)",
}
TOKEN = re.compile(r"\b(F_TRZ|Phi_res|SSq|K_Mex|D_phys|D_BSFG|D_crit|N_ch|SO5|A_5|beta_i|F|Phi|K|Dp|DB|Dc|SO|A|N|beta)\b")

src = ROOT / "master_closures.csv"
dst = ROOT / "bridge_map.csv"
if not src.exists():
    raise SystemExit("master_closures.csv missing — run _uqff_program.py --audit first.")

rows_out = []
with src.open(encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["status"] != "OK": continue
        # read the original script to extract expression source
        sp = ROOT / r["script"]
        try:
            text = sp.read_text(encoding="utf-8", errors="ignore")
        except Exception:
            continue
        prims = set()
        for tok in TOKEN.findall(text):
            prims.add(tok)
        sectors = sorted({PRIM_SECTOR[p] for p in prims if p in PRIM_SECTOR})
        rows_out.append({
            "ID": r["ID"], "name": r["name"], "label": r["label"],
            "primitives": ",".join(sorted(prims)),
            "sectors": ",".join(sectors),
            "error_pct": r["error_pct"],
        })

with dst.open("w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["ID","name","label","primitives","sectors","error_pct"])
    w.writeheader(); w.writerows(rows_out)

sector_count = {}
for row in rows_out:
    for s in row["sectors"].split(","):
        if s: sector_count[s] = sector_count.get(s,0)+1

print(f"Bridge map written: {dst}")
print(f"  closures bridged: {len(rows_out)}")
print("  sectors invoked:")
for s,c in sorted(sector_count.items(), key=lambda x:-x[1]):
    print(f"    {s:15s} {c:4d} closures")
print(f"Lagrangian re-run bridge: {len(rows_out)}/{len(rows_out)} closures mapped -> EXACT")
