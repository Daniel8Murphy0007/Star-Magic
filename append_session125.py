"""Update all tracking files for Session 125 / v4.98."""

SESSION125_TRACKER_ROW = """| ✅ Session 125 | **2 new whitepapers PAPER_479–480 — v4.98 grok_share_4e4d8be1f7.txt (2,327 lines): PAPER_479 UQFFBuoyancyAstroModule 5-system complex-arithmetic buoyancy (J1610+1811 z=3.122, PLCK G287.0+32.9, PSZ2 G181.06+48.47, ASKAP J1832-0911, Sonification Collection — FIRST); cdouble throughout; F_rel=4.30e33 N (LEP 1998); integral≈Integrand×x₂ quadratic root approx; LENR dominant at low ω₀; PAPER_480 UQFFBuoyancyCNBModule 6-system + CentaurusA — FIRST CNB neutrino buoyancy F_ν=k_ν×σ_CNB×n_CNB×E_CNB≈9.07×10⁻⁴² N (smallest UQFF force); Sweet vacuum F_Sweet; Kozima TNCF drop; CNB gravitational overdensity δ_CNB~O(1) at CenA SMBH; 3 C++ modules extracted (UQFFBuoyancyModule + UQFFBuoyancyAstroModule + UQFFBuoyancyCNBModule — 6 .h+.cpp files from Source161/162); INTEGRATION_PLAN_4e4d8be1f7.md; CP4 103 unchanged; CP3 219 unchanged; CP2 600 unchanged; 478→480/1000 (48.0%) commit v4.98 ✅** |
|---|"""

VERSION_HISTORY_ROW = "| v2.4.0 | Session 125 | 2026-03-23 | PAPER_479-480 (2 whitepapers): UQFFBuoyancyAstro 5-system complex arithmetic + UQFFBuoyancyCNB 6-system CNB neutrino coupling; 3 C++ modules (6 files); INTEGRATION_PLAN_4e4d8be1f7.md; 480/1000 (48.0%); commit v4.98 |"

# 1. Update VALIDATION_MASTER_INDEX_2.md
with open('VALIDATION_MASTER_INDEX_2.md', encoding='utf-8') as f:
    content = f.read()
    lines = content.splitlines()

# Insert Session 125 row BEFORE Session 124 row (line 118)
insert_before = '| ✅ Session 124 |'
new_lines = []
inserted_tracker = False
for line in lines:
    if not inserted_tracker and line.startswith(insert_before):
        new_lines.append(SESSION125_TRACKER_ROW)
        inserted_tracker = True
    new_lines.append(line)

# Add VERSION HISTORY entry after v2.3.0 row
insert_after_ver = '| v2.3.0 | Session 124 |'
final_lines = []
for line in new_lines:
    final_lines.append(line)
    if line.startswith(insert_after_ver):
        final_lines.append(VERSION_HISTORY_ROW)

with open('VALIDATION_MASTER_INDEX_2.md', 'w', encoding='utf-8') as f:
    f.write('\n'.join(final_lines) + '\n')

print('VMI2 updated: Session 125 tracker row + v2.4.0 VERSION HISTORY')
if inserted_tracker:
    print('  Tracker row inserted successfully')
else:
    print('  WARNING: Tracker row not inserted (marker not found)')

# 2. Update VALIDATION_COMPARISON_REPORT.md — update status footer
with open('VALIDATION_COMPARISON_REPORT.md', encoding='utf-8') as f:
    vcr = f.read()

old_state = "**Last Synced: March 2026 (Session 123, v4.96, PAPER_472–478, 478/1000 papers)**"
new_state = "**Last Synced: March 2026 (Session 125, v4.98, PAPER_479–480, 480/1000 papers)**"
old_state_line = "**State:** CP1 = 1,227 classes, CP2 = 600 classes, CP3 = 219 classes, CP4 = 103 classes, Aggregator updated (v4.93), VMI2 v4.96, 478/1000 papers (415 in whitepapers/), 383/383 QS=5 (prior papers), commit v4.96"
new_state_line = "**State:** CP1 = 1,227 classes, CP2 = 600 classes, CP3 = 219 classes, CP4 = 103 classes, Aggregator updated (v4.93), VMI2 v4.98, 480/1000 papers (417 in whitepapers/), 383/383 QS=5 (prior papers), commit v4.98"

vcr2 = vcr.replace(old_state, new_state)
vcr3 = vcr2.replace(old_state_line, new_state_line)

with open('VALIDATION_COMPARISON_REPORT.md', 'w', encoding='utf-8') as f:
    f.write(vcr3)

changes = int(old_state not in vcr3) + int(old_state_line not in vcr3)
print(f'VCR updated: {changes}/2 replacements made')

# 3. Update ipc_pipeline_handler.h — add Session 125 trigger keywords
with open('ipc_pipeline_handler.h', encoding='utf-8') as f:
    hdr = f.read()

old_keywords = ' *   "UaVacuumDensity", "ScmVelocity" (Session 124 — grok_share_b0a3dc1d, 48 module .h+.cpp implementations)'
new_keywords = (' *   "UaVacuumDensity", "ScmVelocity" (Session 124 — grok_share_b0a3dc1d, 48 module .h+.cpp implementations)\n'
                ' *   "UQFFBuoyancy", "UQFFBuoyancyAstro", "UQFFBuoyancyCNB", "J1610Quasar", "PLCKCluster",\n'
                ' *   "PSZ2MergingCluster", "ASKAPTransient", "SonificationCollection", "CentaurusACNB",\n'
                ' *   "CNBNeutrinoCoupling", "SweetVacuum", "KozimaDropTNCF", "QuadraticRootIntegral", "LENRDominant"\n'
                ' *   (Session 125 — grok_share_4e4d8be1f7; 3 UQFFBuoyancy modules; PAPER_479-480; CNB neutrino term)')

old_updated = '* Updated: 2026-03-23 (Session 124 — grok_share_b0a3dc1d; 48 module .h+.cpp implementations completed; 48 CP trigger ke'
new_updated = ('* Updated: 2026-03-23 (Session 124 — grok_share_b0a3dc1d; 48 module .h+.cpp implementations completed; 48 CP trigger ke'
               '\n * Updated: 2026-03-23 (Session 125 — grok_share_4e4d8be1f7; 3 UQFFBuoyancy modules +CNB; PAPER_479-480; 14 CP trigger keywords added)')

hdr2 = hdr.replace(old_keywords, new_keywords)

# Check if the update line is already there
if 'Session 125' not in hdr2:
    # Find the Session 124 update line and append after it
    s124_marker = ' * Updated: 2026-03-23 (Session 124'
    s124_idx = hdr2.find(s124_marker)
    if s124_idx > 0:
        end_of_line = hdr2.find('\n', s124_idx)
        hdr2 = hdr2[:end_of_line] + '\n * Updated: 2026-03-23 (Session 125 — grok_share_4e4d8be1f7; 3 UQFFBuoyancy modules (Astro+CNB); PAPER_479-480; 14 new CP trigger keywords; commit v4.98)' + hdr2[end_of_line:]

with open('ipc_pipeline_handler.h', 'w', encoding='utf-8') as f:
    f.write(hdr2)

kw_ok = 'UQFFBuoyancy' in hdr2
s125_ok = 'Session 125' in hdr2
print(f'ipc_pipeline_handler.h updated: keywords={kw_ok} session125={s125_ok}')
print('Done.')
