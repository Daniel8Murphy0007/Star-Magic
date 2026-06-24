# UQFF Star-Magic v5.27.2 — Multi-namespace CLI discovery

**Date:** 2026-06-22
**Type:** Patch release (CLI enhancements, no breaking changes, no physics changes)
**Upgrade:** `pip install --upgrade uqff`

---

## What changed

### `uqff search <substr>` — search ALL closure namespaces

v5.27.1's `uqff list/predict` only saw the 794 `PARADOX_TO_CLOSURE` dispatch
keys. Many headline closures (Holmlid LENR, all 8 Clay Millennium prize
problems, nuclear magic numbers, 248 bucket observables) live in **other**
name spaces and were invisible to the CLI. v5.27.2 fixes this.

```bash
uqff search holmlid          # now finds holmlid_D_minus_1
uqff search millennium       # finds all 8 Clay Millennium aliases
uqff search reactor          # finds star_magic_reactor (COP 555:1)
uqff search alpha            # finds 10 hits across 7 namespaces
```

### `uqff predict <name>` — falls back through 5 sources

```bash
uqff predict holmlid_D_minus_1        # 630 eV LENR from calculate_lenr_full
uqff predict yang_mills               # 5970 GeV mass gap from Millennium
uqff predict magic_numbers            # all 7 magic numbers from nuclear
uqff predict hubble_tension           # 67.4 km/s/Mpc from paradox dispatch
```

Search order: PARADOX_TO_CLOSURE → calculate_lenr_full → calculate_nuclear_magic
→ bucket_observables. Returns the first match with its source labeled.

### `uqff list --all` — list across all 5 namespaces

```bash
uqff list --all --filter holmlid     # filtered across all sources
uqff list --all                      # ~1,300 total closure names
uqff list                            # just paradox dispatch (default)
```

### The 5 namespaces now exposed

| Namespace | Count | Examples |
|---|---:|---|
| `PARADOX_TO_CLOSURE` | 794 | `hubble_tension`, `strong_cp_problem` |
| `PARADOX_TO_MILLENNIUM` (aliases) | 18 | `yang_mills`, `riemann`, `bsd` |
| `calculate_lenr_full` sub-keys | 10 | `holmlid_D_minus_1`, `star_magic_reactor`, `parkhomov_NiH` |
| `calculate_nuclear_magic` sub-keys | 10 | `magic_numbers`, `be_per_a_peak_MeV`, `alpha_binding_MeV` |
| Bucket observables (9 surfaces) | 248 | `alpha (fine structure)`, `m_W (W boson) GeV` |
| **TOTAL discoverable** | **~1,080** | (was 794 in v5.27.1) |

### Backward compatibility

`uqff predict <paradox_key>` still works exactly as before — the new
multi-source fallback only kicks in if the paradox dispatcher returns
`None`. No existing scripts break.

---

## Upgrade

```bash
pip install --upgrade uqff
uqff version            # should print 'uqff 5.27.2'
uqff search holmlid     # confirm Holmlid is now discoverable
```

---

## What didn't change

- All 9 truly-independent primitives, 794 closures, 857/857 fidelity gate, ΔBIC = 94.1
- Physics is identical to v5.27.0 / v5.27.1
- Public API surfaces (34 `calculate_*` functions) unchanged
- License unchanged (AGPL-3.0 + commercial)

---

**Author:** Daniel T. Murphy / Star-Magic Research Program
**Copyright:** © 2025-2026, all rights reserved.
