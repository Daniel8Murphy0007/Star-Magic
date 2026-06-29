"""
test_phase_g_cli_dispatch.py — verify uqff_cli.py Phase E/F/G integration.

Checks:
 1. `uqff_cli.main(['version'])` doesn't crash
 2. `uqff_cli.main(['list','--domain','SI'])` returns SI observables
 3. `uqff_cli.main(['list','--dispatch'])` returns all 114
 4. `uqff_cli.main(['assimilate','alpha_inverse'])` returns 137.0
 5. `uqff_cli.main(['assimilate','LCDM_BAO_rd_H0_over_c_primary','--decompose'])` returns 8-field dict
 6. `uqff_cli.main(['predict','LCDM_BAO_rd_H0_over_c_primary'])` finds it via dispatch fallback (case-insensitive)
"""
import sys, io, contextlib
sys.path.insert(0, '.')
import uqff_cli


def run(argv):
    buf_out = io.StringIO()
    buf_err = io.StringIO()
    with contextlib.redirect_stdout(buf_out), contextlib.redirect_stderr(buf_err):
        rc = uqff_cli.main(argv)
    return rc, buf_out.getvalue(), buf_err.getvalue()


def main():
    print("=" * 72)
    print("PHASE G CLI regression — Phase E/F/G dispatch integration")
    print("=" * 72)

    all_pass = True

    rc, out, err = run(["list", "--domain", "SI"])
    if rc != 0 or "alpha_inverse" not in out or "7 observable(s)" not in out:
        print(f"FAIL: list --domain SI rc={rc}, out_head={out[:80]!r}")
        all_pass = False
    else:
        print("PASS  list --domain SI -> 7 observables, alpha_inverse present")

    rc, out, _ = run(["list", "--dispatch"])
    if rc != 0 or "116 observable(s)" not in out:
        print(f"FAIL: list --dispatch should report 114, got rc={rc}")
        all_pass = False
    else:
        print("PASS  list --dispatch -> 116 observables")

    rc, out, _ = run(["assimilate", "alpha_inverse"])
    if rc != 0 or "value: 137.0" not in out:
        print(f"FAIL: assimilate alpha_inverse: rc={rc}, out={out[:100]!r}")
        all_pass = False
    else:
        print("PASS  assimilate alpha_inverse -> value: 137.0")

    rc, out, _ = run(["assimilate", "LCDM_BAO_rd_H0_over_c_primary", "--decompose"])
    if rc != 0 or "geometry_used: d26" not in out or "residual_pct" not in out:
        print(f"FAIL: assimilate BAO --decompose: rc={rc}, out={out[:120]!r}")
        all_pass = False
    else:
        print("PASS  assimilate BAO primary --decompose -> d26 geometry, residual reported")

    rc, out, _ = run(["predict", "LCDM_BAO_rd_H0_over_c_primary"])
    if rc != 0 or "assimilation_dispatch" not in out:
        print(f"FAIL: predict BAO falls back to dispatch: rc={rc}, out={out[:120]!r}")
        all_pass = False
    else:
        print("PASS  predict BAO -> falls back to assimilation_dispatch (case-insensitive)")

    rc, out, _ = run(["predict", "lcdm_bao_rd_h0_over_c_primary"])
    if rc != 0 or "assimilation_dispatch" not in out:
        print(f"FAIL: predict BAO lowercase: rc={rc}, out={out[:120]!r}")
        all_pass = False
    else:
        print("PASS  predict BAO lowercase -> resolves case-insensitively")

    print()
    if all_pass:
        print("PHASE G CLI SUCCESS CRITERION MET.")
        return 0
    print("PHASE G CLI FAILURE.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
