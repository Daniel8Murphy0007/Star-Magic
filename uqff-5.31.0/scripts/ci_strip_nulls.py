#!/usr/bin/env python3
"""Strip trailing null bytes from uqff_fidelity_tests.py.

Some Windows edit/save cycles append null padding to the test file. This
script removes them in place. Idempotent and safe to run on any platform.

Usage:  python scripts/ci_strip_nulls.py [path]
        (defaults to uqff_fidelity_tests.py in cwd)
"""
import pathlib
import sys


def main() -> int:
    target = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "uqff_fidelity_tests.py")
    if not target.is_file():
        print(f"ERROR: {target} not found")
        return 1
    data = target.read_bytes()
    n = data.count(b"\x00")
    if n:
        target.write_bytes(data.replace(b"\x00", b""))
        print(f"stripped {n} null bytes from {target}")
    else:
        print(f"no null bytes in {target}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
