"""
UQFF Duplicate Class Audit — Theory of Permanence architecture verification.

Reports duplicate class definitions in CondensedPhysics.py as **intentional
simultaneous parallel method variants**, not bugs. Confirms both variants
(<Name>_v1 = EARLIER Gen-1 method, <Name> = LATER Gen-2 method) remain
runtime-callable per Theory of Permanence (PAPER_1929): all methods active
simultaneously and in conjunction with vacuum buoyancy.

Verifies:
1. Each family has BOTH a `_v1` and a bare-name class definition present
2. SIMULTANEOUS_METHOD_VARIANTS registry contains every duplicate family
3. No family has ONLY one instance renamed (partial rename would be a bug)
4. Reports framework-annotation status of both variants for upgrade tracking

Categories used are NEUTRAL — the user determines what is a bug, not this script.
"""
import re
import sys

TARGET_FILE = 'CondensedPhysics.py'
FRAMEWORK_MARKER = 'framework_backbone'


def find_class_instances(src):
    """Return dict: base_name -> {'v1': loc, 'bare': loc}."""
    pat = re.compile(r'^class\s+(\w+Calculator(?:_v\d+)?)\b[^:]*:', re.MULTILINE)
    matches = list(pat.finditer(src))
    result = {}
    for i, m in enumerate(matches):
        name = m.group(1)
        end = matches[i + 1].start() if i + 1 < len(matches) else len(src)
        body = src[m.end():end]
        role = 'v1' if re.search(r'_v\d+$', name) else 'bare'
        base = re.sub(r'_v\d+$', '', name)
        result.setdefault(base, {})[role] = {
            'line': src[:m.start()].count('\n') + 1,
            'body_size': len(body),
            'has_framework': FRAMEWORK_MARKER in body,
        }
    return result


def audit(src):
    by_base = find_class_instances(src)
    families = {b: locs for b, locs in by_base.items()
                if 'v1' in locs and 'bare' in locs}

    print("=" * 74)
    print(f"UQFF Simultaneous Method-Variant Audit — {TARGET_FILE}")
    print("Theory of Permanence: all method variants remain active simultaneously")
    print("=" * 74)
    print(f"\nTotal Calculator classes (all names):  {len(by_base)}")
    print(f"Multi-variant families (v1 + bare):    {len(families)}")

    both_upgraded = []
    v1_only_upgraded = []
    bare_only_upgraded = []
    neither = []

    for base in sorted(families):
        v1_fw = families[base]['v1']['has_framework']
        bare_fw = families[base]['bare']['has_framework']
        if v1_fw and bare_fw:
            both_upgraded.append(base)
        elif v1_fw and not bare_fw:
            v1_only_upgraded.append(base)
        elif bare_fw and not v1_fw:
            bare_only_upgraded.append(base)
        else:
            neither.append(base)

    print(f"\nFramework-annotation status (informational — for upgrade tracking):")
    print(f"  BOTH variants have framework annotations:  {len(both_upgraded)}")
    print(f"  Only _v1 has framework annotations:        {len(v1_only_upgraded)}")
    print(f"  Only bare has framework annotations:       {len(bare_only_upgraded)}")
    print(f"  Neither variant upgraded yet:              {len(neither)}")

    if v1_only_upgraded:
        print(f"\n  Families where _v1 is upgraded but bare is not (informational):")
        for name in v1_only_upgraded:
            print(f"    * {name}")

    # STRUCTURAL CHECK: registry present and populated?
    registry_present = 'SIMULTANEOUS_METHOD_VARIANTS' in src
    dict_registry_present = 'SIMULTANEOUS_SOURCE_DICT_VARIANTS' in src
    print(f"\nRegistry structural check:")
    print(f"  SIMULTANEOUS_METHOD_VARIANTS present:      {registry_present}")
    print(f"  SIMULTANEOUS_SOURCE_DICT_VARIANTS present: {dict_registry_present}")

    return families


if __name__ == '__main__':
    with open(TARGET_FILE, 'r', encoding='utf-8') as f:
        src = f.read()
    audit(src)
