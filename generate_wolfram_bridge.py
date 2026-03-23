"""
generate_wolfram_bridge.py
=========================
Generates wolfram_sources_bridge.cpp which provides:
1. Minimal-modification copies of all 705 unique system-specific PhysicsTerm
   subclasses from source*_wolfram.cpp files, placed in per-file namespaces.
2. A single registerAllWolframSourceTerms(CalculatorCore& core) function
   that registers all unique bridge terms under Batch 21.

Strategy: minimal modification approach
- Copy each unique class body verbatim
- ONLY change:
  a) "public PhysicsTerm" -> "public ::PhysicsTerm" to use the global definition
  b) "double compute(double t) const override" -> add ", const std::map<...>& params"
  c) Strip conflicting #include / forward declarations from the namespace body
- Place each file's classes in a sanitized namespace to avoid name conflicts
- Skip classes that appear in 2+ files (duplicates: DynamicVacuumTerm×33 etc.)
- Skip boilerplate (DynamicVacuumTerm, QuantumCouplingTerm, DarkMatterHaloTerm)
"""

import re
import os
import glob
from collections import defaultdict

# Classes to unconditionally skip (generic boilerplate appearing 2+ times)
BOILERPLATE_CLASSES = {
    "DynamicVacuumTerm",
    "QuantumCouplingTerm",
    "DarkMatterHaloTerm",
}

def sanitize_namespace(filename):
    """
    Convert a source file path to a valid C++ namespace identifier.
    source47_wolfram.cpp -> src47w
    source48-5_wolfram.cpp -> src485w
    source82_wolfram_SMBH_ONLY_FINAL_BACKUP.cpp -> src82smb
    """
    base = os.path.splitext(filename)[0]  # remove .cpp
    # Normalize: lowercase, replace - and . with _
    base = base.replace('-', '_').replace('.', '_')
    m = re.match(r'source(\d+)_wolfram_?([a-z0-9_]*)', base, re.IGNORECASE)
    if m:
        num = m.group(1)
        suffix = re.sub(r'[^a-z0-9]', '', m.group(2).lower())[:4]
        return f"sw{num}{suffix}"
    # Fallback: sanitize the whole name
    safe = re.sub(r'[^a-z0-9_]', '', base.lower())
    return safe[:16] or 'swX'

def count_all_classes(files):
    """Count class appearances to find duplicates across all files."""
    counts = {}
    for fp in files:
        with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        for m in re.finditer(r'^class\s+(\w+Term)\s*:', content, re.MULTILINE):
            cn = m.group(1)
            counts[cn] = counts.get(cn, 0) + 1
    return counts

def extract_class_body(content, class_name):
    """
    Extract full class body text using brace counting.
    Returns (start_pos, end_pos, raw_text) or None.
    """
    pattern = rf'(?:^|\n)(class\s+{re.escape(class_name)}\b[^{{]*\{{)'
    m = re.search(pattern, content)
    if not m:
        return None
    start = m.start(1)
    brace_count = 0
    i = content.index('{', start)
    while i < len(content):
        ch = content[i]
        if ch == '{':
            brace_count += 1
        elif ch == '}':
            brace_count -= 1
            if brace_count == 0:
                end = i + 1
                if end < len(content) and content[end] == ';':
                    end += 1
                return content[start:end]
        i += 1
    return None

def adapt_class_body(class_body):
    """
    Minimally adapt the class body for inclusion in the bridge file:
    1. Change base class to ::PhysicsTerm
    2. Add ',const std::map<std::string,double>&' param to compute(double t)
    3. Remove any #include lines or inline forward declarations of PhysicsTerm
    """
    # 1. Base class
    adapted = re.sub(
        r':\s*public\s+PhysicsTerm\b',
        r': public ::PhysicsTerm',
        class_body
    )
    # 2. compute signature — handle both override cases
    adapted = re.sub(
        r'double\s+compute\s*\(\s*double\s+t\s*\)\s*const\s+override',
        'double compute(double t, const std::map<std::string, double>& /*params*/) const override',
        adapted
    )
    return adapted

def extract_return_str(body, method_name):
    """Extract the first quoted string returned by a no-arg method."""
    pattern = rf'{re.escape(method_name)}\s*\(\s*\).*?return\s+"([^"]+)"'
    m = re.search(pattern, body, re.DOTALL)
    return m.group(1) if m else None

def main():
    files = sorted(glob.glob('source*_wolfram*.cpp'))
    # Exclude infrastructure/bridge files that aren't PhysicsTerm definitions
    skip_files = {
        'source174_wolfram_bridge_embedded.cpp',
        'source175_uqff_wolfram_export.cpp',
        'source176_auto_full_uqff.cpp',
        'source177_wolfram_field_unity.cpp',
        'source200_wolfram.cpp',
    }
    files = [f for f in files if os.path.basename(f) not in skip_files]

    print(f"Processing {len(files)} source wolfram files...")

    # Count class appearances to find duplicates
    print("Scanning for duplicate class names...")
    class_counts = count_all_classes(files)
    dup_classes = {cn for cn, cnt in class_counts.items() if cnt > 1}
    print(f"  Duplicate class names (2+ files): {len(dup_classes)}")

    # Collect all unique classes and their adapted bodies
    bridge_items = []   # (ns_name, class_name, term_name, category, adapted_body, src_file)
    skipped_dup = 0
    skipped_nofind = 0

    for fp in files:
        with open(fp, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()

        fn = os.path.basename(fp)
        ns_name = sanitize_namespace(fn)

        for m in re.finditer(r'^class\s+(\w+Term)\s*:', content, re.MULTILINE):
            class_name = m.group(1)
            if class_name in dup_classes or class_name in BOILERPLATE_CLASSES:
                skipped_dup += 1
                continue

            body = extract_class_body(content, class_name)
            if not body:
                skipped_nofind += 1
                continue

            adapted = adapt_class_body(body)
            term_name = extract_return_str(body, 'getName') or class_name.replace('Term', '')
            category = extract_return_str(body, 'getCategory') or 'Wolfram-Source'

            bridge_items.append((ns_name, class_name, term_name, category, adapted, fn))

    # Deduplicate term names in registration to avoid overwriting
    seen_term_names = {}
    final_items = []
    for item in bridge_items:
        ns_name, class_name, term_name, category, adapted, fn = item
        if term_name in seen_term_names:
            # Make unique by appending namespace suffix
            term_name = f"{term_name}_{ns_name}"
        seen_term_names[term_name] = True
        final_items.append((ns_name, class_name, term_name, category, adapted, fn))

    print(f"  Unique bridge classes: {len(final_items)}")
    print(f"  Skipped (duplicates/boilerplate): {skipped_dup}")
    print(f"  Skipped (could not extract): {skipped_nofind}")

    # Group by namespace
    ns_groups = defaultdict(list)
    for item in final_items:
        ns_groups[item[0]].append(item)

    # Build the bridge file
    out = []
    out.append('// wolfram_sources_bridge.cpp')
    out.append('// AUTO-GENERATED by generate_wolfram_bridge.py')
    out.append(f'// Batch 21: {len(final_items)} unique system-specific PhysicsTerm bridge classes')
    out.append('// from source4-source85_wolfram.cpp companion files.')
    out.append('//')
    out.append('// Minimal modification: only compute() signature adapted for current API.')
    out.append('// DO NOT EDIT BY HAND - regenerate using generate_wolfram_bridge.py')
    out.append('')
    out.append('// Note: PhysicsTerm, CalculatorCore, std::map already defined before this include.')
    out.append('')
    out.append('#include <cmath>')
    out.append('#include <complex>')
    out.append('')

    for ns_name in sorted(ns_groups.keys()):
        items = ns_groups[ns_name]
        src_file = items[0][5]
        out.append(f'// ---- from {src_file} ----')
        out.append(f'namespace {ns_name} {{')
        out.append('')
        for ns_, cn, tn, cat, body, _ in items:
            out.append(body)
            out.append('')
        out.append(f'}} // namespace {ns_name}')
        out.append('')

    # Registration function
    out.append('// ====================================================================')
    out.append('// BATCH 21 REGISTRATION FUNCTION')
    out.append('// ====================================================================')
    out.append('')
    out.append('void registerAllWolframSourceTerms(CalculatorCore& core) {')
    out.append(f'    // Batch 21: {len(final_items)} unique system-specific Wolfram source companion terms')
    out.append('    // Generated from source4-source85_wolfram.cpp companion files')
    for ns_name, class_name, term_name, category, _, _ in final_items:
        safe_t = term_name.replace('"', '\\"')
        safe_c = category.replace('"', '\\"')
        out.append(f'    core.registerPhysicsTerm("{safe_t}", std::make_unique<{ns_name}::{class_name}>(), "{safe_c}");')
    out.append('}')
    out.append('')
    out.append(f'// Bridge terms total: {len(final_items)}')

    with open('wolfram_sources_bridge.cpp', 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))

    print(f'\nGenerated wolfram_sources_bridge.cpp ({len(final_items)} classes)')

    # Summary file
    with open('wolfram_extraction/batch21_summary.txt', 'w', encoding='utf-8') as f:
        f.write(f'Batch 21 Bridge Summary\n========================\n')
        f.write(f'Total bridge classes: {len(final_items)}\n')
        f.write(f'Skipped dups/boilerplate: {skipped_dup}\n')
        f.write(f'Skipped not-found: {skipped_nofind}\n\n')
        f.write('Classes per namespace:\n')
        for ns in sorted(ns_groups.keys()):
            f.write(f'  {ns}: {len(ns_groups[ns])}\n')
    print('Summary: wolfram_extraction/batch21_summary.txt')

if __name__ == '__main__':
    main()
