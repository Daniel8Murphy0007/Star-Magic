import re

# Simpler approach: find each class and insert before the next class definition
RUN_TESTS_TEMPLATE = """
    @staticmethod
    def run_tests() -> dict:
        \"\"\"Run self-diagnostic tests for {class_name}.\"\"\"
        tests = []
        try:
            instance = {class_name}()
            tests.append({{'name': 'Instantiation', 'passed': instance is not None}})
        except Exception as e:
            tests.append({{'name': 'Instantiation', 'passed': False, 'error': str(e)}})
        return {{'class': '{class_name}', 'tests': tests, 'passed': sum(1 for t in tests if t['passed']), 'total': len(tests)}}
"""

missing_classes = [
    'UQFFLogger', 'JSONRPCServer', 'VacuumFluctuationCalculator', 'MagneticCalculator',
    'QuantumCalculator', 'ResonanceCalculator', 'TriadicCalculator', 'BuoyancyCalculator',
    'CosmologicalCalculator', 'SuperconductiveCalculator', 'MUGECalculator', 'EquationResult',
    'FloydSweetVacuumCalculator', 'CosmicEgg26DCalculator', 'HeisenbergVacuumCalculator',
    'NegativeTimeCalculator', 'QuantumWaveMixin', 'NegativeTimeModel', 'AetherVacuumEnergyModel',
    'CosmicEggModel', 'SgrAStarGravityModel', 'RetrocausalModel', 'TRZModel', 'VoidOscillationModel',
    'TimeVaryingVacuumModel', 'TriadicGravityCalculator', 'UnifiedFieldSolver', 'MagnetarMUGECalculator',
    'SgrAStarCalculator', 'Phase2Calculator', 'ConsciousnessCloud', 'StarMagicEnergyStructure',
    'StarMagicBlackHoleInteraction', 'StarMagicVacuumEnergy', 'GravitationalCalculator',
    'CoAnQiCalculator', 'EquationFamily', 'ReferenceSystem', 'ReferenceSystemLibrary',
]

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    lines = f.readlines()

# Find class start lines (0-indexed)
class_starts = {}
for i, line in enumerate(lines):
    match = re.match(r'^class (\w+)', line)
    if match:
        class_starts[match.group(1)] = i

# Find where to insert for each missing class
insertions = []  # (line_before, method_text, class_name)
for class_name in missing_classes:
    if class_name not in class_starts:
        print(f'WARNING: {class_name} not found')
        continue
    
    class_start = class_starts[class_name]
    
    # Find the next class definition or end of file
    next_class_line = None
    for name, start in class_starts.items():
        if start > class_start:
            if next_class_line is None or start < next_class_line:
                next_class_line = start
    
    if next_class_line is None:
        next_class_line = len(lines)
    
    # Find the right insertion point - the last "    def " in this class block
    insert_line = None
    for i in range(next_class_line - 1, class_start, -1):
        # Look for a line that's a method definition or the end of a method docstring/body
        if lines[i].strip() and not lines[i].startswith('#'):
            # Check if we're still inside the class (indented content)
            if lines[i].startswith('    ') or lines[i].strip() == '':
                insert_line = i
                break
    
    if insert_line is None:
        insert_line = next_class_line - 1
    
    # Find the last "return" statement in the last method of this class
    for i in range(next_class_line - 1, class_start, -1):
        if lines[i].strip().startswith('return ') or lines[i].strip().startswith("return {"):
            insert_line = i
            break
    
    method = RUN_TESTS_TEMPLATE.format(class_name=class_name)
    insertions.append((insert_line + 1, method, class_name))

# Sort by line number descending and insert
insertions.sort(key=lambda x: x[0], reverse=True)

added = 0
for line_num, method, class_name in insertions:
    lines.insert(line_num, method + '\n')
    added += 1
    print(f'Added run_tests() to {class_name} after line {line_num}')

with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
    f.writelines(lines)

print(f'\nTotal: {added} run_tests() methods added')
