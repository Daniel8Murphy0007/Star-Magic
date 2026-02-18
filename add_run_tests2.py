import re

RUN_TESTS_TEMPLATE = '''
    @staticmethod
    def run_tests() -> dict:
        """Run self-diagnostic tests for {class_name}."""
        tests = []
        try:
            instance = {class_name}()
            tests.append({{'name': 'Instantiation', 'passed': instance is not None}})
        except Exception as e:
            tests.append({{'name': 'Instantiation', 'passed': False, 'error': str(e)}})
        return {{'class': '{class_name}', 'tests': tests, 'passed': sum(1 for t in tests if t['passed']), 'total': len(tests)}}
'''

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
    content = f.read()

# Find class boundaries properly
class_pattern = r'^class (\w+)[^\n]*:'

class_positions = []
for m in re.finditer(class_pattern, content, re.MULTILINE):
    class_positions.append((m.group(1), m.start(), m.end()))

# Build insertion points - find end of each class (before next class or EOF)
insertions = []
for i, (name, start, _) in enumerate(class_positions):
    if name not in missing_classes:
        continue
        
    # Find where this class ends (start of next class at column 0)
    if i + 1 < len(class_positions):
        next_class_start = class_positions[i + 1][1]
    else:
        next_class_start = len(content)
    
    # Find the last non-blank line before next class
    class_block = content[start:next_class_start]
    
    # Find the last method definition in this class
    method_matches = list(re.finditer(r'\n(    def \w+\([^)]*\).*?)(?=\n    def |\nclass |\n# [=]+|\Z)', class_block, re.DOTALL))
    
    if method_matches:
        last_method_end = start + method_matches[-1].end()
        # Insert after last method
        insert_pos = last_method_end
    else:
        # No methods found, insert before next class
        insert_pos = next_class_start - 1
    
    method = RUN_TESTS_TEMPLATE.format(class_name=name)
    insertions.append((insert_pos, method, name))

# Sort descending by position
insertions.sort(key=lambda x: x[0], reverse=True)

added = 0
for pos, method, class_name in insertions:
    content = content[:pos] + '\n' + method + content[pos:]
    added += 1
    print(f'Added run_tests() to {class_name}')

with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f'\nTotal: {added} run_tests() methods added')
