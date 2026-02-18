import re

# Template for run_tests method
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

# Classes missing run_tests (class_name, line_number)
missing_classes = [
    ('UQFFLogger', 49),
    ('JSONRPCServer', 108),
    ('VacuumFluctuationCalculator', 3755),
    ('MagneticCalculator', 3801),
    ('QuantumCalculator', 3826),
    ('ResonanceCalculator', 3854),
    ('TriadicCalculator', 3905),
    ('BuoyancyCalculator', 3934),
    ('CosmologicalCalculator', 4109),
    ('SuperconductiveCalculator', 4141),
    ('MUGECalculator', 4181),
    ('EquationResult', 4493),
    ('FloydSweetVacuumCalculator', 4506),
    ('CosmicEgg26DCalculator', 4608),
    ('HeisenbergVacuumCalculator', 4697),
    ('NegativeTimeCalculator', 4793),
    ('QuantumWaveMixin', 4903),
    ('NegativeTimeModel', 75631),
    ('AetherVacuumEnergyModel', 75756),
    ('CosmicEggModel', 75875),
    ('SgrAStarGravityModel', 76012),
    ('RetrocausalModel', 76139),
    ('TRZModel', 76262),
    ('VoidOscillationModel', 76381),
    ('TimeVaryingVacuumModel', 76500),
    ('TriadicGravityCalculator', 77081),
    ('UnifiedFieldSolver', 77163),
    ('MagnetarMUGECalculator', 77244),
    ('SgrAStarCalculator', 77386),
    ('Phase2Calculator', 77430),
    ('ConsciousnessCloud', 77516),
    ('StarMagicEnergyStructure', 77576),
    ('StarMagicBlackHoleInteraction', 77635),
    ('StarMagicVacuumEnergy', 77677),
    ('GravitationalCalculator', 77728),
    ('CoAnQiCalculator', 77834),
    ('EquationFamily', 77963),
    ('ReferenceSystem', 78034),
    ('ReferenceSystemLibrary', 78075),
]

with open('CondensedPhysics.py', 'r', encoding='utf-8') as f:
    content = f.read()

# Find all class positions
class_pattern = r'^class (\w+).*?:'
class_positions = [(m.group(1), m.start()) for m in re.finditer(class_pattern, content, re.MULTILINE)]

# Build a map of class name -> next class start position
class_ends = {}
for i, (name, start) in enumerate(class_positions):
    if i + 1 < len(class_positions):
        class_ends[name] = class_positions[i + 1][1]
    else:
        class_ends[name] = len(content)

# Insert run_tests methods (work backwards to preserve positions)
insertions = []
for class_name, _ in missing_classes:
    if class_name in class_ends:
        end_pos = class_ends[class_name]
        # Find insertion point (before next class definition)
        # We need to find the last newline before the next class
        insert_pos = content.rfind('\n\n', 0, end_pos)
        if insert_pos == -1:
            insert_pos = end_pos - 1
        method = RUN_TESTS_TEMPLATE.format(class_name=class_name)
        insertions.append((insert_pos, method, class_name))

# Sort by position descending to insert from end to start
insertions.sort(key=lambda x: x[0], reverse=True)

added = 0
for pos, method, class_name in insertions:
    content = content[:pos] + method + content[pos:]
    added += 1
    print(f'Added run_tests() to {class_name}')

with open('CondensedPhysics.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f'\nTotal: {added} run_tests() methods added')
