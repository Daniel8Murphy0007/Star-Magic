from pathlib import Path

path = Path('Star-MagicProofEngine.py')
text = path.read_text(encoding='utf-8')

for name in ['_prove_refactored_umbilicus_mass_balance','_compute_quantum_variables_placeholder','_compute_heaviside_component_fraction']:
    idx = text.find(f'def {name}')
    if idx == -1:
        print(name, 'NOT FOUND')
        continue
    start = idx
    end = text.find('\ndef ', start + 1)
    if end == -1:
        end = len(text)
    snippet = text[start:end]
    print('===', name, '===')
    print(snippet)
    print('\n')
