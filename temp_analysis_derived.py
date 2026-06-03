from pathlib import Path
import ast
import re

path = Path('Star-MagicProofEngine.py')
text = path.read_text(encoding='utf-8')

# Parse AST of the whole file
module = ast.parse(text)

# Gather callable names referenced in PROOF_DERIVATION_MODES
modes = {}
pattern = re.compile(r"'([^']+)':\s*\{[^}]*?'callable':\s*([^,\n]+)", re.S)
for name, expr in pattern.findall(text[text.index('self.PROOF_DERIVATION_MODES'):text.find('\ndef ', text.index('self.PROOF_DERIVATION_MODES'))]):
    modes[name] = expr.strip()

# Map method names to AST nodes
methods = {node.name: node for node in module.body if isinstance(node, ast.FunctionDef)}
# Also methods inside class
for node in module.body:
    if isinstance(node, ast.ClassDef):
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                methods[item.name] = item

# identify lambda-only entries
lambda_entries = {name: expr for name, expr in modes.items() if expr.startswith('lambda')}

# function to detect arithmetic or computation expressions
def has_computation(node):
    if isinstance(node, ast.BinOp):
        return True
    if isinstance(node, ast.UnaryOp):
        return True
    if isinstance(node, ast.Call):
        return True
    if isinstance(node, ast.Compare):
        return True
    if isinstance(node, ast.BoolOp):
        return True
    if isinstance(node, ast.If):
        return True
    if isinstance(node, ast.For):
        return True
    if isinstance(node, ast.While):
        return True
    if isinstance(node, ast.AugAssign):
        return True
    if isinstance(node, ast.Assign):
        return True
    if isinstance(node, ast.Return):
        return has_computation(node.value) if node.value else False
    if isinstance(node, ast.Dict):
        return any(has_computation(v) for v in node.values)
    if isinstance(node, ast.List):
        return any(has_computation(elt) for elt in node.elts)
    if isinstance(node, ast.Tuple):
        return any(has_computation(elt) for elt in node.elts)
    if isinstance(node, ast.IfExp):
        return True
    if isinstance(node, ast.Attribute):
        return True
    if isinstance(node, ast.Subscript):
        return True
    if isinstance(node, ast.Name):
        return False
    if isinstance(node, ast.Constant):
        return False
    if isinstance(node, ast.Lambda):
        return True
    return False

# Determine which methods are likely not explicit derivations
no_derivation = []
for name, expr in modes.items():
    if expr.startswith('lambda'):
        no_derivation.append((name, 'lambda placeholder', expr))
        continue
    if expr.startswith('self.'):
        method_name = expr.split('.', 1)[1]
        if '(' in method_name:
            method_name = method_name[:method_name.index('(')]
        node = methods.get(method_name)
        if node is None:
            no_derivation.append((name, 'missing method', expr))
            continue
        # inspect body for assignments/call/arith
        comp = False
        for stmt in node.body:
            if has_computation(stmt):
                comp = True
                break
        if not comp:
            no_derivation.append((name, 'method with no computation', method_name))

print('TOTAL MODES', len(modes))
print('LAMBDA / placeholders', len(lambda_entries))
for entry in no_derivation:
    print(entry)

# list method names with no direct computation but still referenced by registry
print('\nMETHODS REFERENCED BUT NO COMPUTATION BODY:')
for name, kind, expr in no_derivation:
    if kind == 'method with no computation':
        print(name, expr)
