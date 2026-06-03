from pathlib import Path
import re
import ast

p = Path('Star-MagicProofEngine.py')
text = p.read_text(encoding='utf-8')
start = text.index('self.PROOF_DERIVATION_MODES')
end = text.find('\ndef ', start)
block = text[start:end]
entries = re.findall(r"'([^']+)':\s*\{[^}]*?'callable':\s*([^,\n]+)", block, re.S)
methods = {}
module = ast.parse(text)
for node in module.body:
    if isinstance(node, ast.FunctionDef):
        methods[node.name] = node
    if isinstance(node, ast.ClassDef):
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                methods[item.name] = item


def analyze(node):
    if not isinstance(node, ast.FunctionDef):
        return 'unknown'
    has_assign_or_if = False
    for stmt in node.body:
        if isinstance(stmt, ast.Return):
            if isinstance(stmt.value, ast.Call):
                if isinstance(stmt.value.func, ast.Attribute) and isinstance(stmt.value.func.value, ast.Name) and stmt.value.func.value.id == 'self':
                    return 'wrapper'
            if isinstance(stmt.value, ast.Constant) and isinstance(stmt.value.value, (str, int, float)):
                return 'constant'
        if isinstance(stmt, (ast.Assign, ast.AnnAssign, ast.AugAssign, ast.If, ast.For, ast.While, ast.With)):
            has_assign_or_if = True
    if has_assign_or_if:
        return 'computed'
    return 'unknown'

incomplete = []
for name, expr in entries:
    expr = expr.strip()
    if expr.startswith('lambda'):
        incomplete.append((name, 'lambda placeholder', expr))
    elif expr.startswith('self.'):
        meth = expr.split('.', 1)[1]
        if '(' in meth:
            meth = meth[:meth.index('(')]
        node = methods.get(meth)
        if node is None:
            incomplete.append((name, 'missing method', expr))
        else:
            kind = analyze(node)
            if kind in ('constant', 'wrapper'):
                incomplete.append((name, kind, meth))
print('INCOMPLETE', len(incomplete))
for item in incomplete:
    print(item)
