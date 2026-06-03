from pathlib import Path
import re
import ast

p = Path('Star-MagicProofEngine.py')
text = p.read_text(encoding='utf-8')
start = text.index('self.PROOF_DERIVATION_MODES')
end = text.find('\ndef ', start)
block = text[start:end]
entries = re.findall(r"'([^']+)':\s*\{[^}]*?'callable':\s*([^,\n]+)", block, re.S)

module = ast.parse(text)
methods = {}
for node in module.body:
    if isinstance(node, ast.FunctionDef):
        methods[node.name] = node
    if isinstance(node, ast.ClassDef):
        for item in node.body:
            if isinstance(item, ast.FunctionDef):
                methods[item.name] = item

class Finder(ast.NodeVisitor):
    def __init__(self):
        self.has_compute = False
        self.has_assign = False
        self.has_call = False
        self.return_wrapper = False
        self.only_constant_return = False

    def visit_Assign(self, node):
        self.has_assign = True
        self.generic_visit(node)

    def visit_AnnAssign(self, node):
        self.has_assign = True
        self.generic_visit(node)

    def visit_AugAssign(self, node):
        self.has_assign = True
        self.generic_visit(node)

    def visit_If(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_For(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_While(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_With(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_Call(self, node):
        self.has_call = True
        self.generic_visit(node)

    def visit_BinOp(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_UnaryOp(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_BoolOp(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_Compare(self, node):
        self.has_compute = True
        self.generic_visit(node)

    def visit_Return(self, node):
        if isinstance(node.value, ast.Call):
            if isinstance(node.value.func, ast.Attribute) and isinstance(node.value.func.value, ast.Name) and node.value.func.value.id == 'self':
                self.return_wrapper = True
        if isinstance(node.value, ast.Constant):
            self.only_constant_return = True
        if isinstance(node.value, ast.Dict):
            self.only_constant_return = True
        self.generic_visit(node)


def analyze(node):
    finder = Finder()
    finder.visit(node)
    if finder.return_wrapper and not finder.has_assign and not finder.has_compute:
        return 'wrapper'
    if finder.only_constant_return and not finder.has_assign and not finder.has_compute and not finder.has_call:
        return 'constant_return'
    if not finder.has_assign and not finder.has_compute and not finder.has_call and finder.return_wrapper:
        return 'wrapper'
    return 'derived'

incomplete=[]
for name, expr in entries:
    expr = expr.strip()
    if expr.startswith('lambda'):
        incomplete.append((name, 'lambda placeholder', expr))
        continue
    if expr.startswith('self.'):
        meth = expr.split('.', 1)[1]
        if '(' in meth:
            meth = meth[:meth.index('(')]
        node = methods.get(meth)
        if node is None:
            incomplete.append((name, 'missing method', expr))
            continue
        kind = analyze(node)
        if kind != 'derived':
            incomplete.append((name, kind, meth))

print('INCOMPLETE', len(incomplete))
for item in incomplete:
    print(item)
