"""Fix overflow-prone math.exp patterns in CP4 new classes."""
with open('CondensedPhysics4.py', 'r', encoding='utf-8') as f:
    content = f.read()

replacements = [
    (
        'growth = math.exp(self.KAPPA * t + self.SSQ * t / 26.0)',
        '_exp_arg = self.KAPPA * t + self.SSQ * t / 26.0; growth = math.exp(_exp_arg) if _exp_arg < 709 else float("inf")',
    ),
    (
        'E0 * math.exp(self.KAPPA * t + self.SSQ * t / 26.0) * s26',
        'E0 * (math.exp(_ea) if (_ea := self.KAPPA * t + self.SSQ * t / 26.0) < 709 else float("inf")) * s26',
    ),
    (
        'rho0 * s26 * math.exp(self.KAPPA * t + self.SSQ * t / 26.0)',
        'rho0 * s26 * (math.exp(_ea) if (_ea := self.KAPPA * t + self.SSQ * t / 26.0) < 709 else float("inf"))',
    ),
    (
        'self.RHO_VAC_SCM * s26 * math.exp(self.KAPPA * t + self.SSQ * t / 26.0)',
        'self.RHO_VAC_SCM * s26 * (math.exp(_ea) if (_ea := self.KAPPA * t + self.SSQ * t / 26.0) < 709 else float("inf"))',
    ),
]

total = 0
for old, new in replacements:
    c = content.count(old)
    content = content.replace(old, new)
    total += c
    print(f'Replaced {c}: {old[:60]}...')

with open('CondensedPhysics4.py', 'w', encoding='utf-8') as f:
    f.write(content)

print(f'Total replacements: {total}')
