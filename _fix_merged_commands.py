"""Fix merged LaTeX commands: \hbaromega -> \hbar\omega, \partialt -> \partial t, etc."""
import os, re

# Merged command patterns: merged_form -> correct_form
# Must be ordered longest-first to avoid partial matches
FIXES = [
    # \hbar + Greek
    (r'\hbarvarepsilon', r'\hbar\varepsilon'),
    (r'\hbarepsilon', r'\hbar\epsilon'),
    (r'\hbaromega', r'\hbar\omega'),
    (r'\hbarOmega', r'\hbar\Omega'),
    (r'\hbarsigma', r'\hbar\sigma'),
    (r'\hbarlambda', r'\hbar\lambda'),
    (r'\hbarDelta', r'\hbar\Delta'),
    (r'\hbaralpha', r'\hbar\alpha'),
    (r'\hbargamma', r'\hbar\gamma'),
    (r'\hbarGamma', r'\hbar\Gamma'),
    (r'\hbarTheta', r'\hbar\Theta'),
    (r'\hbartheta', r'\hbar\theta'),
    (r'\hbarkappa', r'\hbar\kappa'),
    (r'\hbarbeta', r'\hbar\beta'),
    (r'\hbareta', r'\hbar\eta'),
    (r'\hbarphi', r'\hbar\phi'),
    (r'\hbarPhi', r'\hbar\Phi'),
    (r'\hbarpsi', r'\hbar\psi'),
    (r'\hbarPsi', r'\hbar\Psi'),
    (r'\hbarzeta', r'\hbar\zeta'),
    (r'\hbarrho', r'\hbar\rho'),
    (r'\hbartau', r'\hbar\tau'),
    (r'\hbarchi', r'\hbar\chi'),
    (r'\hbarpi', r'\hbar\pi'),
    (r'\hbarPi', r'\hbar\Pi'),
    (r'\hbarmu', r'\hbar\mu'),
    (r'\hbarnu', r'\hbar\nu'),
    (r'\hbarxi', r'\hbar\xi'),
    # \hbar + variable
    (r'\hbarc', r'\hbar c'),
    (r'\hbark', r'\hbar k'),
    # \nabla + commands
    (r'\nablacdot', r'\nabla\cdot'),
    (r'\nablatimes', r'\nabla\times'),
    (r'\nablaTimes', r'\nabla\times'),
    (r'\nablaomega', r'\nabla\omega'),
    (r'\nablasigma', r'\nabla\sigma'),
    (r'\nablarho', r'\nabla\rho'),
    (r'\nablaphi', r'\nabla\phi'),
    (r'\nablaPhi', r'\nabla\Phi'),
    (r'\nablapsi', r'\nabla\psi'),
    (r'\nablamu', r'\nabla\mu'),
    # \partial + Greek (longer first)
    (r'\partialvarepsilon', r'\partial\varepsilon'),
    (r'\partialomega', r'\partial\omega'),
    (r'\partialsigma', r'\partial\sigma'),
    (r'\partiallambda', r'\partial\lambda'),
    (r'\partialtheta', r'\partial\theta'),
    (r'\partialDelta', r'\partial\Delta'),
    (r'\partialalpha', r'\partial\alpha'),
    (r'\partialgamma', r'\partial\gamma'),
    (r'\partialkappa', r'\partial\kappa'),
    (r'\partialbeta', r'\partial\beta'),
    (r'\partialeta', r'\partial\eta'),
    (r'\partialphi', r'\partial\phi'),
    (r'\partialrho', r'\partial\rho'),
    (r'\partialtau', r'\partial\tau'),
    (r'\partialmu', r'\partial\mu'),
    (r'\partialnu', r'\partial\nu'),
    # \partial + variable  
    (r'\partialt', r'\partial t'),
    (r'\partialr', r'\partial r'),
    (r'\partialx', r'\partial x'),
    (r'\partialy', r'\partial y'),
    (r'\partialz', r'\partial z'),
    # \sqrt merges
    (r'\sqrtG', r'\sqrt{G}'),
    (r'\sqrtpi', r'\sqrt{\pi}'),
    # \frac merges  
    (r'\fracdelta', r'\frac{\delta}'),
    (r'\fracpartial', r'\frac{\partial}'),
    # \vec merges
    (r'\vecnabla', r'\vec{\nabla}'),
    # \dot merges
    (r'\dotomega', r'\dot{\omega}'),
    (r'\dotphi', r'\dot{\phi}'),
    (r'\dottheta', r'\dot{\theta}'),
    (r'\dotr', r'\dot{r}'),
    (r'\dota', r'\dot{a}'),
]

fixed_count = 0
fix_details = {}

for fn in sorted(os.listdir('whitepapers')):
    if not fn.endswith('.md'):
        continue
    path = os.path.join('whitepapers', fn)
    with open(path, encoding='utf-8') as f:
        content = f.read()
    
    original = content
    for merged, correct in FIXES:
        if merged in content:
            content = content.replace(merged, correct)
            fix_details[merged] = fix_details.get(merged, 0) + 1
    
    if content != original:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        fixed_count += 1

print(f"Fixed merged commands in {fixed_count} papers")
if fix_details:
    print("\nDetails:")
    for k, v in sorted(fix_details.items(), key=lambda x: -x[1]):
        print(f"  {k} -> {v} papers")
