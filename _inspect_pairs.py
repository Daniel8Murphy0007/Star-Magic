from pathlib import Path

for n in range(1109, 1126):
    g = Path(f'whitepapers/PAPER_{n}.md')
    named_list = list(Path('whitepapers').glob(f'PAPER_{n}_*.md'))
    if g.exists() and named_list:
        nf = named_list[0]
        g_text = g.read_text(encoding='utf-8', errors='replace')
        n_text = nf.read_text(encoding='utf-8', errors='replace')
        print(f'PAPER_{n}: generic={len(g_text.splitlines())}L  named={len(n_text.splitlines())}L  named_file={nf.name}')
    elif g.exists():
        print(f'PAPER_{n}: generic EXISTS, no named file found')
    else:
        print(f'PAPER_{n}: generic MISSING')
