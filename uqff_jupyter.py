"""UQFF Star-Magic Jupyter / IPython integration.

Provides rich HTML rendering of UQFF closure results in Jupyter notebooks,
plus an `%uqff` line magic for quick inline lookups.

Install + activate:
    pip install 'uqff[jupyter]'

In a Jupyter notebook:
    import uqff_jupyter
    uqff_jupyter.enable()                    # activate rich display

    %load_ext uqff_jupyter                   # register the %uqff magic
    %uqff predict hubble_tension             # inline closure fetch
    %uqff search holmlid                     # inline multi-namespace search

After enable() is called, any returned dict from a `calculate_*` surface
renders as a formatted HTML table instead of raw `repr()`.
"""
import html as _html
import json
from typing import Any

import uqff_pure_calculator as u

try:
    from uqff_cli import (
        _try_paradox, _try_lenr_full, _try_nuclear, _try_bucket_observable,
        _all_paradox_keys, _all_millennium_keys, _all_lenr_keys,
        _all_nuclear_keys, _all_bucket_observables,
    )
except ImportError:
    _try_paradox = lambda n: u.calculate_paradox({"paradox": n}).get("value")
    _try_lenr_full = lambda n: None
    _try_nuclear = lambda n: None
    _try_bucket_observable = lambda n: None
    def _all_paradox_keys(): return sorted(u.PARADOX_TO_CLOSURE.keys())
    def _all_millennium_keys(): return []
    def _all_lenr_keys(): return []
    def _all_nuclear_keys(): return []
    def _all_bucket_observables(): return {}


_VERSION = "5.59.0"


# ---------------------------------------------------------------------------
# HTML rendering for closure results
# ---------------------------------------------------------------------------

_CSS = """
<style>
.uqff-card {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif;
  border: 1px solid #d0d7de; border-radius: 6px; padding: 12px;
  margin: 8px 0; background: #f6f8fa; color: #1f2328;
}
.uqff-card .uqff-title { font-weight: 600; font-size: 14px; color: #0969da; margin-bottom: 8px; }
.uqff-card .uqff-source { font-size: 11px; color: #57606a; text-transform: uppercase; letter-spacing: 0.05em; margin-bottom: 6px; }
.uqff-card table { border-collapse: collapse; width: 100%; font-size: 13px; }
.uqff-card td, .uqff-card th { padding: 4px 8px; border-bottom: 1px solid #d8dee4; text-align: left; }
.uqff-card th { background: #ddf4ff; color: #0550ae; font-weight: 600; }
.uqff-card td.num { text-align: right; font-variant-numeric: tabular-nums; }
.uqff-card .uqff-pill { display: inline-block; padding: 1px 6px; border-radius: 10px; font-size: 10px; font-weight: 500; background: #dafbe1; color: #1a7f37; margin-left: 4px; }
.uqff-card .uqff-residual-bad { color: #cf222e; }
.uqff-card .uqff-residual-good { color: #1a7f37; }
</style>
""".strip()


def _fmt_num(v: Any) -> str:
    if isinstance(v, (int, float)) and not isinstance(v, bool):
        if abs(v) > 0 and (abs(v) < 1e-3 or abs(v) > 1e6):
            return f"{v:.6e}"
        return f"{v:.6g}"
    return str(v)[:120]


def _render_closure_dict(d: dict, title: str = "", source: str = "") -> str:
    rows = []
    for k, v in d.items():
        if isinstance(v, (dict, list)):
            preview = json.dumps(v, default=str)[:80]
            preview = _html.escape(preview)
            val_html = f'<code style="font-size:11px;color:#57606a">{preview}{"..." if len(json.dumps(v, default=str)) > 80 else ""}</code>'
        else:
            val_html = _html.escape(_fmt_num(v))
            if 'residual' in k.lower():
                try:
                    pct = float(v)
                    css_class = 'uqff-residual-good' if pct < 1.0 else 'uqff-residual-bad'
                    val_html = f'<span class="{css_class}">{pct:.4f}%</span>'
                except (TypeError, ValueError):
                    pass
        rows.append(f'<tr><td>{_html.escape(str(k))}</td><td class="num">{val_html}</td></tr>')

    rows_html = "\n".join(rows) if rows else "<tr><td colspan=2><em>(empty)</em></td></tr>"
    src_html = f'<div class="uqff-source">{_html.escape(source)}</div>' if source else ""
    title_html = f'<div class="uqff-title">{_html.escape(title)}</div>' if title else ""
    return (
        f'{_CSS}\n<div class="uqff-card">{title_html}{src_html}'
        f'<table>{rows_html}</table></div>'
    )


def _render_observables(obs_list: list, title: str = "") -> str:
    cols = ["observable", "uqff_derived", "anchor", "residual_pct"]
    header = '<tr>' + ''.join(f'<th>{c}</th>' for c in cols) + '</tr>'
    rows = []
    for o in obs_list[:50]:
        if not isinstance(o, dict):
            continue
        row_cells = []
        for c in cols:
            v = o.get(c, '')
            cell_class = ' class="num"' if c != 'observable' else ''
            txt = _html.escape(_fmt_num(v))
            if c == 'residual_pct':
                try:
                    pct = float(v)
                    css = 'uqff-residual-good' if pct < 1.0 else 'uqff-residual-bad'
                    txt = f'<span class="{css}">{pct:.4f}%</span>'
                except (TypeError, ValueError):
                    pass
            row_cells.append(f'<td{cell_class}>{txt}</td>')
        rows.append('<tr>' + ''.join(row_cells) + '</tr>')

    note = f'<div class="uqff-source">{len(obs_list)} observables; showing first 50</div>' if len(obs_list) > 50 else ''
    title_html = f'<div class="uqff-title">{_html.escape(title)}</div>' if title else ""
    return (
        f'{_CSS}\n<div class="uqff-card">{title_html}{note}'
        f'<table>{header}{"".join(rows)}</table></div>'
    )


def render_uqff_result(value: Any, source: str = "", name: str = "") -> str:
    """Convert any UQFF result into rich HTML for Jupyter display."""
    title = f"UQFF closure: {name}" if name else "UQFF closure"
    if isinstance(value, dict):
        if "observables" in value and isinstance(value["observables"], list):
            return _render_observables(value["observables"], title=title)
        return _render_closure_dict(value, title=title, source=source)
    if isinstance(value, list):
        if value and isinstance(value[0], dict) and 'observable' in value[0]:
            return _render_observables(value, title=title)
        items = "".join(f'<li>{_html.escape(_fmt_num(v))}</li>' for v in value[:50])
        return f'{_CSS}\n<div class="uqff-card"><div class="uqff-title">{_html.escape(title)}</div><ul>{items}</ul></div>'
    return f'{_CSS}\n<div class="uqff-card"><div class="uqff-title">{_html.escape(title)}</div><code>{_html.escape(_fmt_num(value))}</code></div>'


# ---------------------------------------------------------------------------
# enable() — install display hook
# ---------------------------------------------------------------------------

def enable():
    """Activate rich display rendering for UQFF closure results in this kernel."""
    try:
        from IPython import get_ipython
        from IPython.display import HTML, display
    except ImportError as e:
        raise SystemExit(
            "IPython not installed. Install with: pip install 'uqff[jupyter]'"
        ) from e

    shell = get_ipython()
    if shell is None:
        print("uqff_jupyter.enable(): not in an IPython kernel; no-op.")
        return

    html_formatter = shell.display_formatter.formatters['text/html']

    def _uqff_dict_formatter(obj):
        if not isinstance(obj, dict):
            return None
        if 'value' in obj and len(obj) == 1:
            return render_uqff_result(obj['value'])
        return None

    html_formatter.for_type(dict, _uqff_dict_formatter)
    print(f"UQFF Jupyter rich display ENABLED (v{_VERSION}).")
    print("  Any calculate_*({}) result will now render as a styled HTML table.")
    print("  Use `%load_ext uqff_jupyter` to enable the `%uqff` line magic.")


# ---------------------------------------------------------------------------
# %uqff line magic
# ---------------------------------------------------------------------------

def load_ipython_extension(ipython):
    """Called by `%load_ext uqff_jupyter`."""
    from IPython.core.magic import register_line_magic
    from IPython.display import HTML, display

    @register_line_magic
    def uqff(line):
        """%uqff <subcommand> [args]

        Subcommands:
          predict <name>     fetch a closure value
          search <substr>    multi-namespace search
          list [--filter X]  list dispatch keys
          status             production status summary
          version            version + metrics
        """
        parts = line.strip().split()
        if not parts:
            display(HTML(render_uqff_result({"usage": "%uqff <subcommand> [args]"}, name="uqff magic")))
            return

        cmd = parts[0].lower()
        args = parts[1:]

        if cmd == "predict" and args:
            name = args[0].lower()
            for source_name, fn in [
                ("PARADOX_TO_CLOSURE", _try_paradox),
                ("calculate_lenr_full", _try_lenr_full),
                ("calculate_nuclear_magic", _try_nuclear),
                ("bucket_observables", _try_bucket_observable),
            ]:
                value = fn(name)
                if value is not None:
                    display(HTML(render_uqff_result(value, source=source_name, name=name)))
                    return
            display(HTML(render_uqff_result({"error": f"closure '{name}' not found"})))
        elif cmd == "search" and args:
            needle = args[0].lower()
            results = {}
            for source, keys in [
                ("PARADOX_TO_CLOSURE", _all_paradox_keys()),
                ("calculate_lenr_full", _all_lenr_keys()),
                ("calculate_nuclear_magic", _all_nuclear_keys()),
            ]:
                hits = [k for k in keys if needle in k.lower()]
                if hits:
                    results[source] = hits[:20]
            for surf, names in _all_bucket_observables().items():
                hits = [n for n in names if needle in n.lower()]
                if hits:
                    results[surf] = hits[:20]
            display(HTML(render_uqff_result(results, name=f"search '{needle}'")))
        elif cmd == "status":
            try:
                summary = u.calculate_status_report({})["value"]["summary"]
                display(HTML(render_uqff_result(summary, name="status")))
            except Exception as e:
                display(HTML(render_uqff_result({"error": str(e)})))
        elif cmd == "version":
            display(HTML(render_uqff_result({
                "uqff_version": _VERSION,
                "closures": len(u.PARADOX_TO_CLOSURE),
            }, name="version")))
        elif cmd == "list":
            keys = _all_paradox_keys()
            flt = None
            if "--filter" in args:
                i = args.index("--filter")
                if i + 1 < len(args):
                    flt = args[i+1].lower()
            if flt:
                keys = [k for k in keys if flt in k.lower()]
            display(HTML(render_uqff_result({"count": len(keys), "first_50": keys[:50]}, name="list")))
        else:
            display(HTML(render_uqff_result(
                {"usage": "%uqff <predict|search|list|status|version> [args]"})))

    print("uqff_jupyter loaded. Use %uqff <subcommand> for quick lookups.")
