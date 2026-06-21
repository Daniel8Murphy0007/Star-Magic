"""Sphinx configuration for UQFF Star-Magic documentation."""
import os
import sys

# Add repo root so autodoc can import uqff_pure_calculator
sys.path.insert(0, os.path.abspath(".."))

# -- Project information -----------------------------------------------------
project = "UQFF Star-Magic"
copyright = "2025-2026, Daniel T. Murphy / Star-Magic Research Program"
author = "Daniel T. Murphy"
release = "5.27.0"
version = "5.27"

# -- General configuration ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# Suppress warnings for the calculator's no-docstring style (Rule 3 compliance)
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": False,
}
autodoc_mock_imports = []
autosummary_generate = True
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# -- HTML output -------------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_title = "UQFF Star-Magic v5.27 Documentation"
html_short_title = "UQFF"

html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "includehidden": True,
    "titles_only": False,
    "prev_next_buttons_location": "both",
}

html_context = {
    "display_github": True,
    "github_user": "Daniel8Murphy0007",
    "github_repo": "Star-Magic",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# -- Cross-references --------------------------------------------------------
intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# -- LaTeX output (for PDF builds) ------------------------------------------
latex_elements = {
    "papersize": "letterpaper",
    "pointsize": "10pt",
    "preamble": "",
    "figure_align": "htbp",
}

latex_documents = [
    ("index", "UQFF.tex", "UQFF Star-Magic Documentation",
     "Daniel T. Murphy", "manual"),
]
