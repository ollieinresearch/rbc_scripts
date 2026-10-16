"""Sphinx configuration for the documentation site.

Build locally (from the repository root)::

    pip install -e ".[docs]"
    sphinx-build -b html docs docs/_build/html      # then open docs/_build/html/index.html

or ``sphinx-autobuild docs docs/_build/html`` for a live preview that rebuilds on save.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "sim"))     # the simulation script, for its options
# The example gallery runs the plotting scripts on synthetic data, never on real results.
os.environ["RBC_DATA"] = str(ROOT / "docs" / "examples" / "data")
os.environ.setdefault("MPLBACKEND", "Agg")

import logging  # noqa: E402

import rbc  # noqa: E402

logging.getLogger("fontTools").setLevel(logging.ERROR)   # font-subsetting chatter in gallery output

# Edit to the repository's address once it is on GitHub.
GITHUB_URL = "https://github.com/ollieinresearch/rbc_scripts"

project = "Rayleigh–Bénard Convection Simulation Scripts"
author = "Oliver Brown"
copyright = "2020–2026, Oliver Brown and contributors"
release = version = rbc.__version__

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "numpydoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.mathjax",
    "sphinx_design",
    "sphinx_copybutton",
    "sphinxarg.ext",
    "sphinx_tippy",
    "sphinx_gallery.gen_gallery",
]

source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
exclude_patterns = ["_build", "examples/README.rst", "examples/data", "**.ipynb_checkpoints",
                    "gallery/auto_examples/*.ipynb", "gallery/auto_examples/*.py",
                    "gallery/auto_examples/*.zip", "gallery/auto_examples/*.md5",
                    "gallery/auto_examples/*.json"]
templates_path = ["_templates"]

# MyST (Markdown) features: $math$, ::: fences for directives, definition lists, {attrs}.
myst_enable_extensions = ["dollarmath", "amsmath", "colon_fence", "deflist", "attrs_inline",
                          "attrs_block", "substitution", "fieldlist"]
myst_heading_anchors = 3
myst_dmath_double_inline = True

# API reference from the docstrings. Modules needing MPI, Dedalus or VTK are imported
# with those packages replaced by stand-ins, so the docs build anywhere.
autosummary_generate = True
autosummary_imported_members = False
autodoc_mock_imports = ["dedalus", "mpi4py", "pyvista", "vtkmodules"]
autodoc_default_options = {"member-order": "bysource"}
autodoc_typehints = "none"
numpydoc_show_class_members = False
numpydoc_class_members_toctree = False
numpydoc_xref_param_type = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
    "scipy": ("https://docs.scipy.org/doc/scipy", None),
    "matplotlib": ("https://matplotlib.org/stable", None),
    "pandas": ("https://pandas.pydata.org/docs", None),
    "h5py": ("https://docs.h5py.org/en/stable", None),
    "mpi4py": ("https://mpi4py.readthedocs.io/en/stable", None),
    "dedalus": ("https://dedalus-project.readthedocs.io/en/latest", None),
}
intersphinx_timeout = 10

# Hover previews of cross-references, glossary terms and equations.
tippy_anchor_parent_selector = "article.bd-article"
tippy_enable_mathjax = True
tippy_enable_wikitips = False
tippy_props = {"placement": "auto-start", "maxWidth": 500, "interactive": False, "delay": [200, 100]}

copybutton_prompt_text = r"\$ |>>> "
copybutton_prompt_is_regexp = True

sphinx_gallery_conf = {
    "examples_dirs": "examples",
    "gallery_dirs": "gallery/auto_examples",
    "filename_pattern": r"/plot_",
    "download_all_examples": False,
    "remove_config_comments": True,
    "show_signature": False,
    "reset_modules": ("matplotlib",),
    "abort_on_example_error": False,
}

html_theme = "pydata_sphinx_theme"
html_title = "RBC Simulation Scripts"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_show_sourcelink = False
html_theme_options = {
    "logo": {"text": "RBC Simulation Scripts"},
    "icon_links": [{"name": "GitHub", "url": GITHUB_URL, "icon": "fa-brands fa-github"}],
    "header_links_before_dropdown": 6,
    "navbar_align": "left",
    "show_toc_level": 2,
    "navigation_depth": 3,
    "secondary_sidebar_items": ["page-toc"],
    "footer_start": ["copyright"],
    "footer_end": ["sphinx-version", "theme-version"],
    "show_prev_next": True,
}
html_context = {"default_mode": "auto"}
