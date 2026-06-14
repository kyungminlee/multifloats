# Sphinx configuration for the multifloats documentation.
#
# Build (after `doxygen Doxyfile` has produced the XML):
#     doc/.venv/bin/sphinx-build -b html doc doc/_build/html
# or use the Makefile / build.sh in this directory.

import os.path

# -- Project information -----------------------------------------------------
project = "multifloats"
author = "Kyungmin Lee"
copyright = "2026, Kyungmin Lee"

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",        # Markdown (MyST) authoring
    "breathe",            # Doxygen XML -> Sphinx
    "sphinx_copybutton",  # copy button on code blocks
    "sphinx_design",      # grids / cards / tabs
]

# Files and directories ignored by the builder.
exclude_patterns = [
    "_build",
    ".venv",
    "_doxygen",
    "Doxyfile",
    "requirements.txt",
    "README.md",          # build instructions, not a doc page
    "developer",          # developer-only notes, not part of the user docs
]

# -- MyST (Markdown) ---------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",   # ::: fenced directives
    "deflist",       # definition lists
    "dollarmath",    # $...$ / $$...$$ math
    "fieldlist",
    "substitution",
    "tasklist",
]
myst_heading_anchors = 3  # auto cross-reference anchors for h1..h3

# -- Breathe (C / C++ API extraction) ----------------------------------------
breathe_projects = {
    "multifloats": os.path.join(os.path.dirname(__file__), "_doxygen", "xml"),
}
breathe_default_project = "multifloats"
breathe_default_members = ("members",)
breathe_domain_by_extension = {"h": "cpp"}  # treat .h as C++

# -- HTML output -------------------------------------------------------------
html_theme = "furo"
html_title = "multifloats"
html_static_path = ["_static"]
