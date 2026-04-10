"""Sphinx configuration for cht_cyclones."""

import os
import sys

# Make the package importable for autodoc
sys.path.insert(0, os.path.abspath(".."))

project = "cht_cyclones"
copyright = "2024, Deltares"
author = "Maarten van Ormondt, Kees Nederhoff"
release = "1.0.3"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "geopandas": ("https://geopandas.org/en/stable/", None),
}

napoleon_numpy_docstring = True
napoleon_google_docstring = False
autodoc_member_order = "bysource"
autodoc_typehints = "description"
