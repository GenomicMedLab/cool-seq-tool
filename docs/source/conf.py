# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'cool-seq-tool'
copyright = "2023, Wagner Lab, Nationwide Children's Hospital"
author = "Wagner Lab, Nationwide Children's Hospital"
html_title = "cool-seq-tool"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.linkcode",
    "sphinx_copybutton",
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']

# -- autodoc things ----------------------------------------------------------
import os  # noqa: E402
import sys  # noqa: E402
sys.path.insert(0, os.path.abspath("../../cool_seq_tool"))

autodoc_preserve_defaults = True

# -- get version -------------------------------------------------------------
# from cool_seq_tool.version import __version__  # noqa: E402
# version = __version__
# release = version

# -- linkcode ----------------------------------------------------------------
def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None
    filename = info["module"].replace(".", "/")
    return f"https://github.com/GenomicMedLab/cool-seq-tool/blob/main/{filename}.py"

# -- code block style --------------------------------------------------------
pygments_style = "default"
pygements_dark_style = "monokai"
