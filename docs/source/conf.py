# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "Cool-Seq-Tool"
copyright = "2021-2024, Wagner Lab"
author = "Wagner Lab"
html_title = "Cool-Seq-Tool"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx_autodoc_typehints",
    "sphinx.ext.linkcode",
    "sphinx.ext.autosummary",
    "sphinx_copybutton",
    "sphinx_github_changelog",
]

templates_path = ["_templates"]
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "furo"
html_static_path = []
html_css_files = [
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/fontawesome.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/solid.min.css",
    "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/brands.min.css",
]
html_theme_options = {
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/genomicmedlab/cool-seq-tool",
            "html": "",
            "class": "fa-brands fa-solid fa-github",
        },
        {
            "name": "Wagner Lab",
            "url": "https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab",
            "html": "",
            "class": "fa-solid fa-house",
        },
    ],
}
# -- autodoc things ----------------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.abspath("../../cool_seq_tool"))
autodoc_preserve_defaults = True

# -- get version -------------------------------------------------------------
from cool_seq_tool.version import __version__

version = __version__
release = version


# -- linkcode ----------------------------------------------------------------
def linkcode_resolve(domain, info):
    if domain != "py":
        return None
    if not info["module"]:
        return None
    filename = info["module"].replace(".", "/")
    return f"https://github.com/genomicmedlab/cool-seq-tool/blob/main/{filename}.py"  # noqa: E501


# -- code block style --------------------------------------------------------
pygments_style = "default"
pygements_dark_style = "monokai"

# -- changelog ---------------------------------------------------------------
# import os
#
# sphinx_github_changelog_token = os.environ.get("GENOMICMEDLAB_GH_DOCS_API_TOKEN", "failed -- make sure to set API token")
