"""Module for package and distribution."""
from setuptools import setup

exec(open("uta_tools/version.py").read())
setup(version=__version__)  # noqa: F821
