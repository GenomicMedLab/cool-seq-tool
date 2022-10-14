"""Module for package and distribution."""
from setuptools import setup

exec(open("cool_seq_tool/version.py").read())
setup(version=__version__)  # noqa: F821
