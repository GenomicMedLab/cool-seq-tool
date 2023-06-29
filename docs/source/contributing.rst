Contributing to Cool Seq Tool
=============================

Bug reports and feature requests
--------------------------------

Bugs and new feature requests can be submitted to the `issue tracker on GitHub <https://github.com/GenomicMedLab/cool-seq-tool/issues>`_. See `this StackOverflow post <https://stackoverflow.com/help/minimal-reproducible-example>`_ for tips on how to craft a helpful bug report.

Development prerequisites
-------------------------
For a development install, we recommend using Pipenv. See the `pipenv docs <https://pipenv-fork.readthedocs.io/en/latest/#install-pipenv-today>`_ for direction on installing Pipenv in your environment.

Setup
-----
Clone the repository: ::

    git clone https://github.com/GenomicMedLab/cool-seq-tool
    cd cool-seq-tool

Then initialize the Pipenv environment: ::

    pipenv update
    pipenv install --dev
    pipenv shell

Alternatively, use a virtual environment and install all dependency groups: ::

    python3 -m venv venv
    source venv/bin/activate
    python3 -m pip install -e ".[tests,dev,docs]"

We use `pre-commit <https://pre-commit.com/#usage>`_ to run conformance tests before commits. This provides checks for:

* Code style
* Added large files
* Private keys

Before your first commit, run: ::

    pre-commit install

When running the web server, enable hot-reloading on new code changes: ::

    uvicorn cool_seq_tool.api:app --reload

Style
-----

Code style is managed by `flake8 <https://github.com/PyCQA/flake8>`_ and should be checked via pre-commit hook before commits. Final QC is applied with GitHub Actions to every pull request.

Tests
-----

Tests are executed with `pytest <https://docs.pytest.org/en/7.1.x/getting-started.html>`_: ::

    pytest

Documentation
-------------

The documentation is built with Sphinx, which is included as part of the developer dependencies. To build a local copy, ensure that the Gene Normalizer is installed in your current Python environment, then navigate to the `docs/` subdirectory and use `make` to build the HTML version: ::

    pipenv shell
    cd docs
    make html

See the `Sphinx documentation <https://www.sphinx-doc.org/en/master/>`_ for more information.
