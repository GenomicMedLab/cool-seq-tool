.. _installation:

Installation
============

Prerequisites
-------------

* Python 3.8+
* A UNIX-like environment (e.g. MacOS, WSL, Ubuntu)
* A recent version of PostgreSQL (ideally at least 11+)

Install Cool-Seq-Tool
---------------------

Install Cool-Seq-Tool from `PyPI <https://pypi.org/project/cool-seq-tool/>`_:

    pip install cool-seq-tool

.. _dependency-groups:

.. note::

   Cool-Seq-Tool provides extra dependency groups for development and testing purposes. Most users won't need to install them.

   * ``dev`` includes packages for linting and performing other kinds of quality checks
   * ``tests`` includes packages for running tests
   * ``docs`` includes packages for writing and building documentation

Set up UTA
----------

Cool-Seq-Tool requires an available instance of the Universal Transcript Archive (UTA) database. Complete installation instructions (via Docker or a local server) are available at the `UTA GitHub repository <https://github.com/biocommons/uta>`_. For local usage, we recommend the following:

.. long-term, it would be best to move this over to the UTA repo to avoid duplication

.. code-block::

   createuser -U postgres uta_admin
   createuser -U postgres anonymous
   createdb -U postgres -O uta_admin uta

   export UTA_VERSION=uta_20210129b.pgd.gz  # most recent as of 2023/12/05
   curl -O https://dl.biocommons.org/uta/$UTA_VERSION
   gzip -cdq ${UTA_VERSION} | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5433

By default, Cool-Seq-Tool expects to connect to the UTA database via a PostgreSQL connection served local on port 5433, under the PostgreSQL username ``uta_admin`` and the schema ``uta_20210129b``.

Set up SeqRepo
--------------

Cool-Seq-Tool requires access to `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ data. In general, we recommend the following for local setup:

.. long-term, it would be best to move this over to seqrepo to avoid duplication

.. code-block::

   pip install seqrepo
   export SEQREPO_VERSION=2021-01-29  # or newer if available -- check `seqrepo list-remote-instances`
   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i $SEQREPO_VERSION

If you encounter a permission error similar to the one below:

.. code-block::

   PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2021-01-29._fkuefgd' -> '/usr/local/share/seqrepo/2021-01-29'

Try moving data manually with ``sudo``:

.. code-block::

   sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION

See `mirroring documentation <https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/mirror.rst>`_ on the SeqRepo GitHub repo for instructions and additional troubleshooting.
