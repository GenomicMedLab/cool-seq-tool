.. _installation:

Installation
============

Prerequisites
-------------

* Python 3.10+
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

   export UTA_VERSION=uta_20241220.pgd.gz  # most recent as of 2025/03/10
   curl -O https://dl.biocommons.org/uta/$UTA_VERSION
   gzip -cdq ${UTA_VERSION} | psql -h localhost -U uta_admin --echo-errors --single-transaction -v ON_ERROR_STOP=1 -d uta -p 5432

By default, Cool-Seq-Tool expects to connect to the UTA database via a PostgreSQL connection served local on port 5432, under the PostgreSQL username ``uta_admin`` and the schema ``uta_20241220``.

Set up SeqRepo
--------------

Cool-Seq-Tool requires access to `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ data. In general, we recommend the following for local setup:

.. long-term, it would be best to move this over to seqrepo to avoid duplication

.. code-block::

   pip install seqrepo
   export SEQREPO_VERSION=2024-12-20
   sudo mkdir /usr/local/share/seqrepo
   sudo chown $USER /usr/local/share/seqrepo
   seqrepo pull -i $SEQREPO_VERSION

.. note::

   Our lab typically uses the latest SeqRepo release, which is ``2024-12-20`` as of this commit. To check for the presence of newer snapshots, use the ``seqrepo list-remote-instances`` CLI command.

While this should no longer occur with the latest SeqRepo release, some users in the past have reported experiencing the following error:

.. code-block::

   PermissionError: [Error 13] Permission denied: '/usr/local/share/seqrepo/2024-12-20._fkuefgd' -> '/usr/local/share/seqrepo/2024-12-20'

Try moving data manually with ``sudo``:

.. code-block::

   sudo mv /usr/local/share/seqrepo/$SEQREPO_VERSION.* /usr/local/share/seqrepo/$SEQREPO_VERSION

See `mirroring documentation <https://github.com/biocommons/biocommons.seqrepo/blob/main/docs/mirror.rst>`_ on the SeqRepo GitHub repo for instructions and additional troubleshooting.

Check data availability
-----------------------

The :py:meth:`check_status <cool_seq_tool.resources.status.check_status>` method is provided to check that data dependencies like SeqRepo and UTA are available and configured correctly:

.. code-block:: pycon

   >>> from cool_seq_tool.resources.status import check_status
   >>> import asyncio
   >>> asyncio.run(check_status())
   ResourceStatus(uta=True, seqrepo=True, transcript_mappings=True, mane_summary=True, lrg_refseqgene=True, liftover=True)
