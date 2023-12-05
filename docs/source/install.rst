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

Acquire genome and transcript data
----------------------------------

Cool-Seq-Tool also makes use of external data extracted from a handful of sources. This data is fetched automatically at runtime if it's not available, so no additional action is necessary by default -- however, if necessary, it can be configured and placed manually.

LRG_RefSeqGene
++++++++++++++

Cool-Seq-Tool fetches the latest version of ``LRG_RefSeqGene`` if the environment variable ``LRG_REFSEQGENE_PATH`` is not set. When ``LRG_REFSEQGENE_PATH`` is set, Cool-Seq-Tool will look at this path and expect the ``LRG_RefSeqGene`` file. This file is found can be found `here <https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene>`_.

MANE summary
++++++++++++

Cool-Seq-Tool fetches the latest version of ``MANE.GRCh38.*.summary.txt.gz`` if the environment variable ``MANE_SUMMARY_PATH`` is not set. When ``MANE_SUMMARY_PATH`` is set, Cool-Seq-Tool will look at this path and expect the MANE Summary Data file. This file is found can be found `here <https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/>`_.

Transcript mappings
+++++++++++++++++++

Cool-Seq-Tool is packaged with transcript mapping data acquired from `Ensembl BioMart <http://www.ensembl.org/biomart/martview>`_. If the environment variable ``TRANSCRIPT_MAPPINGS_PATH`` is not set, Cool-Seq-Tool will use the built-in file. When ``TRANSCRIPT_MAPPINGS_PATH`` is set, Cool-Seq-Tool will look at this path and expect to find the transcript mapping TSV file. To manually acquire this file, see the :ref:`contributor instructions <build_transcript_mappings_tsv>`.