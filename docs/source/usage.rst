Usage
=====

Cool-Seq-Tool provides easy access to, and useful operations on, a selection of important genomic resources. Modules are divided into three groups:

* :ref:`Data sources <sources_modules_api_index>`, for basic acquisition and setup for a data source via Python

* :ref:`Data handlers <handlers_modules_api_index>`, for additional operations on top of existing sources

* :ref:`Data mappers <mappers_modules_api_index>`, for functions that incorporate multiple sources/handlers to produce output

A core :py:class:`CoolSeqTool <cool_seq_tool.app.CoolSeqTool>` class encapsulates all of their functions and can be used for easy initialization:

.. code-block:: pycon

   >>> from cool_seq_tool.app import CoolSeqTool
   >>> cst = CoolSeqTool()

Descriptions and examples of functions can be found in the :ref:`API Reference <api_reference>` section.

REST server
-----------

TODO Possibly staged for deletion?

.. _configuration:

Environment configuration
-------------------------

Individual classes will accept arguments upon initialization to set parameters regarding data sources. In general, these parameters are also configurable via environment variables, e.g. in a cloud deployment.

.. list-table::
   :widths: 25 100
   :header-rows: 1

   * - Variable
     - Description
   * - ``LRG_REFSEQGENE_PATH``
     - Path to LRG_RefSeqGene file. Used in :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>` to provide mappings between gene symbols and RefSeq/Ensembl transcript accessions. If not defined, defaults to the most recent version (formatted as ``data/LRG_RefSeqGene_YYYYMMDD``) within the Cool-Seq-Tool library directory.
   * - ``TRANSCRIPT_MAPPINGS_PATH``
     - Path to transcript mapping file generated from `Ensembl BioMart <http://www.ensembl.org/biomart/martview>`_. Used in :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>`. If not defined, uses a copy of the file that is bundled within the Cool-Seq-Tool installation. See the :ref:`contributor instructions <build_transcript_mappings_tsv>` for information on manually rebuilding it.
   * - ``MANE_SUMMARY_PATH``
     - Path to MANE Summary file. Used in :py:class:`MANETranscriptMappings <cool_seq_tool.sources.mane_transcript_mappings.MANETranscriptMappings>` to provide MANE transcript annotations. If not defined, defaults to the most recent version (formatted as ``data/MANE.GRCh38vX.X.summary.txt``) within the Cool-Seq-Tool library directory.
   * - ``SEQREPO_ROOT_DIR``
     - Path to SeqRepo directory (i.e. contains ``aliases.sqlite3`` database file, and ``sequences`` directory). Used by :py:class:`SeqRepoAccess <cool_seq_tool.handlers.seqrepo_access.SeqRepoAccess`. If not defined, defaults to ``/usr/local/share/seqrepo/latest``.
   * - ``UTA_DB_URL``
     - A `libpq connection string <https://www.postgresql.org/docs/current/libpq-connect.html#LIBPQ-CONNSTRING>`_, i.e. of the form ``postgresql://<user>:<password>@<host>:<port>/<database>/<schema>``, used by the :py:class:`cool_seq_tool.sources.uta_database.UTADatabase` class. By default, it is set to ``postgresql://uta_admin:uta@localhost:5433/uta/uta_20210129b``.
   * - ``LIFTOVER_CHAIN_37_TO_38``
     - A path to a `chainfile <https://genome.ucsc.edu/goldenPath/help/chain.html>`_ for lifting from GRCh37 to GRCh38. Used by :py:class:`cool_seq_tool.sources.uta_database.UTADatabase` as input to `pyliftover <https://pypi.org/project/pyliftover/>`_. If not provided, pyliftover will fetch it automatically from UCSC.
   * - ``LIFTOVER_CHAIN_38_TO_37``
     - A path to a `chainfile <https://genome.ucsc.edu/goldenPath/help/chain.html>`_ for lifting from GRCh38 to GRCh37. Used by :py:class:`cool_seq_tool.sources.uta_database.UTADatabase` as input to `pyliftover <https://pypi.org/project/pyliftover/>`_. If not provided, pyliftover will fetch it automatically from UCSC.

Schema support
--------------

Many genomic data objects produced by Cool-Seq-Tool are structured in conformance with the `Variation Representation Specification <https://vrs.ga4gh.org/en/stable/>`_, courtesy of the `VRS-Python <https://github.com/ga4gh/vrs-python>` library.
