.. _usage:

Usage
=====

Cool-Seq-Tool provides easy access to, and useful operations on, a selection of important genomic resources. Modules are divided into three groups:

* :ref:`Data sources <sources_modules_api_index>`, for basic acquisition and setup for a data source via Python

* :ref:`Data handlers <handlers_modules_api_index>`, for additional operations on top of existing sources

* :ref:`Data mappers <mappers_modules_api_index>`, for functions that incorporate multiple sources/handlers to produce output

The core :py:class:`CoolSeqTool <cool_seq_tool.app.CoolSeqTool>` class encapsulates all of their functions and can be used for easy initialization and access:

.. code-block:: pycon

   >>> from cool_seq_tool import CoolSeqTool
   >>> cst = CoolSeqTool()
   >>> cst.seqrepo_access.translate_alias("NM_002529.3")[0][-1]
   'ga4gh:SQ.RSkww1aYmsMiWbNdNnOTnVDAM3ZWp1uA'
   >>> cst.transcript_mappings.ensembl_protein_for_gene_symbol["BRAF"][0]
   'ENSP00000419060'
   >>> await cst.uta_db.get_ac_from_gene("BRAF")
   ['NC_000007.14', 'NC_000007.13']

Descriptions and examples of functions can be found in the :ref:`API Reference <api_reference>` section.

.. _async_note:

.. note::

   Many component classes in Cool-Seq-Tool, including :py:class:`UtaDatabase <cool_seq_tool.sources.uta_database.UtaDatabase>`, :py:class:`ExonGenomicCoordsMapper <cool_seq_tool.mappers.exon_genomic_coords.ExonGenomicCoordsMapper>`, and :py:class:`ManeTranscript <cool_seq_tool.mappers.mane_transcript>`, define public methods as ``async``. This means that, when used inside another function, they must be called with ``await``:

   .. code-block:: python

       from cool_seq_tool import CoolSeqTool

       async def do_thing():
           mane_mapper = CoolSeqTool().mane_transcript
           result = mane_mapper.g_to_grch38("NC_000001.11", 100, 200)
           print(type(result))
           # <class 'coroutine'>
           awaited_result = await result
           print(awaited_result)
           # {'ac': 'NC_000001.11', 'pos': (100, 200)}

   In a REPL, ``asyncio.run()`` can be used to call coroutines outside of functions. Many of our docstring examples will use this pattern.

   .. code-block:: pycon

      >>> import asyncio
      >>> from cool_seq_tool import cool_seq_tool
      >>> mane_mapper = CoolSeqTool().mane_transcript
      >>> result = asyncio.run(mane_mapper.g_to_grch38("NC_000001.11", 100, 200))
      >>> print(result)
      {'ac': 'NC_000001.11', 'pos': (100, 200)}

   See the `asyncio module documentation <https://docs.python.org/3/library/asyncio.html>`_ for more information.

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
     - Path to LRG_RefSeqGene file. Used in :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>` to provide mappings between gene symbols and RefSeq/Ensembl transcript accessions. If not defined, uses `wags-tails <https://wags-tails.readthedocs.io/stable/usage.html#configuration>`_ to fetch the latest version, downloading it from the NCBI server if necessary.
   * - ``TRANSCRIPT_MAPPINGS_PATH``
     - Path to transcript mapping file generated from `Ensembl BioMart <http://www.ensembl.org/biomart/martview>`_. Used in :py:class:`TranscriptMappings <cool_seq_tool.sources.transcript_mappings.TranscriptMappings>`. If not defined, uses a copy of the file that is bundled within the Cool-Seq-Tool installation. See the :ref:`contributor instructions <build_transcript_mappings_tsv>` for information on manually rebuilding it.
   * - ``MANE_SUMMARY_PATH``
     - Path to MANE Summary file. Used in :py:class:`ManeTranscriptMappings <cool_seq_tool.sources.mane_transcript_mappings.ManeTranscriptMappings>` to provide MANE transcript annotations. If not defined, uses `wags-tails <https://wags-tails.readthedocs.io/stable/usage.html#configuration>`_ to fetch the latest version, downloading it from the NCBI server if necessary.
   * - ``SEQREPO_ROOT_DIR``
     - Path to SeqRepo directory (i.e. contains ``aliases.sqlite3`` database file, and ``sequences`` directory). Used by :py:class:`SeqRepoAccess <cool_seq_tool.handlers.seqrepo_access.SeqRepoAccess>`. If not defined, defaults to ``/usr/local/share/seqrepo/latest``.
   * - ``UTA_DB_URL``
     - A `libpq connection string <https://www.postgresql.org/docs/current/libpq-connect.html#LIBPQ-CONNSTRING>`_, i.e. of the form ``postgresql://<user>:<password>@<host>:<port>/<database>/<schema>``, used by the :py:class:`UtaDatabase <cool_seq_tool.sources.uta_database.UtaDatabase>` class. By default, it is set to ``postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220``.
   * - ``LIFTOVER_CHAIN_37_TO_38``
     - A path to a `chainfile <https://genome.ucsc.edu/goldenPath/help/chain.html>`_ for lifting from GRCh37 to GRCh38. Used by the :py:class:`LiftOver <cool_seq_tool.mappers.liftover.LiftOver>` class as input to `agct <https://pypi.org/project/agct/>`_. If not provided, agct will fetch it automatically from UCSC.
   * - ``LIFTOVER_CHAIN_38_TO_37``
     - A path to a `chainfile <https://genome.ucsc.edu/goldenPath/help/chain.html>`_ for lifting from GRCh38 to GRCh37. Used by the :py:class:`LiftOver <cool_seq_tool.mappers.liftover.LiftOver>` class as input to `agct <https://pypi.org/project/agct/>`_. If not provided, agct will fetch it automatically from UCSC.
