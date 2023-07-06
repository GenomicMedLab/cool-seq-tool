Using cool-seq-tool
===================

.. note::

   Sequence and transcript conversion methods are defined as ``async``, making use of Python's `coroutines library <https://docs.python.org/3/library/asyncio-task.html>`_. This means they must be ``await``\ ed. Users working in IPython or Jupyter should also see the `relevant documentation <https://ipython.readthedocs.io/en/stable/interactive/autoawait.html>`_ for more information.

Initialization
--------------

The ``CoolSeqTool`` class can be initialized with no arguments to use default locations for all data sources:

.. code-block:: python

   from cool_seq_tool import CoolSeqTool
   cst = CoolSeqTool()





Transcript to genomic coordinates
---------------------------------

reference to :py:meth:`ref <cool_seq_tool.app.CoolSeqTool.transcript_to_genomic_coordinates>`

.. code-block:: python

   from cool_seq_tool import CoolSeqTool
   cst = CoolSeqTool()
   result = await cst.transcript_to_genomic_coordinates(gene="BRAF", transcript="refseq:NM_004333.4",
                                                        exon_start=1, exon_end=5)
   print(result.genomic_data.chr)  # 'NC_000007.14'
   print(result.genomic_data.start)  # 140924764
   print(result.strand)  # -1

.. code-block::

   { 'genomic_data': { 'chr': 'NC_000007.14',
                       'end': 140807959,
                       'exon_end': 5,
                       'exon_end_offset': 0,
                       'exon_start': 1,
                       'exon_start_offset': 0,
                       'gene': 'BRAF',
                       'start': 140924764,
                       'strand': -1,
                       'transcript': 'NM_004333.4'},
     'service_meta': { 'name': 'cool_seq_tool',
                       'response_datetime': datetime.datetime(2023, 7, 6, 9, 29, 47, 359864),
                       'url': 'https://github.com/GenomicMedLab/cool-seq-tool',
                       'version': '0.1.14-dev0'},
     'warnings': []}


Genomic to transcript coordinates
---------------------------------

FASTA generation
----------------


.. _residue-mode:

Coordinate input
----------------

residue mode
