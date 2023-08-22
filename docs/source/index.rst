cool-seq-tool |version|
=======================

``cool-seq-tool`` (**C**\ ommon **O**\ perations **o**\ n **L**\ ots of **Seq**\ uences **Tool**) provides utilities for sequence conversion and retrieval. Pooling data from sources like `the Universal Transcript Archive <#>`_, `SeqRepo <#>`_, `the Gene Normalizer <#>`_, and `others <#>`_, it enables consistent and unambiguous conversions between genomic coordinates and high-quality transcript annotations from the `MANE project <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_.

.. code-block:: pycon

   >>> from cool_seq_tool import CoolSeqTool
   >>> cst = CoolSeqTool()
   >>> result = await cst.transcript_to_genomic_coordinates(transcript="NM_002529.3", exon_start=1)
   >>> (result.genomic_data.gene, result.genomic_data.chr, result.genomic_data.start, result.genomic_data.strand)
   ('NTRK1', 'NC_000001.11', 156861146, 1)

``cool-seq-tool`` is a library created to support efforts to `normalize variation descriptions <https://github.com/cancervariants/variation-normalization/>`_ and `model gene fusions <https://cancervariants.org/projects/fusions>`_ under the mantle of the `Variation Interpretation for Cancer Consortium (VICC) <https://cancervariants.org>`_. It is developed primarily by the `Wagner Lab <https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab>`_. Full source code is available on `GitHub <https://github.com/GenomicMedLab/cool-seq-tool>`_.

.. toctree::
   :hidden:
   :maxdepth: 2

    Installation<install>
    Usage<usage>
    REST service<rest_service>
    Transcript policy<transcript_policy>
    Data sources<data_sources/index>
    API reference<api/index>
    Contributing<contributing>
    License<license>
