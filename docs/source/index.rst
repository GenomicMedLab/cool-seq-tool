Cool-Seq-Tool |version|
=======================

The **CoolSeqTool** provides:

* A Pythonic API on top of sequence data of interest to tertiary analysis tools, including mappings between gene names and transcripts, `MANE transcript <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ descriptions, and the `Universal Transcript Archive <https://github.com/biocommons/uta>`_
* Augmented access to the `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ database, including multiple additional methods and tools
* Mapping tools that combine the above to support translation between various references sequences and annotation layers, and to MANE-designated transcripts

See the :ref:`Installation <installation>` and :ref:`Usage <usage>` pages for information on getting started. Individual classes and methods are documented within the :ref:`API reference <api_reference>`.

CoolSeqTool was created to support the `Knowledgebase Integration Project <https://cancervariants.org/projects/integration/>`_ of the `Variant Interpretation for Cancer Consortium (VICC) <https://cancervariants.org/>`_. It is developed primarily by the `Wagner Lab <https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab>`_. Full source code is available on `GitHub <https://github.com/genomicmedlab/cool-seq-tool>`_.

.. toctree::
   :hidden:
   :maxdepth: 2

    Installation<install>
    Usage<usage>
    Transcript Selection<transcript_selection>
    API Reference<reference/index>
    Contributing<contributing>
    Changelog<changelog>
    License<license>
