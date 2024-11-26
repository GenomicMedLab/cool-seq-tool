Cool-Seq-Tool |version|
=======================

.. image:: https://img.shields.io/pypi/v/cool-seq-tool.svg
   :alt: PyPI version
   :target: https://pypi.python.org/pypi/cool-seq-tool

.. image:: https://img.shields.io/pypi/l/cool-seq-tool.svg
   :alt: License
   :target: https://github.com/genomicmedlab/cool-seq-tool/blob/main/LICENSE

.. image:: https://img.shields.io/pypi/pyversions/cool-seq-tool?color=gr
   :alt: PyPI - supported Python versions

.. image:: https://github.com/genomicmedlab/cool-seq-tool/actions/workflows/checks.yaml/badge.svg
   :alt: tests status
   :target: https://github.com/genomicmedlab/cool-seq-tool/actions/workflows/checks.yaml

The Common Operations On Lots-Of Sequences Tool, **Cool-Seq-Tool**, provides:

* A Pythonic API on top of sequence data of interest to tertiary analysis tools, including mappings between gene names and transcripts, `MANE transcript <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ descriptions, and transcript alignment data from the `Universal Transcript Archive <https://github.com/biocommons/uta>`_
* Augmented access to the `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_ database, including multiple additional methods and tools
* Mapping tools, including a transcript selection algorithm for selecting a representative transcript defined :ref:`here <transcript_selection_policy>`, that combine the above to support translation between various references sequences and annotation layers, and transcripts

See the :ref:`Installation <installation>` and :ref:`Usage <usage>` pages for information on getting started. Individual classes and methods are documented within the :ref:`API reference <api_reference>`.

Cool-Seq-Tool was created to support the `Knowledgebase Integration Project <https://cancervariants.org/projects/integration/>`_ of the `Variant Interpretation for Cancer Consortium (VICC) <https://cancervariants.org/>`_. It is developed primarily by the `Wagner Lab <https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab>`_. Full source code is available on `GitHub <https://github.com/genomicmedlab/cool-seq-tool>`_.

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
