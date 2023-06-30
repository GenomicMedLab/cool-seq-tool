cool-seq-tool
=============

``cool-seq-tool`` (**C**\ ommon **O**\ perations **o**\ n **L**\ ots of **Seq**\ uences **Tool**) provides:

 - Transcript alignment data from the `Universal Transcript Archive <https://github.com/biocommons/uta>`_
 - Fast access to sequence data using `SeqRepo <https://github.com/biocommons/biocommons.seqrepo>`_
 - Liftover between assemblies (GRCh38 <--> GRCh37) from `PyLiftover <https://github.com/konstantint/pyliftover>`_
 - Lifting over to a preferred `MANE <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_ compatible transcript.

.. code-block:: pycon

   >>> import cool_seq_tool
   >>> print("example")
   "example"

``cool-seq-tool`` is a library created to support efforts to `normalize variation descriptions <https://github.com/cancervariants/variation-normalization/>`_ and `model gene fusions <https://github.com/cancervariants/fusion-curation/>`_. It is developed primarily by the `Wagner Lab <https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab>`_. Full source code is available on `GitHub <https://github.com/GenomicMedLab/cool-seq-tool>`_.

.. TODO link to the live fusion instance instead?

.. toctree::
   :hidden:
   :maxdepth: 2

    Installation<install>
    Usage<usage>
    Transcript Policy<transcript_policy>
    API Reference<api/index>
    Contributing<contributing>
    License<license>
