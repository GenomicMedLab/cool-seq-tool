Transcript Selection Policy
===========================

Often, user transcript queries match to multiple known transcripts. This section describes our decisionmaking policy for selecting a representative transcript using sequence attributes and MANE annotations from the :py:mod:`cool_seq_tool.data_sources.mane_transcript` module.

For more information on the MANE project, see the `RefSeq documentation <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_.

Representative transcript priority
----------------------------------

We evaluate all compatible transcripts against each of the below criteria, and select as representative the transcript which meets the earliest criterion.

1. Transcript is annotated as a MANE Select transcript
2. Transcript is annotated as a MANE Plus Clinical transcript
3. Longest compatible remaining

   a. If there is a tie, choose the first-published transcript (lowest-numbered accession for RefSeq/Ensembl) among those transcripts meeting this criterion

.. note::

   .. TODO explain more?

   We want the most recent version of a transcript associated with an assembly.

Compatible Transcripts
----------------------

Compatible transcripts are those that pass validation checks. The checks that we make are:

* Validating that the position exists on an accession
* Validating against reference sequences
* Validating exon structure
