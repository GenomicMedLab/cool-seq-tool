.. _transcript-policy:

Transcript selection policy
===========================

Often, user transcript queries match to multiple known transcripts. This section describes our decisionmaking policy for selecting a representative transcript using sequence attributes and MANE annotations from the :py:mod:`cool_seq_tool.data_sources.mane_transcript` module.

For more information on the MANE project, see the `RefSeq documentation <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_.

Transcript compatibility
------------------------

A compatible transcript must pass the following validation checks:

.. TODO more descriptive for each

* Validating that the position exists on an accession
* Validating against reference sequences
* Validating exon structure

Representative transcript priority
----------------------------------

.. TODO reword

We evaluate all compatible transcripts against each of the below criteria, and select as representative the transcript which meets the earliest criterion.

1. Transcript is annotated as a MANE Select transcript
2. Transcript is annotated as a MANE Plus Clinical transcript
3. Longest compatible remaining

   a. If there is a tie, choose the first-published transcript (lowest-numbered accession for RefSeq/Ensembl) among those transcripts meeting this criterion

.. note::

   .. TODO clarify?

   We want the most recent version of a transcript associated with an assembly.
