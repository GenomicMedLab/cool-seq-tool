.. _transcript_selection_policy:

Transcript Selection
====================

One of the core uses of Cool-Seq-Tool is to acquire and use consensus-based, representative transcripts in performing genomic analysis. Here, we describe the selection processes, programmed in the :py:class:`ManeTranscript <cool_seq_tool.mappers.mane_transcript.ManeTranscript>` class, for choosing the best available transcripts that are compatible with requested data.

We rely heavily on transcripts annotated under the `Matched Annotation from NCBI and EMBL-EBI (MANE) Transcripts` project. For more information on the MANE project, see the `NCBI MANE page <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_.

.. _transcript_compatibility:

Transcript compatibility
------------------------

The following validation checks are performed to determine compatibility of a transcript position:

* The position exists on an accession
* The sequence matches the expected reference sequence for a given accession or position
* Exon numbering matches known exon structure

A transcript that fails to pass any of these checks is discarded as incompatible.

Representative transcript priority
----------------------------------

All compatible transcripts are evaluated and ordered against the below criteria. The candidate transcript which meets the earliest criterion is chosen as representative.

#. Transcript is annotated as a `MANE Select` transcript
#. Transcript is annotated as a `MANE Plus Clinical` transcript
#. Transcript is the longest-compatible remaining transcript

   #. If there is a tie, choose the first-published (lowest-numbered RefSeq/Ensembl accession) transcript

   .. note::

      We always prefer the most recent version of a transcript associated with an assembly.
