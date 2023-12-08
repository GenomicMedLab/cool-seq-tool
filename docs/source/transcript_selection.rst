Transcript Selection
====================

One of the core uses of Cool-Seq-Tool is to acquire and use consensus, representative transcripts in performing genomic analysis. Here, we describe the selection processes, programmed in the :py:class:`MANETranscript <cool_seq_tool.mappers.mane_transcript.MANETranscript>` class, for choosing the best available transcripts that are compatible with requested data.

We rely heavily on transcripts annotated under the `Matched Annotation from NCBI and EMBL-EBI (MANE)` Transcripts project. For more information on the MANE project, see the `NCBI MANE page <https://www.ncbi.nlm.nih.gov/refseq/MANE/>`_.

Transcript compatibility
------------------------

The following validation checks are performed to determine compatibility of a transcript position:

* The position exists on an accession
* "Validating reference sequences"  # TODO what does this mean?
* "Validating exon structure"  # TODO what does this mean?

A transcript that fails to pass any of these checks is discarded as incompatible.

Representative transcript priority
----------------------------------

All compatible transcripts are evaluated and ordered against the below criteria. The candidate transcript which meets the earliest criterion is chosen as representative.

#. Transcript is annotated as a `MANE Select` transcript
#. Transcript is annotated as a `MANE Plus Clinical` transcript
#. Transcript is the longest-compatible remaining transcript
#. Transcript is the first-published (lowest-numbered RefSeq/Ensembl accession) remaining transcript

.. note::

   We prefer the most recent version of a transcript associated with an assembly.
