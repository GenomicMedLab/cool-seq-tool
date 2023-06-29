Transcript Selection Policy
===========================

This document contains information on the Transcript Selection Policy. We use this policy for selecting a representative transcript using sequence attributes and MANE annotations from the `mane_transcript` module.

More information on MANE can be found [here](https://www.ncbi.nlm.nih.gov/refseq/MANE/).

Representative transcript priority
----------------------------------

We evaluate all compatible transcripts against each of the below criteria, and select the transcript which meets the earliest criterion as representative.

1. Transcript is annotated as a MANE Select transcript
2. Transcript is annotated as a MANE Plus Clinical transcript
3. Longest Compatible Remaining
   a. If there is a tie, choose the first-published transcript (lowest-numbered accession for RefSeq/Ensembl) among those transcripts meeting this criterion

.. note::

   We want the most recent version of a transcript associated with an **assembly**.

Compatible Transcripts
----------------------

Compatible transcripts are those that pass validation checks. The checks that we make are:

   - Validating the position exists on an accession
   - Validating reference sequences
   - Validating exon structure

