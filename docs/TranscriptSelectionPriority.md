# Transcript Selection Priority

This document contains information on the Transcript Selection Priority. We use this priority when retrieving MANE data in the `mane_transcript` module.

More information on MANE can be found [here](https://www.ncbi.nlm.nih.gov/refseq/MANE/).

## Priority

1. MANE Select
2. MANE Plus Clinical
3. Longest Compatible Remaining\
   a. If there is a tie, choose the first-published transcript (lowest-numbered accession for RefSeq/Ensembl) among those transcripts meeting criterion\
   _Note: We want the most recent version of a transcript associated with an assembly_

## Compatible Transcripts

Compatible transcripts are those that pass validation checks. The checks that we make are:

   - Validating the position exists on an accession
   - Validating reference sequences
   - Validating exon structure
