"""Provide paths to shared files, and trigger data acquisition if unavailable."""
from os import environ
from pathlib import Path

from cool_seq_tool.data import DataDownload


d = DataDownload()

provided_transcript_mappings_path = environ.get("TRANSCRIPT_MAPPINGS_PATH")
if provided_transcript_mappings_path:
    TRANSCRIPT_MAPPINGS_PATH = Path(provided_transcript_mappings_path)
else:
    TRANSCRIPT_MAPPINGS_PATH = d._data_dir / "transcript_mapping.tsv"

provided_mane_summary_path = environ.get("MANE_SUMMARY_PATH", "")
if provided_mane_summary_path:
    MANE_SUMMARY_PATH = Path(provided_mane_summary_path)
else:
    MANE_SUMMARY_PATH = d.get_mane_summary()

provided_lrg_refseq_path = environ.get("LRG_REFSEQGENE_PATH", "")
if provided_lrg_refseq_path:
    LRG_REFSEQGENE_PATH = Path(provided_lrg_refseq_path)
else:
    LRG_REFSEQGENE_PATH = d.get_lrg_refseq_gene_data()

SEQREPO_ROOT_DIR = environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")
