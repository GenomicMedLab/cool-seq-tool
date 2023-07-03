"""Provide paths to shared files, and trigger data acquisition if unavailable."""
from os import environ
from pathlib import Path


APP_ROOT = Path(__file__).resolve().parents[0]

TRANSCRIPT_MAPPINGS_PATH = Path(environ.get("TRANSCRIPT_MAPPINGS_PATH",
                                            f"{APP_ROOT}/data/transcript_mapping.tsv"))

MANE_SUMMARY_PATH = environ.get("MANE_SUMMARY_PATH", "")
LRG_REFSEQGENE_PATH = environ.get("LRG_REFSEQGENE_PATH", "")
if not all((MANE_SUMMARY_PATH, LRG_REFSEQGENE_PATH)):
    from cool_seq_tool.data import DataDownload  # noqa: E402, I202
    d = DataDownload()

    if not MANE_SUMMARY_PATH:
        MANE_SUMMARY_PATH = d._mane_summary_path

    if not LRG_REFSEQGENE_PATH:
        LRG_REFSEQGENE_PATH = d._lrg_refseqgene_path
MANE_SUMMARY_PATH = Path(MANE_SUMMARY_PATH)
LRG_REFSEQGENE_PATH = Path(LRG_REFSEQGENE_PATH)

SEQREPO_ROOT_DIR = environ.get("SEQREPO_ROOT_DIR", "/usr/local/share/seqrepo/latest")
