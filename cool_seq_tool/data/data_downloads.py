"""Module for downloading data files."""
from ftplib import FTP
from os import remove
import gzip
import shutil
import datetime

from dateutil import parser
from apybiomart import query

from cool_seq_tool import APP_ROOT


class DataDownload:
    """Class for downloading data files."""

    def __init__(self) -> None:
        """Initialize DataDownload."""
        self._make_data_dir()

    def _make_data_dir(self) -> None:
        """Make data directory"""
        self._data_dir = APP_ROOT / "data"
        self._data_dir.mkdir(exist_ok=True, parents=True)

    def download_transcript_mappings_data(self) -> None:
        """Acquire transcript mapping data from Ensembl BioMart."""
        transcript_file_path = self._data_dir / "transcript_mapping.tsv"
        result = query(
            dataset="hsapiens_gene_ensembl",
            filters={},
            attributes=[
                "ensembl_gene_id",
                "ensembl_gene_id_version",
                "ensembl_transcript_id",
                "ensembl_transcript_id_version",
                "ensembl_peptide_id",
                "ensembl_peptide_id_version",
                "transcript_mane_select",
                "external_gene_name"
            ]
        )
        result.to_csv(transcript_file_path, sep="\t")

    def _download_mane_summary(self) -> None:
        """Download latest MANE summary data file."""
        with FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            ftp.login()
            ftp.cwd("/refseq/MANE/MANE_human/current")
            files = ftp.nlst()
            mane_summary_file = \
                [f for f in files if f.endswith(".summary.txt.gz")]
            if not mane_summary_file:
                raise Exception("Unable to download MANE summary data")
            mane_summary_file = mane_summary_file[0]
            self._mane_summary_path = \
                self._data_dir / mane_summary_file[:-3]
            mane_data_path = self._data_dir / mane_summary_file
            if not self._mane_summary_path.exists():
                with open(mane_data_path, "wb") as fp:
                    ftp.retrbinary(f"RETR {mane_summary_file}", fp.write)
                with gzip.open(mane_data_path, "rb") as f_in:
                    with open(self._mane_summary_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(mane_data_path)

    def _download_lrg_refseq_gene_data(self) -> None:
        """Download latest LRG_RefSeqGene file."""
        with FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            ftp.login()
            lrg_refseqgene_file = "LRG_RefSeqGene"
            ftp_dir_path = "/refseq/H_sapiens/RefSeqGene/"
            ftp_file_path = f"{ftp_dir_path}{lrg_refseqgene_file}"
            timestamp = ftp.voidcmd(f"MDTM {ftp_file_path}")[4:].strip()
            date = str(parser.parse(timestamp)).split()[0]
            version = datetime.datetime.strptime(date,
                                                 "%Y-%m-%d").strftime("%Y%m%d")

            fn_versioned = f"{lrg_refseqgene_file}_{version}"
            lrg_refseqgene_path = self._data_dir / lrg_refseqgene_file
            self._lrg_refseqgene_path = self._data_dir / fn_versioned
            if not self._lrg_refseqgene_path.exists():
                ftp.cwd(ftp_dir_path)
                with open(lrg_refseqgene_path, "wb") as fp:
                    ftp.retrbinary(f"RETR {lrg_refseqgene_file}", fp.write)
                with open(lrg_refseqgene_path, "rb") as f_in:
                    with open(self._lrg_refseqgene_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(lrg_refseqgene_path)
