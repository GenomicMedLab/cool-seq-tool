"""Handle acquisition of external data."""
import datetime
import gzip
import logging
import shutil
from ftplib import FTP
from os import remove
from pathlib import Path

from dateutil import parser

from cool_seq_tool import APP_ROOT

logger = logging.getLogger("cool_seq_tool")


class DataDownload:
    """Manage downloadable data files. Responsible for checking if files are available
    under expected locations, and fetching them if not.
    """

    def __init__(self) -> None:
        """Initialize downloadable data locations."""
        self._data_dir = APP_ROOT / "data"

    def get_mane_summary(self) -> Path:
        """Identify latest MANE summary data. If unavailable locally, download from
        `NCBI FTP server <https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/>`_.

        :return: path to MANE summary file
        """
        with FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            ftp.login()
            ftp.cwd("/refseq/MANE/MANE_human/current")
            files = ftp.nlst()
            mane_summary_file = [f for f in files if f.endswith(".summary.txt.gz")]
            if not mane_summary_file:
                raise Exception("Unable to download MANE summary data")
            mane_summary_file = mane_summary_file[0]
            self._mane_summary_path = self._data_dir / mane_summary_file[:-3]
            mane_data_path = self._data_dir / mane_summary_file
            if not self._mane_summary_path.exists():
                logger.info("Downloading MANE summary file from NCBI.")
                with open(mane_data_path, "wb") as fp:
                    ftp.retrbinary(f"RETR {mane_summary_file}", fp.write)
                with gzip.open(mane_data_path, "rb") as f_in:
                    with open(self._mane_summary_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(mane_data_path)
                logger.info("MANE summary file download complete.")
        return self._mane_summary_path

    def get_lrg_refseq_gene_data(self) -> Path:
        """Identify latest LRG RefSeq Gene file. If unavailable locally, download from
        `NCBI FTP server <https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/RefSeqGene/>`_.

        :return: path to acquired LRG RefSeq Gene data file
        """
        with FTP("ftp.ncbi.nlm.nih.gov") as ftp:
            ftp.login()
            lrg_refseqgene_file = "LRG_RefSeqGene"
            ftp_dir_path = "/refseq/H_sapiens/RefSeqGene/"
            ftp_file_path = f"{ftp_dir_path}{lrg_refseqgene_file}"
            timestamp = ftp.voidcmd(f"MDTM {ftp_file_path}")[4:].strip()
            date = str(parser.parse(timestamp)).split()[0]
            version = datetime.datetime.strptime(date, "%Y-%m-%d").strftime("%Y%m%d")
            fn_versioned = f"{lrg_refseqgene_file}_{version}"
            lrg_refseqgene_path = self._data_dir / lrg_refseqgene_file
            self._lrg_refseqgene_path = self._data_dir / fn_versioned
            if not self._lrg_refseqgene_path.exists():
                logger.info("Downloading LRG RefSeq data from NCBI.")
                ftp.cwd(ftp_dir_path)
                with open(lrg_refseqgene_path, "wb") as fp:
                    ftp.retrbinary(f"RETR {lrg_refseqgene_file}", fp.write)
                with open(lrg_refseqgene_path, "rb") as f_in:
                    with open(self._lrg_refseqgene_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
                remove(lrg_refseqgene_path)
                logger.info("LRG RefSeq data download complete.")
        return self._lrg_refseqgene_path
