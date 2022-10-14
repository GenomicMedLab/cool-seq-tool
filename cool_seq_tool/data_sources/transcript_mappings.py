"""The module for Transcript Mappings."""
import csv
from typing import Dict, List, Optional

from cool_seq_tool import TRANSCRIPT_MAPPINGS_PATH, LRG_REFSEQGENE_PATH


class TranscriptMappings:
    """The transcript mappings class."""

    def __init__(self, transcript_file_path: str = TRANSCRIPT_MAPPINGS_PATH,
                 lrg_refseqgene_path: str = LRG_REFSEQGENE_PATH) -> None:
        """Initialize the transcript mappings class.

        :param str transcript_file_path: Path to transcript mappings file
        :param str lrg_refseqgene_path: Path to LRG RefSeqGene file
        """
        # ENSP <-> Gene Symbol
        self.ensembl_protein_version_for_gene_symbol: Dict[str, List[str]] = {}
        self.ensembl_protein_version_to_gene_symbol: Dict[str, str] = {}
        self.ensembl_protein_for_gene_symbol: Dict[str, List[str]] = {}
        self.ensembl_protein_to_gene_symbol: Dict[str, str] = {}

        # Gene Symbol <-> ENST
        self.ensembl_transcript_version_for_gene_symbol: \
            Dict[str, List[str]] = {}
        self.ensembl_transcript_version_to_gene_symbol: Dict[str, str] = {}
        self.ensembl_transcript_for_gene_symbol: Dict[str, List[str]] = {}
        self.ensembl_transcript_to_gene_symbol: Dict[str, str] = {}

        # NP_ <-> Gene Symbol
        self.refseq_protein_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_protein_to_gene_symbol: Dict[str, str] = {}

        # NM_ <-> Gene Symbol
        self.refseq_rna_version_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_rna_version_to_gene_symbol: Dict[str, str] = {}
        self.refseq_rna_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_rna_to_gene_symbol: Dict[str, str] = {}

        # LRG <-> Gene Symbol
        self.refseq_lrg_for_gene_symbol: Dict[str, List[str]] = {}
        self.refseq_lrg_to_gene_symbol: Dict[str, str] = {}

        # NP -> NM
        self.np_to_nm: Dict[str, str] = {}

        # ENSP -> ENST
        self.ensp_to_enst: Dict[str, str] = {}

        self._load_transcript_mappings_data(transcript_file_path)
        self._load_refseq_gene_symbol_data(lrg_refseqgene_path)

    def _load_transcript_mappings_data(self,
                                       transcript_file_path: str) -> None:
        """Load transcript mappings file to dictionaries.

        :param str transcript_file_path: Path to transcript mappings file
        """
        with open(transcript_file_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row["Gene name"]
                if gene:
                    versioned_protein_transcript = \
                        row["Protein stable ID version"]
                    if versioned_protein_transcript:
                        self.ensembl_protein_version_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(versioned_protein_transcript)
                        self.ensembl_protein_version_to_gene_symbol[
                            versioned_protein_transcript] = gene
                    protein_transcript = row["Protein stable ID"]
                    if protein_transcript:
                        self.ensembl_protein_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(protein_transcript)
                        self.ensembl_protein_to_gene_symbol[
                            protein_transcript] = gene
                    versioned_transcript = \
                        row["Transcript stable ID version"]
                    if versioned_transcript:
                        self.ensembl_transcript_version_for_gene_symbol \
                            .setdefault(gene, []) \
                            .append(versioned_transcript)
                        self.ensembl_transcript_version_to_gene_symbol[
                            versioned_transcript] = gene
                    transcript = row["Transcript stable ID"]
                    if transcript:
                        self.ensembl_transcript_for_gene_symbol\
                            .setdefault(gene, []) \
                            .append(transcript)
                        self.ensembl_transcript_to_gene_symbol[
                            transcript] = gene
                    if versioned_transcript and versioned_protein_transcript:
                        self.ensp_to_enst[versioned_protein_transcript] = \
                            versioned_transcript

    def _load_refseq_gene_symbol_data(self, lrg_refseqgene_path: str) -> None:
        """Load data from RefSeq Gene Symbol file to dictionaries.

        :param str lrg_refseqgene_path: Path to LRG RefSeqGene file
        """
        with open(lrg_refseqgene_path) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row["Symbol"]
                if gene:
                    refseq_transcript = row["Protein"]
                    if refseq_transcript:
                        self.refseq_protein_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(refseq_transcript)
                        self.refseq_protein_to_gene_symbol[
                            refseq_transcript] = gene
                    rna_transcript = row["RNA"]
                    if rna_transcript:
                        self.refseq_rna_version_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(rna_transcript)
                        self.refseq_rna_version_to_gene_symbol[
                            rna_transcript] = gene
                        if "." in rna_transcript:
                            rna_t = rna_transcript.split(".")[0]
                            self.refseq_rna_for_gene_symbol.\
                                setdefault(gene, []).\
                                append(rna_t)
                            self.refseq_rna_to_gene_symbol[
                                rna_t] = gene
                    if refseq_transcript and rna_transcript:
                        self.np_to_nm[refseq_transcript] = rna_transcript
                    lrg = row["LRG"]
                    if lrg:
                        self.refseq_lrg_for_gene_symbol.\
                            setdefault(gene, []).\
                            append(lrg)
                        self.refseq_lrg_to_gene_symbol[lrg] = gene
                        t = row["t"]
                        p = row["p"]
                        if t:
                            t_lrg = lrg + t
                            self.refseq_lrg_for_gene_symbol. \
                                setdefault(gene, []). \
                                append(t_lrg)
                            self.refseq_lrg_to_gene_symbol[t_lrg] = gene
                        if p:
                            p_lrg = lrg + p
                            self.refseq_lrg_for_gene_symbol. \
                                setdefault(gene, []). \
                                append(p_lrg)
                            self.refseq_lrg_to_gene_symbol[p_lrg] = gene

    def protein_transcripts(self, identifier: str) -> List[str]:
        """Return a list of protein transcripts for a gene symbol.

        :param str identifier: Gene identifier to get protein transcripts for
        :return: Protein transcripts for a gene symbol
        """
        protein_transcripts = list()
        protein_transcripts += \
            self.ensembl_protein_version_for_gene_symbol.get(
                identifier, "")
        protein_transcripts += \
            self.ensembl_protein_for_gene_symbol.get(identifier, "")
        protein_transcripts += \
            self.refseq_protein_for_gene_symbol.get(identifier, "")
        return list(set(protein_transcripts))

    def coding_dna_transcripts(self, identifier: str) -> List[str]:
        """Return transcripts from a coding dna refseq for a gene symbol.

        :param str identifier: Gene identifier to find transcripts for
        :return: cDNA transcripts for a gene symbol
        """
        genomic_transcripts = list()
        genomic_transcripts += \
            self.ensembl_transcript_version_for_gene_symbol.get(identifier,
                                                                "")
        genomic_transcripts += \
            self.refseq_rna_version_for_gene_symbol.get(identifier, "")
        genomic_transcripts += \
            self.refseq_rna_version_for_gene_symbol.get(identifier, "")
        return list(set(genomic_transcripts))

    def get_gene_symbol_from_ensembl_protein(self, q: str) -> Optional[str]:
        """Return the gene symbol for a Ensembl Protein.

        :param str q: ensembl protein accession
        :return: Gene symbol
        """
        gene_symbol = self.ensembl_protein_version_to_gene_symbol.get(q)
        if not gene_symbol:
            if "." in q:
                q = q.split(".")[0]
                gene_symbol = self.ensembl_protein_to_gene_symbol.get(q)
        return gene_symbol

    def get_gene_symbol_from_refeq_protein(self, q: str) -> Optional[str]:
        """Return the gene symbol for a Refseq Protein.

        :param str q: RefSeq protein accession
        :return: Gene symbol
        """
        return self.refseq_protein_to_gene_symbol.get(q)

    def get_gene_symbol_from_refseq_rna(self, q: str) -> Optional[str]:
        """Return gene symbol for a Refseq RNA Transcript.

        :param str q: RefSeq RNA transcript accession
        :return: Gene symbol
        """
        gene_symbol = self.refseq_rna_version_to_gene_symbol.get(q)
        if not gene_symbol:
            if "." in q:
                q = q.split(".")[0]
                gene_symbol = self.refseq_rna_to_gene_symbol.get(q)
        return gene_symbol

    def get_gene_symbol_from_ensembl_transcript(self, q: str) -> Optional[str]:
        """Return gene symbol for an Ensembl Transcript.

        :param str q: Ensembl transcript accession
        :return: Gene symbol
        """
        gene_symbol = self.ensembl_transcript_version_to_gene_symbol.get(q)
        if not gene_symbol:
            if "." in q:
                q = q.split(".")[0]
                gene_symbol = self.ensembl_transcript_to_gene_symbol.get(q)
        return gene_symbol

    def get_gene_symbol_from_lrg(self, q: str) -> Optional[str]:
        """Return gene symbol for LRG.

        :param str q: LRG accession
        :return: Gene symbol
        """
        return self.refseq_lrg_to_gene_symbol.get(q)
