"""Provide mappings between gene symbols and RefSeq + Ensembl transcript accessions."""

import csv
from pathlib import Path

from cool_seq_tool.resources.data_files import DataFile, get_data_file


class TranscriptMappings:
    """Provide mappings between gene symbols and RefSeq + Ensembl transcript accessions.

    Uses ``LRG_RefSeqGene`` and ``transcript_mappings.tsv``, which will automatically
    be acquired if they aren't already available. See the
    :ref:`configuration <configuration>` section in the documentation for information
    about manual acquisition of data.

    In general, this class's methods expect to receive NCBI gene symbols, so users
    should be careful about the sourcing of their input in cases where terms are
    conflicted or ambiguous (which, to be fair, should be relatively rare).
    """

    def __init__(
        self,
        transcript_file_path: Path | None = None,
        lrg_refseqgene_path: Path | None = None,
        from_local: bool = False,
    ) -> None:
        """Initialize the transcript mappings class.

        :param transcript_file_path: Path to transcript mappings file
        :param lrg_refseqgene_path: Path to LRG RefSeqGene file
        :param from_local: if ``True``, don't check for or acquire latest version --
            just provide most recent locally available file, if possible, and raise
            error otherwise
        """
        # ENSP <-> Gene Symbol
        self.ensembl_protein_version_for_gene_symbol: dict[str, list[str]] = {}
        self.ensembl_protein_version_to_gene_symbol: dict[str, str] = {}
        self.ensembl_protein_for_gene_symbol: dict[str, list[str]] = {}
        self.ensembl_protein_to_gene_symbol: dict[str, str] = {}

        # Gene Symbol <-> ENST
        self.ensembl_transcript_version_for_gene_symbol: dict[str, list[str]] = {}
        self.ensembl_transcript_version_to_gene_symbol: dict[str, str] = {}
        self.ensembl_transcript_for_gene_symbol: dict[str, list[str]] = {}
        self.ensembl_transcript_to_gene_symbol: dict[str, str] = {}

        # NP_ <-> Gene Symbol
        self.refseq_protein_for_gene_symbol: dict[str, list[str]] = {}
        self.refseq_protein_to_gene_symbol: dict[str, str] = {}

        # NM_ <-> Gene Symbol
        self.refseq_rna_version_for_gene_symbol: dict[str, list[str]] = {}
        self.refseq_rna_version_to_gene_symbol: dict[str, str] = {}
        self.refseq_rna_for_gene_symbol: dict[str, list[str]] = {}
        self.refseq_rna_to_gene_symbol: dict[str, str] = {}

        # NP -> NM
        self.np_to_nm: dict[str, str] = {}

        # ENSP -> ENST
        self.ensp_to_enst: dict[str, str] = {}

        self._load_transcript_mappings_data(
            transcript_file_path
            or get_data_file(DataFile.TRANSCRIPT_MAPPINGS, from_local)
        )
        self._load_refseq_gene_symbol_data(
            lrg_refseqgene_path or get_data_file(DataFile.LRG_REFSEQGENE, from_local)
        )

    def _load_transcript_mappings_data(self, transcript_file_path: Path) -> None:
        """Load transcript mappings file to dictionaries.

        :param transcript_file_path: Path to transcript mappings file
        """
        with transcript_file_path.open() as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row["Gene name"]
                if gene:
                    versioned_protein_transcript = row["Protein stable ID version"]
                    if versioned_protein_transcript:
                        self.ensembl_protein_version_for_gene_symbol.setdefault(
                            gene, []
                        ).append(versioned_protein_transcript)
                        self.ensembl_protein_version_to_gene_symbol[
                            versioned_protein_transcript
                        ] = gene
                    protein_transcript = row["Protein stable ID"]
                    if protein_transcript:
                        self.ensembl_protein_for_gene_symbol.setdefault(
                            gene, []
                        ).append(protein_transcript)
                        self.ensembl_protein_to_gene_symbol[protein_transcript] = gene
                    versioned_transcript = row["Transcript stable ID version"]
                    if versioned_transcript:
                        self.ensembl_transcript_version_for_gene_symbol.setdefault(
                            gene, []
                        ).append(versioned_transcript)
                        self.ensembl_transcript_version_to_gene_symbol[
                            versioned_transcript
                        ] = gene
                    transcript = row["Transcript stable ID"]
                    if transcript:
                        self.ensembl_transcript_for_gene_symbol.setdefault(
                            gene, []
                        ).append(transcript)
                        self.ensembl_transcript_to_gene_symbol[transcript] = gene
                    if versioned_transcript and versioned_protein_transcript:
                        self.ensp_to_enst[versioned_protein_transcript] = (
                            versioned_transcript
                        )

    def _load_refseq_gene_symbol_data(self, lrg_refseqgene_path: Path) -> None:
        """Load data from RefSeq Gene Symbol file to dictionaries.

        :param Path lrg_refseqgene_path: Path to LRG RefSeqGene file
        """
        with lrg_refseqgene_path.open() as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                gene = row["Symbol"]
                if gene:
                    refseq_transcript = row["Protein"]
                    if refseq_transcript:
                        self.refseq_protein_for_gene_symbol.setdefault(gene, []).append(
                            refseq_transcript
                        )
                        self.refseq_protein_to_gene_symbol[refseq_transcript] = gene
                    rna_transcript = row["RNA"]
                    if rna_transcript:
                        self.refseq_rna_version_for_gene_symbol.setdefault(
                            gene, []
                        ).append(rna_transcript)
                        self.refseq_rna_version_to_gene_symbol[rna_transcript] = gene
                        if "." in rna_transcript:
                            rna_t = rna_transcript.split(".")[0]
                            self.refseq_rna_for_gene_symbol.setdefault(gene, []).append(
                                rna_t
                            )
                            self.refseq_rna_to_gene_symbol[rna_t] = gene
                    if refseq_transcript and rna_transcript:
                        self.np_to_nm[refseq_transcript] = rna_transcript

    def protein_transcripts(self, identifier: str) -> list[str]:
        """Return a list of protein transcripts for a gene symbol.

        >>> from cool_seq_tool.sources import TranscriptMappings
        >>> braf_txs = TranscriptMappings().protein_transcripts("BRAF")
        >>> braf_txs.sort()
        >>> braf_txs[-1]
        'NP_004324.2'

        :param identifier: Gene identifier to get protein transcripts for
        :return: Protein transcripts for a gene symbol
        """
        protein_transcripts = []
        protein_transcripts += self.ensembl_protein_version_for_gene_symbol.get(
            identifier, ""
        )
        protein_transcripts += self.ensembl_protein_for_gene_symbol.get(identifier, "")
        protein_transcripts += self.refseq_protein_for_gene_symbol.get(identifier, "")
        return list(set(protein_transcripts))

    def coding_dna_transcripts(self, identifier: str) -> list[str]:
        """Return transcripts from a coding dna refseq for a gene symbol.

        :param identifier: Gene identifier to find transcripts for
        :return: cDNA transcripts for a gene symbol
        """
        genomic_transcripts = []
        genomic_transcripts += self.ensembl_transcript_version_for_gene_symbol.get(
            identifier, ""
        )
        genomic_transcripts += self.refseq_rna_version_for_gene_symbol.get(
            identifier, ""
        )
        genomic_transcripts += self.refseq_rna_version_for_gene_symbol.get(
            identifier, ""
        )
        return list(set(genomic_transcripts))

    def get_gene_symbol_from_ensembl_protein(self, q: str) -> str | None:
        """Return the gene symbol for a Ensembl Protein.

        :param q: ensembl protein accession
        :return: Gene symbol
        """
        gene_symbol = self.ensembl_protein_version_to_gene_symbol.get(q)
        if not gene_symbol and "." in q:
            q = q.split(".")[0]
            gene_symbol = self.ensembl_protein_to_gene_symbol.get(q)
        return gene_symbol

    def get_gene_symbol_from_refeq_protein(self, q: str) -> str | None:
        """Return the gene symbol for a Refseq Protein.

        :param q: RefSeq protein accession
        :return: Gene symbol
        """
        return self.refseq_protein_to_gene_symbol.get(q)

    def get_gene_symbol_from_refseq_rna(self, q: str) -> str | None:
        """Return gene symbol for a Refseq RNA Transcript.

        :param q: RefSeq RNA transcript accession
        :return: Gene symbol
        """
        gene_symbol = self.refseq_rna_version_to_gene_symbol.get(q)
        if not gene_symbol and "." in q:
            q = q.split(".")[0]
            gene_symbol = self.refseq_rna_to_gene_symbol.get(q)
        return gene_symbol

    def get_gene_symbol_from_ensembl_transcript(self, q: str) -> str | None:
        """Return gene symbol for an Ensembl Transcript.

        :param q: Ensembl transcript accession
        :return: Gene symbol
        """
        gene_symbol = self.ensembl_transcript_version_to_gene_symbol.get(q)
        if not gene_symbol and "." in q:
            q = q.split(".")[0]
            gene_symbol = self.ensembl_transcript_to_gene_symbol.get(q)
        return gene_symbol
