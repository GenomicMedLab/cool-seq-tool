"""Provide transcript lookup and metadata tools via the UTA database.

In an asyncio runtime:

    >>> from cool_seq_tool.sources.uta_database import (
    ...     create_uta_connection_pool,
    ...     UtaDatabase,
    ...     UtaRepository,
    ... )
    >>> pool = await create_uta_connection_pool()
    >>> uta_db = UtaDatabase(pool)
    >>> async with uta_db.repository() as repo:
    ...     braf_exists = await repo.gene_exists("BRAF")
    >>> braf_exists
    True

"""

import ast
import logging
import os
import warnings
from collections.abc import AsyncIterator, Mapping, Sequence
from contextlib import asynccontextmanager
from typing import Literal
from urllib.parse import ParseResult as UrlLibParseResult
from urllib.parse import urlparse, urlunparse

import boto3
import polars as pl
from botocore.exceptions import ClientError
from psycopg import AsyncConnection, AsyncCursor
from psycopg_pool import AsyncConnectionPool
from pydantic import Field, StrictInt, StrictStr

from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    BaseModelForbidExtra,
    GenomicTxData,
    GenomicTxMetadata,
    Strand,
)

_logger = logging.getLogger(__name__)


class GenomicAlnData(BaseModelForbidExtra):
    """Represent genomic alignment data from UTA tx_exon_aln_mv view"""

    hgnc: StrictStr = Field(..., description="HGNC gene symbol.")
    ord: StrictInt = Field(..., description="Exon number. 0-based.")
    alt_ac: StrictStr = Field(..., description="RefSeq genomic accession.")
    alt_start_i: StrictInt = Field(
        ...,
        description="`alt_ac`'s start index of the exon using inter-residue coordinates.",
    )
    alt_end_i: StrictInt = Field(
        ...,
        description="`alt_ac`'s end index of the exon using inter-residue coordinates.",
    )
    alt_strand: Strand = Field(..., description="Strand.")


class TxExonAlnData(GenomicAlnData):
    """Represent data from UTA tx_exon_aln_mv view"""

    tx_ac: StrictStr = Field(..., description="Transcript accession.")
    tx_start_i: StrictInt = Field(
        ...,
        description="`tx_ac`'s start index of the exon using inter-residue coordinates.",
    )
    tx_end_i: StrictInt = Field(
        ...,
        description="`tx_ac`'s end index of the exon using inter-residue coordinates.",
    )
    alt_aln_method: StrictStr = Field(
        ..., description="The alignment method used to compare sequences."
    )
    tx_exon_id: StrictInt = Field(..., description="`tx_ac` exon identifier.")
    alt_exon_id: StrictInt = Field(..., description="`alt_ac` exon identifier.")


class NoMatchingAlignmentError(Exception):
    """Raise for failure to find alignment matching user parameters"""


class UtaRepository:
    """Connection-scoped repository for issuing queries against UTA.

    This class encapsulates predefined UTA queries and related result parsing.
    It operates on an active psycopg async connection provided at initialization
    time and does not manage connection lifecycle or pooling.

    Instances are intended to be short-lived and used within the scope of a
    checked-out connection (e.g., from a connection pool).
    """

    def __init__(self, conn: AsyncConnection) -> None:
        """Initialize the repository with an active database connection.

        :param conn: Active psycopg async connection to a UTA database.
            The caller is responsible for connection lifecycle management.
        """
        self._conn = conn

    async def execute_query(
        self, q: str, params: Sequence | Mapping | None = None
    ) -> AsyncCursor:
        """Execute an arbitrary query against the UTA DB

        This method is marked as public so that downstream applications can run custom
        queries using the same DB connection. However, that means they are responsible
        for managing the cursor themselves

        :param q: raw query. May need to specify schema depending on connection context.
        :param params: query variables, if needed. These should not be hard-coded into the query.
        :return: query result cursor
        """
        return await self._conn.execute(q, params)

    async def create_genomic_table(self) -> None:
        """Create table containing genomic accession information."""
        check_table_exists = """
            SELECT EXISTS (
               SELECT FROM information_schema.tables
               WHERE table_name = 'genomic'
            );
            """
        genomic_table_exists = await (
            await self.execute_query(check_table_exists)
        ).fetchone()
        if genomic_table_exists is None:
            _logger.critical(
                "SELECT EXISTS query in UtaDatabase._create_genomic_table "
                "returned invalid response"
            )
            msg = "SELECT EXISTS query returned invalid response"
            raise ValueError(msg)
        if not genomic_table_exists[0]:
            create_genomic_table = """
                CREATE TABLE genomic AS
                    SELECT t.hgnc, aes.alt_ac, aes.alt_aln_method,
                        aes.alt_strand, ae.start_i AS alt_start_i,
                        ae.end_i AS alt_end_i
                    FROM (((((transcript t
                        JOIN exon_set tes ON (((t.ac = tes.tx_ac)
                            AND (tes.alt_aln_method = 'transcript'::text))))
                        JOIN exon_set aes ON (((t.ac = aes.tx_ac)
                            AND (aes.alt_aln_method <> 'transcript'::text))))
                        JOIN exon te ON
                            ((tes.exon_set_id = te.exon_set_id)))
                        JOIN exon ae ON
                            (((aes.exon_set_id = ae.exon_set_id)
                            AND (te.ord = ae.ord))))
                        LEFT JOIN exon_aln ea ON
                            (((te.exon_id = ea.tx_exon_id) AND
                            (ae.exon_id = ea.alt_exon_id))));
                """
            await self.execute_query(create_genomic_table)

            indexes = [
                "CREATE INDEX alt_pos_index ON genomic (alt_ac, alt_start_i, alt_end_i);",
                "CREATE INDEX gene_alt_index ON genomic (hgnc, alt_ac);",
                "CREATE INDEX alt_ac_index ON genomic (alt_ac);",
            ]
            for create_index in indexes:
                await self.execute_query(create_index)

    async def get_alt_ac_start_or_end(
        self, tx_ac: str, tx_exon_start: int, tx_exon_end: int, gene: str | None
    ) -> GenomicAlnData:
        """Get genomic data for related transcript exon start or end.

        :param tx_ac: Transcript accession
        :param tx_exon_start: Transcript's exon start coordinate
        :param tx_exon_end: Transcript's exon end coordinate
        :param gene: HGNC gene symbol, if available
        :return: Genomic alignment data if match found
        :raise NoMatchingAlignmentError: if unable to find alignment matching given params
        """
        query = """
            SELECT
                T.hgnc,
                T.alt_ac,
                T.alt_start_i,
                T.alt_end_i,
                T.alt_strand,
                T.ord
            FROM _cds_exons_fp_v AS C
            JOIN tx_exon_aln_mv AS T ON T.tx_ac = C.tx_ac
            WHERE T.tx_ac = %(tx_ac)s
              AND (%(gene)s::text IS NULL OR T.hgnc = %(gene)s::text)
              AND %(tx_exon_start)s BETWEEN T.tx_start_i AND T.tx_end_i
              AND %(tx_exon_end)s BETWEEN T.tx_start_i AND T.tx_end_i
              AND T.alt_aln_method = 'splign'
              AND T.alt_ac LIKE 'NC_00%%'
            ORDER BY CAST(
                SUBSTR(
                    T.alt_ac,
                    POSITION('.' IN T.alt_ac) + 1,
                    LENGTH(T.alt_ac)
                ) AS INT
            ) DESC;
        """

        params = {
            "tx_ac": tx_ac,
            "tx_exon_start": tx_exon_start,
            "tx_exon_end": tx_exon_end,
            "gene": gene,
        }

        cur = await self.execute_query(query, params)
        row = await cur.fetchone()

        if not row:
            msg = (
                f"Unable to find a result where {tx_ac} has transcript "
                f"coordinates {tx_exon_start} and {tx_exon_end} between "
                f"an exon's start and end coordinates"
            )
            _logger.warning(msg)
            raise NoMatchingAlignmentError(msg)

        return GenomicAlnData(
            hgnc=row[0],
            alt_ac=row[1],
            alt_start_i=row[2],
            alt_end_i=row[3],
            alt_strand=row[4],
            ord=row[5],
        )

    async def get_cds_start_end(self, tx_ac: str) -> tuple[int, int] | None:
        """Return CDS start/end coordinates for a transcript.

        Strips version from Ensembl accessions (``ENS*``) since UTA stores them
        unversioned.

        :param tx_ac: Transcript accession
        :return: (cds_start_i, cds_end_i) if both exist, else None
        """
        # As of 2026-03, Ensembl transcripts in UTA are unversioned, so we need to drop
        # the version specifier
        tx_ac = tx_ac.split(".", 1)[0] if tx_ac.startswith("ENS") else tx_ac
        query = """
            SELECT cds_start_i, cds_end_i
            FROM transcript
            WHERE ac=%(ac)s;
        """
        cds_start_end = await (
            await self.execute_query(query, {"ac": tx_ac})
        ).fetchone()
        if cds_start_end:
            if cds_start_end[0] is not None and cds_start_end[1] is not None:
                return cds_start_end[0], cds_start_end[1]
        else:
            _logger.warning(
                "Unable to get coding start/end site for accession: %s", tx_ac
            )
        return None

    async def get_newest_assembly_ac(self, ac: str) -> list[str]:
        """Return newest accession versions matching the given prefix

        If the accession is Ensembl (``EN`` prefix), results are ordered lexicographically.
        Otherwise, RefSeq-style accessions are ordered by version number in descending order.

        :param ac: Accession (versioned or unversioned)
        :return: List of matching accessions, newest version first
        """
        prefix = ac.split(".", 1)[0]

        if ac.startswith("EN"):
            query = """
                SELECT ac
                FROM _seq_anno_most_recent
                WHERE ac LIKE %(ac_prefix)s
                  AND (descr IS NULL OR descr = '')
                ORDER BY ac;
            """
        else:
            query = """
                SELECT ac
                FROM _seq_anno_most_recent
                WHERE ac LIKE %(ac_prefix)s
                  AND (descr IS NULL OR descr = '')
                ORDER BY
                    SUBSTR(ac, 0, POSITION('.' IN ac)),
                    CAST(SUBSTR(ac, POSITION('.' IN ac) + 1, LENGTH(ac)) AS INT) DESC;
            """

        params = {
            "ac_prefix": f"{prefix}%",
        }
        results = await (await self.execute_query(query, params)).fetchall()
        return [r[0] for r in results]

    async def validate_genomic_ac(self, ac: str) -> bool:
        """Return whether or not genomic accession exists.

        :param ac: Genomic accession
        :return: ``True`` if genomic accession exists. ``False`` otherwise.
        """
        query = """
            SELECT EXISTS(
                SELECT ac
                FROM _seq_anno_most_recent
                WHERE ac = %(ac)s
            );
            """
        cursor = await self.execute_query(query, {"ac": ac})
        result = await cursor.fetchone()
        return result[0]

    ### working here

    async def gene_exists(self, gene: str) -> bool:
        """Return whether or not a gene symbol exists in UTA gene table

        :param gene: Gene symbol
        :return ``True`` if gene symbol exists in UTA, ``False`` if not
        """
        query = """
            SELECT EXISTS(
                SELECT hgnc
                FROM gene
                WHERE hgnc = %(gene)s
            );
            """
        cursor = await self.execute_query(query, {"gene": gene})
        result = await cursor.fetchone()
        return result[0]

    async def transcript_exists(self, transcript: str) -> bool:
        """Return whether or not a transcript exists in the UTA ``tx_exon_aln_mv`` table

        :param transcript: A transcript accession
        :return: ``True`` if transcript exists in UTA, ``False`` if not
        """
        query = """
            SELECT EXISTS(
                SELECT tx_ac
                FROM tx_exon_aln_mv
                WHERE tx_ac = %(tx_ac)s
            );
            """
        cursor = await self.execute_query(query, {"tx_ac": transcript})
        result = await cursor.fetchone()
        return result[0]

    async def get_ac_descr(self, ac: str) -> str | None:
        """Return free-text accession description

        This is typically available only for accessions from older (pre-GRCh38) builds.

        >>> async with uta.repository() as repo:
        ...     result = await repo.get_ac_descr("NC_000001.10")
        >>> result
        'Homo sapiens chromosome 1, GRCh37.p13 Primary Assembly'

        :param ac: chromosome accession, e.g. ``"NC_000001.10"``
        :return: Free-text description provided by source, generally containing assembly and chromosome
        """
        query = """
            SELECT descr
            FROM _seq_anno_most_recent
            WHERE ac = %(ac)s;
            """
        cursor = await self.execute_query(query, {"ac": ac})
        result = await cursor.fetchone()
        if not result:
            _logger.warning("No description entry found for accession %s", ac)
            return None
        result = result[0]
        if result == "":
            result = None
        return result

    async def get_tx_exon_aln_data(
        self,
        tx_ac: str,
        start_pos: int,
        end_pos: int,
        alt_ac: str | None = None,
        use_tx_pos: bool = True,
        like_tx_ac: bool = False,
    ) -> list[TxExonAlnData]:
        """Get alignments between exons and reference sequences.

        This is a direct query against the UTA ``tx_exon_aln_mv`` view.

        :param tx_ac: accession on c. coordinate
        :param start_pos: Start position change
        :param end_pos: End position change
        :param alt_ac: accession on g. coordinate
        :param use_tx_pos: ``True`` if querying on transcript position. This means
            ``start_pos`` and ``end_pos`` are on the c. coordinate
            ``False`` if querying on genomic position. This means ``start_pos`` and
            ``end_pos`` are on the g. coordinate
        :param like_tx_ac: ``True`` if tx_ac condition should be a like statement.
            This is used when you want to query an accession regardless of its version
            ``False`` if tx_condition will be exact match
        :return: List of transcript exon alignment data
        """
        params: dict = {"start_pos": start_pos, "end_pos": end_pos}
        if tx_ac.startswith("EN"):
            params["tx_ac"] = tx_ac.split(".")[0]
            params["alt_aln_method"] = "genebuild"
        else:
            params["tx_ac"] = tx_ac
            params["alt_aln_method"] = "splign"

        if like_tx_ac:
            params["tx_ac"] = f"{params['tx_ac']}%"
            tx_q = "WHERE tx_ac LIKE %(tx_ac)s"
        else:
            tx_q = "WHERE tx_ac=%(tx_ac)s"

        order_by_cond = "ORDER BY CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1, LENGTH(alt_ac)) AS INT)"
        if alt_ac:
            alt_ac_q = "AND alt_ac = %(alt_ac)s"
            params["alt_ac"] = alt_ac
            if alt_ac.startswith("EN"):
                order_by_cond = "ORDER BY alt_ac"
        else:
            alt_ac_q = "AND alt_ac LIKE 'NC_00%%'"

        if use_tx_pos:
            pos_q = """tx_start_i AND tx_end_i"""
        else:
            pos_q = """alt_start_i AND alt_end_i"""

        query = f"""
            SELECT hgnc, tx_ac, tx_start_i, tx_end_i, alt_ac, alt_start_i,
                alt_end_i, alt_strand, alt_aln_method, ord, tx_exon_id, alt_exon_id
            FROM tx_exon_aln_mv
            {tx_q}
            {alt_ac_q}
            AND alt_aln_method = %(alt_aln_method)s
            AND %(start_pos)s BETWEEN {pos_q}
            AND %(end_pos)s BETWEEN {pos_q}
            {order_by_cond}
            """  # noqa: S608
        cursor = await self.execute_query(query, params)
        results = await cursor.fetchall()
        if not results:
            _logger.warning("Unable to find transcript alignment for query: %s", query)
            return []
        if alt_ac and not use_tx_pos and len(results) > 1:
            _logger.debug(
                "Found more than one match for tx_ac %s and alt_ac = %s",
                params["tx_ac"],
                alt_ac,
            )
        return [
            TxExonAlnData(
                hgnc=r[0],
                tx_ac=r[1],
                tx_start_i=r[2],
                tx_end_i=r[3],
                alt_ac=r[4],
                alt_start_i=r[5],
                alt_end_i=r[6],
                alt_strand=r[7],
                alt_aln_method=r[8],
                ord=r[9],
                tx_exon_id=r[10],
                alt_exon_id=r[11],
            )
            for r in results
        ]

    @staticmethod
    def data_from_result(result: TxExonAlnData) -> GenomicTxData | None:
        """Return data found from result.

        :param result: Transcript exon alignment data
        :return: Aligned genomic / transcript exon data
        """
        tx_pos_range = result.tx_start_i, result.tx_end_i
        alt_pos_range = result.alt_start_i, result.alt_end_i

        if (tx_pos_range[1] - tx_pos_range[0]) != (alt_pos_range[1] - alt_pos_range[0]):
            _logger.warning(
                "tx_pos_range %s is not the same length as alt_pos_range %s.",
                tx_pos_range,
                alt_pos_range,
            )
            return None

        return GenomicTxData(
            gene=result.hgnc,
            strand=Strand(result.alt_strand),
            tx_pos_range=tx_pos_range,
            alt_pos_range=alt_pos_range,
            alt_aln_method=result.alt_aln_method,
            tx_exon_id=result.tx_exon_id,
            alt_exon_id=result.alt_exon_id,
        )

    async def get_mane_c_genomic_data(
        self, ac: str, alt_ac: str | None, start_pos: int, end_pos: int
    ) -> GenomicTxMetadata | None:
        """Get MANE transcript and genomic data. Used when going from g. to MANE c.
        representation. This function parses queried data from the tx_exon_aln_mv
        table, and sorts the queried data by the most recent genomic build

        >>> async with uta_db.repository() as repo:
        ...     result = await repo.get_mane_c_genomic_data(
        ...         "NM_004333.6",
        ...         None,
        ...         140753335,
        ...         140753335,
        ...     )
        >>> result.alt_ac
        'NC_000007.14'

        :param ac: MANE transcript accession
        :param alt_ac: NC accession. Used to triangulate on correct genomic data. Can
            be set to ``None`` if unavailable.
        :param start_pos: Genomic start position
        :param end_pos: Genomic end position change
        :return: Metadata for MANE genomic and transcript accessions results if
            successful
        """
        results = await self.get_tx_exon_aln_data(
            tx_ac=ac,
            start_pos=start_pos,
            end_pos=end_pos,
            alt_ac=alt_ac,
            use_tx_pos=False,
        )
        if not results:
            return None

        # Sort by most recent chromosomal accession
        result = results[-1]

        genomic_tx_data = self.data_from_result(result)
        if not genomic_tx_data:
            return None

        coding_start_site = await self.get_cds_start_end(ac)
        if coding_start_site is None:
            _logger.warning("Accession %s not found in UTA", ac)
            return None

        coding_start_site, coding_end_site = coding_start_site

        if genomic_tx_data.strand == Strand.NEGATIVE:
            alt_pos_change_range = (end_pos, start_pos)
            pos_change = (
                genomic_tx_data.alt_pos_range[1] - alt_pos_change_range[0],
                alt_pos_change_range[1] - genomic_tx_data.alt_pos_range[0],
            )
        else:
            alt_pos_change_range = (start_pos, end_pos)
            pos_change = (
                alt_pos_change_range[0] - genomic_tx_data.alt_pos_range[0],
                genomic_tx_data.alt_pos_range[1] - alt_pos_change_range[1],
            )

        return GenomicTxMetadata(
            **genomic_tx_data.model_dump(),
            pos_change=pos_change,
            tx_ac=result.tx_ac,
            alt_ac=result.alt_ac,
            coding_start_site=coding_start_site,
            coding_end_site=coding_end_site,
            alt_pos_change_range=alt_pos_change_range,
        )

    async def get_genomic_tx_data(
        self,
        tx_ac: str,
        pos: tuple[int, int],
        annotation_layer: Literal[AnnotationLayer.CDNA]
        | Literal[AnnotationLayer.GENOMIC] = AnnotationLayer.CDNA,
        alt_ac: str | None = None,
        target_genome_assembly: Assembly = Assembly.GRCH38,
    ) -> GenomicTxMetadata | None:
        """Get transcript mapping to genomic data.

        :param tx_ac: Accession on c. coordinate
        :param pos: (start pos, end pos). These must describe the inter-residue
            coordinates that are being examined.
        :param annotation_layer: Annotation layer for ``ac`` and ``pos``
        :param alt_ac: Accession on g. coordinate
        :param target_genome_assembly: Genome assembly to get genomic data for.
            If ``alt_ac`` is provided, it will return the associated assembly.
        :return: Metadata for genomic and transcript accessions
        """
        results = await self.get_tx_exon_aln_data(
            tx_ac,
            pos[0],
            pos[1],
            use_tx_pos=annotation_layer == AnnotationLayer.CDNA,
            alt_ac=alt_ac,
        )
        if not results:
            return None

        if alt_ac or target_genome_assembly == Assembly.GRCH38:
            result = results[-1]
        else:
            result = results[0]

        genomic_tx_data = self.data_from_result(result)
        if not genomic_tx_data:
            return None

        pos_change = (
            pos[0] - genomic_tx_data.tx_pos_range[0],
            genomic_tx_data.tx_pos_range[1] - pos[1],
        )

        if annotation_layer == AnnotationLayer.CDNA:
            if genomic_tx_data.strand == Strand.NEGATIVE:
                alt_pos_change_range = (
                    genomic_tx_data.alt_pos_range[1] - pos_change[0],
                    genomic_tx_data.alt_pos_range[0] + pos_change[1],
                )
            else:
                alt_pos_change_range = (
                    genomic_tx_data.alt_pos_range[0] + pos_change[0],
                    genomic_tx_data.alt_pos_range[1] - pos_change[1],
                )
        elif genomic_tx_data.strand == Strand.NEGATIVE:
            alt_pos_change_range = (pos[1], pos[0])
        else:
            alt_pos_change_range = pos

        return GenomicTxMetadata(
            **genomic_tx_data.model_dump(),
            tx_ac=result.tx_ac,
            alt_ac=result.alt_ac,
            pos_change=pos_change,
            alt_pos_change_range=alt_pos_change_range,
        )

    async def get_ac_from_gene(self, gene: str) -> list[str]:
        """Return genomic accession(s) associated to a gene.

        :param gene: Gene symbol
        :return: List of genomic accessions, sorted in desc order
        """
        query = """
            SELECT DISTINCT alt_ac
            FROM genomic
            WHERE hgnc = %(gene)s
            AND alt_ac LIKE 'NC_00%%'
            ORDER BY alt_ac;
            """

        cursor = await self.execute_query(query, {"gene": gene})
        results = await cursor.fetchall()
        alt_acs = [r[0] for r in results]
        alt_acs.sort(key=lambda x: int(x.split(".")[-1]), reverse=True)
        return alt_acs

    async def get_gene_from_ac(
        self, ac: str, start_pos: int, end_pos: int | None
    ) -> list[str] | None:
        """Get gene(s) within the provided coordinate range

        >>> async with uta_db.repository() as repo:
        ...     result = await repo.get_gene_from_ac("NC_000017.11", 43044296, 43045802)
        >>> result
        ['BRCA1']

        :param ac: NC accession, e.g. ``"NC_000001.11"``
        :param start_pos: Start position change
        :param end_pos: End position change
        :return: List of HGNC gene symbols
        """
        if end_pos is None:
            end_pos = start_pos
        query = """
            SELECT DISTINCT hgnc
            FROM genomic
            WHERE alt_ac = %(ac)s
            AND %(start_pos)s BETWEEN alt_start_i AND alt_end_i
            AND %(end_pos)s BETWEEN alt_start_i AND alt_end_i;
        """
        cursor = await self.execute_query(
            query, {"ac": ac, "start_pos": start_pos, "end_pos": end_pos}
        )
        results = await cursor.fetchall()
        if not results:
            _logger.warning(
                "Unable to find gene between %s and %s on %s", start_pos, end_pos, ac
            )
            return None
        if len(results) > 1:
            _logger.info(
                "Found more than one gene between %s and %s on %s",
                start_pos,
                end_pos,
                ac,
            )

        return [r[0] for r in results]

    async def get_transcripts(
        self,
        start_pos: int | None = None,
        end_pos: int | None = None,
        gene: str | None = None,
        use_tx_pos: bool = True,
        alt_ac: str | None = None,
    ) -> pl.DataFrame:
        """Get transcripts for a given ``gene`` or ``alt_ac`` related to optional positions.

        :param start_pos: Start position change
            If not provided and ``end_pos`` not provided, all transcripts associated with
            the gene and/or accession will be returned
        :param end_pos: End position change
            If not provided and ``start_pos`` not provided, all transcripts associated
            with the gene and/or accession will be returned
        :param gene: HGNC gene symbol
        :param use_tx_pos: ``True`` if querying on transcript position. This means
            ``start_pos`` and ``end_pos`` are c. coordinate positions. ``False`` if querying
            on genomic position. This means ``start_pos`` and ``end_pos`` are g. coordinate
            positions
        :param alt_ac: Genomic accession. If not provided, must provide ``gene``
        :return: Data Frame containing transcripts associated with a gene.
            Transcripts are ordered by most recent NC accession, then by
            descending transcript length
        """
        schema = ["pro_ac", "tx_ac", "alt_ac", "cds_start_i"]
        if not gene and not alt_ac:
            return pl.DataFrame([], schema=schema)

        pos_cond = ""
        if start_pos is not None and end_pos is not None:
            if use_tx_pos:
                pos_cond = """
                    AND %(start_pos)s + T.cds_start_i
                        BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                    AND %(end_pos)s + T.cds_start_i
                        BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                    """
            else:
                pos_cond = """
                    AND %(start_pos)s BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                    AND %(end_pos)s BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                    """

        order_by_cond = """
        ORDER BY SUBSTR(ALIGN.alt_ac, 0, position('.' in ALIGN.alt_ac)),
        CAST(SUBSTR(ALIGN.alt_ac, position('.' in ALIGN.alt_ac) + 1,
            LENGTH(ALIGN.alt_ac)) AS INT) DESC,
        ALIGN.tx_end_i - ALIGN.tx_start_i DESC;
        """
        if alt_ac:
            alt_ac_cond = "AND ALIGN.alt_ac = %(alt_ac)s"
            if alt_ac.startswith("EN"):
                order_by_cond = "ORDER BY ALIGN.alt_ac;"
        else:
            alt_ac_cond = "AND ALIGN.alt_ac LIKE 'NC_00%%'"

        gene_cond = "AND T.hgnc = %(gene)s" if gene else ""

        query = f"""
            SELECT AA.pro_ac, AA.tx_ac, ALIGN.alt_ac, T.cds_start_i
            FROM associated_accessions as AA
            JOIN transcript as T ON T.ac = AA.tx_ac
            JOIN tx_exon_aln_mv as ALIGN ON T.ac = ALIGN.tx_ac
            WHERE ALIGN.alt_aln_method = 'splign'
            {gene_cond}
            {alt_ac_cond}
            {pos_cond}
            {order_by_cond}
            """  # noqa: S608
        cursor = await self.execute_query(
            query,
            {
                "start_pos": start_pos,
                "end_pos": end_pos,
                "gene": gene,
                "alt_ac": alt_ac,
            },
        )
        results = await cursor.fetchall()
        results = [(r[0], r[1], r[2], r[3]) for r in results]
        results_df = pl.DataFrame(results, schema=schema, orient="row")
        if results:
            results_df = results_df.unique()
        return results_df

    async def get_chr_assembly(self, ac: str) -> tuple[str, Assembly] | None:
        """Get chromosome and assembly for NC accession if not in GRCh38.

        >>> async with uta_db.repository() as repo:
        ...     result = await repo.get_chr_assembly("NC_000007.13")
        >>> result
        ('chr7', <Assembly.GRCH37: 'GRCh37'>)

        Returns ``None`` if unable to find (either unrecognized/invalid, or
        a GRCh38 accession).

        :param ac: RefSeq NC accession, eg ``"NC_000007.13"``
        :return: Chromosome and assembly that accession is on, if available.
        """
        descr = await self.get_ac_descr(ac)
        if not descr:
            # Already GRCh38 Assembly
            return None
        descr = descr.split(",")
        chromosome = f"chr{descr[0].split()[-1]}"
        assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

        try:
            assembly = Assembly(assembly)
        except ValueError:
            _logger.exception("Unable to parse %s as an Assembly", assembly)
            return None

        return chromosome, assembly

    async def p_to_c_ac(self, p_ac: str) -> list[str]:
        """Return cDNA reference sequence accession from protein reference sequence
        accession (i.e. ``p.`` to ``c.`` in HGVS syntax)

        :param p_ac: Protein accession
        :return: List of rows containing c. accessions that are associated with the
            given p. accession. In ascending order.
        """
        # Ensembl accessions do not have versions
        if p_ac.startswith("EN"):
            order_by_cond = "ORDER BY tx_ac;"
        else:
            order_by_cond = """
            ORDER BY SUBSTR(tx_ac, 0, position('.' in tx_ac)),
            CAST(SUBSTR(tx_ac, position('.' in tx_ac) + 1, LENGTH(tx_ac)) AS INT);
            """

        query = f"""
            SELECT tx_ac
            FROM associated_accessions
            WHERE pro_ac = %(p_ac)s
            {order_by_cond}
            """  # noqa: S608
        cursor = await self.execute_query(query, {"p_ac": p_ac})
        result = await cursor.fetchall()
        return [r[0] for r in result]

    async def get_transcripts_from_genomic_pos(
        self, alt_ac: str, g_pos: int
    ) -> list[str]:
        """Get transcripts associated to a genomic ac and position.

        :param alt_ac: Genomic accession
        :param g_pos: Genomic position
        :return: RefSeq transcripts on c. coordinate
        """
        query = """
           SELECT distinct tx_ac
           FROM tx_exon_aln_mv
           WHERE alt_ac = %(alt_ac)s
           AND %(g_pos)s BETWEEN alt_start_i AND alt_end_i
           AND tx_ac LIKE 'NM_%%';
       """
        cursor = await self.execute_query(query, {"alt_ac": alt_ac, "g_pos": g_pos})
        results = await cursor.fetchall()
        return [item for sublist in results for item in sublist]


class ParseResult(UrlLibParseResult):
    """Subclass of url.ParseResult that adds database and schema methods, and provides stringification.

    Source: https://github.com/biocommons/hgvs
    """

    def __new__(cls, pr):  # noqa: ANN001, ANN204
        """Create new instance."""
        return super(ParseResult, cls).__new__(cls, *pr)  # noqa: UP008

    @property
    def database(self) -> str | None:
        """Create database property."""
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self) -> str | None:
        """Create schema property."""
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None  # noqa: PLR2004

    @property
    def sanitized_url(self) -> str:
        """Sanitized DB URL with the password masked"""
        netloc = ""
        if self.username:
            netloc += self.username
            if self.password is not None and self.password != "":
                netloc += ":***"
            netloc += "@"
        if self.hostname:
            netloc += f"{self.hostname}"
        if self.port:
            netloc += f":{self.port}"

        return urlunparse(
            (
                self.scheme,
                netloc,
                self.path,
                self.params,
                self.query,
                self.fragment,
            )
        )


def _get_secret_args() -> tuple[str, str]:
    """Get secrets connection args for UTA DB instances. Used for deployment on AWS.

    This function is tightly coupled to our internal deployment policies;
    it is subject to change or removal in the future.

    :raises ClientError: If unable to retrieve secret value due to decryption
        decryption failure, internal service error, invalid parameter, invalid
        request, or resource not found.
    """
    warnings.warn(
        "Deprecated; subject to change in future releases",
        DeprecationWarning,
        stacklevel=2,
    )

    secret_name = os.environ["UTA_DB_SECRET"]
    region_name = "us-east-2"

    session = boto3.session.Session()
    client = session.client(service_name="secretsmanager", region_name=region_name)
    try:
        get_secret_value_response = client.get_secret_value(SecretId=secret_name)
    except ClientError:
        # For a list of exceptions thrown, see
        # https://docs.aws.amazon.com/secretsmanager/latest/apireference/API_GetSecretValue.html
        _logger.exception("Encountered AWS client error fetching UTA DB secret")
        raise
    secret_val = get_secret_value_response["SecretString"]
    secret = ast.literal_eval(secret_val)

    username, password = secret["username"], secret["password"]
    port, host = secret["port"], secret["host"]
    database = secret["dbname"]
    schema = secret["schema"]
    db_uri = f"postgresql://{username}{':' + password if password else ''}@{host}:{port}/{database}"

    return db_uri, schema


DEFAULT_UTA_DB_URI = "postgresql://uta_admin@localhost:5432/uta"
DEFAULT_UTA_SCHEMA = "uta_20241220"


async def create_uta_connection_pool(
    db_uri: str | None = None,
    uta_schema: str | None = None,
) -> AsyncConnectionPool:
    """Create and initialize a UTA connection pool.

    Connection parameters are resolved in the following order:

    1. If the ``UTA_DB_PROD`` environment variable is set, credentials and schema
       are retrieved from a secret manager via ``_get_secret_args()``.
    2. Otherwise, explicitly provided ``db_uri`` and ``uta_schema`` arguments are used.
    3. If not provided, fall back to environment variables (``UTA_DB_URI``,
       ``UTA_DB_SCHEMA``), then to default constants.

    The resulting connection string is augmented with a ``search_path`` option so that
    all connections automatically target the specified UTA schema.

    After opening the pool, a one-time initialization step is performed to ensure that
    required genomic tables are present.

    :param db_uri: PostgreSQL connection URI (e.g., ``postgresql://user@host:port/db``).
        If not provided, resolved from environment or defaults.
    :param uta_schema: Target UTA schema name (e.g., ``uta_20241220``). Used to set
        the PostgreSQL ``search_path``. If not provided, resolved from environment
        or defaults.
    :return: An open ``AsyncConnectionPool`` configured for the UTA database
    """
    if "UTA_DB_PROD" in os.environ:
        db_uri, uta_schema = _get_secret_args()
    else:
        if db_uri is None:
            db_uri = os.environ.get("UTA_DB_URI", DEFAULT_UTA_DB_URI)
        if uta_schema is None:
            uta_schema = os.environ.get("UTA_DB_SCHEMA", DEFAULT_UTA_SCHEMA)
    _logger.info(
        "Creating connection pool with db_uri=%s and schema=%s",
        ParseResult(urlparse(db_uri)).sanitized_url,
        uta_schema,
    )
    dsn = f"{db_uri}?options=-csearch_path%3D{uta_schema},public"
    pool = AsyncConnectionPool(conninfo=dsn, open=False)
    await pool.open()
    async with pool.connection() as conn:
        await UtaRepository(conn).create_genomic_table()
    return pool


class UtaDatabase:
    """Provide pooled access to connection-scoped UTA repositories.

    This class owns or borrows an async psycopg connection pool and yields
    ``UtaRepository`` instances bound to checked-out connections.
    """

    def __init__(self, pool: AsyncConnectionPool | None = None) -> None:
        """Initialize access wrapper.

        :param pool: Existing async connection pool to use. If omitted, a default
            pool is created lazily on first use.
        """
        self._connection_pool = pool

    async def open(self) -> None:
        """Initialize connection"""
        if self._connection_pool is None:
            self._connection_pool = await create_uta_connection_pool()

    @asynccontextmanager
    async def repository(self) -> AsyncIterator[UtaRepository]:
        """Yield a ``UtaRepository`` backed by a pooled connection.

        If no pool has been provided yet, a default one is created on first use.

        :yield: Repository bound to an active pooled connection
        """
        await self.open()

        async with self._connection_pool.connection() as conn:
            yield UtaRepository(conn)

    async def close(self) -> None:
        """Close the owned connection pool, if present."""
        if self._connection_pool is None:
            _logger.info("Attempted to close nonexistent UTA access connection pool")
            return

        await self._connection_pool.close()
        self._connection_pool = None
