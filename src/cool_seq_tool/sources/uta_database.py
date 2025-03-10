"""Provide transcript lookup and metadata tools via the UTA database."""

import ast
import logging
from os import environ
from typing import Any, Literal, TypeVar
from urllib.parse import ParseResult as UrlLibParseResult
from urllib.parse import quote, unquote, urlparse

import asyncpg
import boto3
import polars as pl
from asyncpg.exceptions import InterfaceError, InvalidAuthorizationSpecificationError
from botocore.exceptions import ClientError
from pydantic import Field, StrictInt, StrictStr

from cool_seq_tool.schemas import (
    AnnotationLayer,
    Assembly,
    BaseModelForbidExtra,
    GenomicTxData,
    GenomicTxMetadata,
    Strand,
)

# use `bound` to upper-bound UtaDatabase or child classes
UTADatabaseType = TypeVar("UTADatabaseType", bound="UtaDatabase")

UTA_DB_URL = environ.get(
    "UTA_DB_URL", "postgresql://uta_admin:uta@localhost:5432/uta/uta_20241220"
)

_logger = logging.getLogger(__name__)


class DbConnectionArgs(BaseModelForbidExtra):
    """Represent database connection arguments"""

    host: str
    port: int
    user: str
    password: str
    database: str


class GenomicAlnData(BaseModelForbidExtra):
    """Represent genomic alignment data from UTA tx_exon_aln_v view"""

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
    """Represent data from UTA tx_exon_aln_v view"""

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


class UtaDatabase:
    """Provide transcript lookup and metadata tools via the Universal Transcript Archive
    (UTA) database.

    Users should use the ``create()`` method to construct a new instance. Note that
    almost all public methods are defined as ``async`` -- see the :ref:`Usage section <async_note>`
    for more information.

    >>> import asyncio
    >>> from cool_seq_tool.sources.uta_database import UtaDatabase
    >>> uta_db = asyncio.run(UtaDatabase.create())
    """

    def __init__(self, db_url: str = UTA_DB_URL) -> None:
        """Initialize DB class. Should only be used by ``create()`` method, and not
        be called directly by a user.

        :param db_url: PostgreSQL connection URL
            Format: ``driver://user:password@host/database/schema``
        """
        self.schema = None
        self._connection_pool = None
        original_pwd = db_url.split("//")[-1].split("@")[0].split(":")[-1]
        self.db_url = db_url.replace(original_pwd, quote(original_pwd))
        self.args = self._get_conn_args()

    def _get_conn_args(self) -> DbConnectionArgs:
        """Return connection arguments.

        :param db_url: raw connection URL
        :return: Database connection arguments
        """
        if "UTA_DB_PROD" in environ:
            secret = ast.literal_eval(self.get_secret())

            password = secret["password"]
            username = secret["username"]
            port = secret["port"]
            host = secret["host"]
            database = secret["dbname"]
            schema = secret["schema"]
            self.schema = schema

            environ["PGPASSWORD"] = password
            environ["UTA_DB_URL"] = (
                f"postgresql://{username}@{host}:{port}/{database}/{schema}"
            )
            return DbConnectionArgs(
                host=host,
                port=int(port),
                database=database,
                user=username,
                password=password,
            )

        url = ParseResult(urlparse(self.db_url))
        self.schema = url.schema
        password = unquote(url.password) if url.password else ""
        return DbConnectionArgs(
            host=url.hostname,
            port=url.port,
            database=url.database,
            user=url.username,
            password=password,
        )

    async def create_pool(self) -> None:
        """Create connection pool if not already created."""
        if not self._connection_pool:
            self.args = self._get_conn_args()
            try:
                self._connection_pool = await asyncpg.create_pool(
                    min_size=1,
                    max_size=10,
                    max_inactive_connection_lifetime=3,
                    command_timeout=60,
                    host=self.args.host,
                    port=self.args.port,
                    user=self.args.user,
                    password=self.args.password,
                    database=self.args.database,
                )
            except InterfaceError as e:
                _logger.error(
                    "While creating connection pool, encountered exception %s", e
                )
                msg = "Could not create connection pool"
                raise Exception(msg) from e

    @classmethod
    async def create(
        cls: type[UTADatabaseType], db_url: str = UTA_DB_URL
    ) -> UTADatabaseType:
        """Manufacture a fully-initialized class instance (a la factory pattern). This
        method should be used instead of calling the class directly to create a new
        instance.

        >>> import asyncio
        >>> from cool_seq_tool.sources.uta_database import UtaDatabase
        >>> uta_db = asyncio.run(UtaDatabase.create())

        :param cls: supplied implicitly
        :param db_url: PostgreSQL connection URL
            Format: ``driver://user:password@host/database/schema``
        :return: UTA DB access class instance
        """
        self = cls(db_url)
        await self._create_genomic_table()
        await self.create_pool()
        return self

    async def execute_query(self, query: str) -> Any:  # noqa: ANN401
        """Execute a query and return its result.

        :param query: Query to make on database
        :return: Query's result
        """

        async def _execute_query(q: str) -> Any:  # noqa: ANN401
            async with (
                self._connection_pool.acquire() as connection,
                connection.transaction(),
            ):
                return await connection.fetch(q)

        if not self._connection_pool:
            await self.create_pool()
        try:
            return await _execute_query(query)
        except InvalidAuthorizationSpecificationError:
            self._connection_pool = None
            await self.create_pool()
            return await _execute_query(query)

    async def _create_genomic_table(self) -> None:
        """Create table containing genomic accession information."""
        check_table_exists = f"""
            SELECT EXISTS (
               SELECT FROM information_schema.tables
               WHERE table_schema = '{self.schema}'
               AND table_name = 'genomic'
            );
            """  # noqa: S608
        genomic_table_exists = await self.execute_query(check_table_exists)
        genomic_table_exists = genomic_table_exists[0].get("exists")
        if genomic_table_exists is None:
            _logger.critical(
                "SELECT EXISTS query in UtaDatabase._create_genomic_table "
                "returned invalid response"
            )
            msg = "SELECT EXISTS query returned invalid response"
            raise ValueError(msg)
        if not genomic_table_exists:
            create_genomic_table = f"""
                CREATE TABLE {self.schema}.genomic AS
                    SELECT t.hgnc, aes.alt_ac, aes.alt_aln_method,
                        aes.alt_strand, ae.start_i AS alt_start_i,
                        ae.end_i AS alt_end_i
                    FROM ((((({self.schema}.transcript t
                        JOIN {self.schema}.exon_set tes ON (((t.ac = tes.tx_ac)
                            AND (tes.alt_aln_method = 'transcript'::text))))
                        JOIN {self.schema}.exon_set aes ON (((t.ac = aes.tx_ac)
                            AND (aes.alt_aln_method <> 'transcript'::text))))
                        JOIN {self.schema}.exon te ON
                            ((tes.exon_set_id = te.exon_set_id)))
                        JOIN {self.schema}.exon ae ON
                            (((aes.exon_set_id = ae.exon_set_id)
                            AND (te.ord = ae.ord))))
                        LEFT JOIN {self.schema}.exon_aln ea ON
                            (((te.exon_id = ea.tx_exon_id) AND
                            (ae.exon_id = ea.alt_exon_id))));
                """  # noqa: S608
            await self.execute_query(create_genomic_table)

            indexes = [
                f"""CREATE INDEX alt_pos_index ON {self.schema}.genomic (alt_ac, alt_start_i, alt_end_i);""",
                f"""CREATE INDEX gene_alt_index ON {self.schema}.genomic (hgnc, alt_ac);""",
                f"""CREATE INDEX alt_ac_index ON {self.schema}.genomic (alt_ac);""",
            ]
            for create_index in indexes:
                await self.execute_query(create_index)

    @staticmethod
    def _transform_list(li: list) -> list[list[Any]]:
        """Transform list to only contain field values

        :param li: List of asyncpg.Record objects
        :return: List of list of objects
        """
        return [list(i) for i in li]

    async def get_alt_ac_start_or_end(
        self, tx_ac: str, tx_exon_start: int, tx_exon_end: int, gene: str | None
    ) -> tuple[GenomicAlnData | None, str | None]:
        """Get genomic data for related transcript exon start or end.

        :param tx_ac: Transcript accession
        :param tx_exon_start: Transcript's exon start coordinate
        :param tx_exon_end: Transcript's exon end coordinate
        :param gene: HGNC gene symbol
        :return: Genomic alignment data and warnings if found
        """
        gene_query = f"AND T.hgnc = '{gene}'" if gene else ""

        query = f"""
            SELECT T.hgnc, T.alt_ac, T.alt_start_i, T.alt_end_i, T.alt_strand, T.ord
            FROM {self.schema}._cds_exons_fp_v as C
            JOIN {self.schema}.tx_exon_aln_v as T ON T.tx_ac = C.tx_ac
            WHERE T.tx_ac = '{tx_ac}'
            {gene_query}
            AND {tx_exon_start} BETWEEN T.tx_start_i AND T.tx_end_i
            AND {tx_exon_end} BETWEEN T.tx_start_i AND T.tx_end_i
            AND T.alt_aln_method = 'splign'
            AND T.alt_ac LIKE 'NC_00%'
            ORDER BY CAST(SUBSTR(T.alt_ac, position('.' in T.alt_ac) + 1,
                LENGTH(T.alt_ac)) AS INT) DESC;
            """  # noqa: S608
        result = await self.execute_query(query)
        if not result:
            msg = (
                f"Unable to find a result where {tx_ac} has transcript "
                f"coordinates {tx_exon_start} and {tx_exon_end} between "
                f"an exon's start and end coordinates"
            )
            if gene_query:
                msg += f" on gene {gene}"
            _logger.warning(msg)
            return None, msg
        return GenomicAlnData(**result[0]), None

    async def get_cds_start_end(self, tx_ac: str) -> tuple[int, int] | None:
        """Get coding start and end site

        :param tx_ac: Transcript accession
        :return: [Coding start site, Coding end site]
        """
        if tx_ac.startswith("ENS"):
            tx_ac = tx_ac.split(".")[0]
        query = f"""
            SELECT cds_start_i, cds_end_i
            FROM {self.schema}.transcript
            WHERE ac='{tx_ac}';
            """  # noqa: S608
        cds_start_end = await self.execute_query(query)
        if cds_start_end:
            cds_start_end = cds_start_end[0]
            if cds_start_end[0] is not None and cds_start_end[1] is not None:  # noqa: RET503
                return cds_start_end[0], cds_start_end[1]
        else:
            _logger.warning(
                "Unable to get coding start/end site for accession: %s", tx_ac
            )
            return None

    async def get_newest_assembly_ac(self, ac: str) -> list[str]:
        """Find accession associated to latest genomic assembly

        :param ac: Accession
        :return: List of accessions associated to latest genomic assembly. Order by
            desc
        """
        # Ensembl accessions do not have versions
        if ac.startswith("EN"):
            order_by_cond = "ORDER BY ac;"
        else:
            order_by_cond = (
                "ORDER BY SUBSTR(ac, 0, position('.' in ac)),"
                "CAST(SUBSTR(ac, position('.' in ac) + 1, LENGTH(ac)) AS INT) DESC;"
            )

        query = f"""
            SELECT ac
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac LIKE '{ac.split('.')[0]}%'
            AND ((descr IS NULL) OR (descr = ''))
            {order_by_cond}
            """  # noqa: S608
        results = await self.execute_query(query)
        if not results:
            return []

        return [r["ac"] for r in results]

    async def validate_genomic_ac(self, ac: str) -> bool:
        """Return whether or not genomic accession exists.

        :param ac: Genomic accession
        :return: ``True`` if genomic accession exists. ``False`` otherwise.
        """
        query = f"""
            SELECT EXISTS(
                SELECT ac
                FROM {self.schema}._seq_anno_most_recent
                WHERE ac = '{ac}'
            );
            """  # noqa: S608
        result = await self.execute_query(query)
        return result[0][0]

    async def gene_exists(self, gene: str) -> bool:
        """Return whether or not a gene symbol exists in UTA gene table

        :param gene: Gene symbol
        :return ``True`` if gene symbol exists in UTA, ``False`` if not
        """
        query = f"""
            SELECT EXISTS(
                SELECT hgnc
                FROM {self.schema}.gene
                WHERE hgnc = '{gene}'
            );
            """  # noqa: S608
        result = await self.execute_query(query)
        return result[0][0]

    async def transcript_exists(self, transcript: str) -> bool:
        """Return whether or not a transcript exists in the UTA tx_exon_aln_v table

        :param transcript: A transcript accession
        :return ``True`` if transcript exists in UTA, ``False`` if not
        """
        query = f"""
            SELECT EXISTS(
                SELECT tx_ac
                FROM {self.schema}.tx_exon_aln_v
                WHERE tx_ac = '{transcript}'
            );
            """  # noqa: S608
        result = await self.execute_query(query)
        return result[0][0]

    async def get_ac_descr(self, ac: str) -> str | None:
        """Return accession description. This is typically available only for accessions
        from older (pre-GRCh38) builds.

        >>> import asyncio
        >>> from cool_seq_tool.sources.uta_database import UtaDatabase
        >>> async def describe():
        ...     uta_db = await UtaDatabase.create()
        ...     result = await uta_db.get_ac_descr("NC_000001.10")
        ...     return result
        >>> asyncio.run(describe())
        'Homo sapiens chromosome 1, GRCh37.p13 Primary Assembly'

        :param ac: chromosome accession, e.g. ``"NC_000001.10"``
        :return: Description containing assembly and chromosome
        """
        query = f"""
            SELECT descr
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac = '{ac}';
            """  # noqa: S608
        result = await self.execute_query(query)
        if not result:
            _logger.warning("Accession %s does not have a description", ac)
            return None
        result = result[0][0]
        if result == "":
            result = None
        return result

    async def get_tx_exon_aln_v_data(
        self,
        tx_ac: str,
        start_pos: int,
        end_pos: int,
        alt_ac: str | None = None,
        use_tx_pos: bool = True,
        like_tx_ac: bool = False,
    ) -> list[TxExonAlnData]:
        """Return queried data from tx_exon_aln_v table.

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
        if tx_ac.startswith("EN"):
            temp_ac = tx_ac.split(".")[0]
            aln_method = f"AND alt_aln_method='genebuild'"  # noqa: F541
        else:
            temp_ac = tx_ac
            aln_method = f"AND alt_aln_method='splign'"  # noqa: F541

        if like_tx_ac:
            tx_q = f"WHERE tx_ac LIKE '{temp_ac}%'"
        else:
            tx_q = f"WHERE tx_ac='{temp_ac}'"

        order_by_cond = "ORDER BY CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1, LENGTH(alt_ac)) AS INT)"
        if alt_ac:
            alt_ac_q = f"AND alt_ac = '{alt_ac}'"
            if alt_ac.startswith("EN"):
                order_by_cond = "ORDER BY alt_ac"
        else:
            alt_ac_q = f"AND alt_ac LIKE 'NC_00%'"  # noqa: F541

        if use_tx_pos:
            pos_q = f"""tx_start_i AND tx_end_i"""  # noqa: F541
        else:
            pos_q = f"""alt_start_i AND alt_end_i"""  # noqa: F541

        query = f"""
            SELECT hgnc, tx_ac, tx_start_i, tx_end_i, alt_ac, alt_start_i,
                alt_end_i, alt_strand, alt_aln_method, ord, tx_exon_id, alt_exon_id
            FROM {self.schema}.tx_exon_aln_v
            {tx_q}
            {alt_ac_q}
            {aln_method}
            AND {start_pos} BETWEEN {pos_q}
            AND {end_pos} BETWEEN {pos_q}
            {order_by_cond}
            """  # noqa: S608
        result = await self.execute_query(query)
        if not result:
            _logger.warning("Unable to find transcript alignment for query: %s", query)
            return []
        if alt_ac and not use_tx_pos and len(result) > 1:
            _logger.debug(
                "Found more than one match for tx_ac %s and alt_ac = %s",
                temp_ac,
                alt_ac,
            )
        return [TxExonAlnData(**r) for r in result]

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
        representation.

        >>> import asyncio
        >>> from cool_seq_tool.sources import UtaDatabase
        >>> async def get_braf_mane():
        ...     uta_db = await UtaDatabase.create()
        ...     result = await uta_db.get_mane_c_genomic_data(
        ...         "NM_004333.6",
        ...         None,
        ...         140753335,
        ...         140753335,
        ...     )
        ...     return result
        >>> braf = asyncio.run(get_braf_mane())
        >>> braf["alt_ac"]
        'NC_000007.14'

        :param ac: MANE transcript accession
        :param alt_ac: NC accession. Used to triangulate on correct genomic data. Can
            be set to ``None`` if unavailable.
        :param start_pos: Genomic start position
        :param end_pos: Genomic end position change
        :return: Metadata for MANE genomic and transcript accessions results if
            successful
        """
        results = await self.get_tx_exon_aln_v_data(
            ac, start_pos, end_pos, alt_ac=alt_ac, use_tx_pos=False
        )
        if not results:
            return None
        result = results[0]

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
        :param pos: (start pos, end pos)
        :param annotation_layer: Annotation layer for ``ac`` and ``pos``
        :param alt_ac: Accession on g. coordinate
        :param target_genome_assembly: Genome assembly to get genomic data for.
            If ``alt_ac`` is provided, it will return the associated assembly.
        :return: Metadata for genomic and transcript accessions
        """
        results = await self.get_tx_exon_aln_v_data(
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
        else:
            if genomic_tx_data.strand == Strand.NEGATIVE:
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
        query = f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.genomic
            WHERE hgnc = '{gene}'
            AND alt_ac LIKE 'NC_00%'
            ORDER BY alt_ac;
            """  # noqa: S608

        records = await self.execute_query(query)
        if not records:
            return []

        alt_acs = [r["alt_ac"] for r in records]
        alt_acs.sort(key=lambda x: int(x.split(".")[-1]), reverse=True)
        return alt_acs

    async def get_gene_from_ac(
        self, ac: str, start_pos: int, end_pos: int
    ) -> list[str] | None:
        """Get gene(s) within the provided coordinate range

        >>> import asyncio
        >>> from cool_seq_tool.sources import UtaDatabase
        >>> async def get_gene():
        ...     uta_db = await UtaDatabase.create()
        ...     result = await uta_db.get_gene_from_ac(
        ...         "NC_000017.11", 43044296, 43045802
        ...     )
        ...     return result
        >>> asyncio.run(get_gene())
        ['BRCA1']

        :param ac: NC accession, e.g. ``"NC_000001.11"``
        :param start_pos: Start position change
        :param end_pos: End position change
        :return: List of HGNC gene symbols
        """
        if end_pos is None:
            end_pos = start_pos
        query = f"""
            SELECT DISTINCT hgnc
            FROM {self.schema}.genomic
            WHERE alt_ac = '{ac}'
            AND {start_pos} BETWEEN alt_start_i AND alt_end_i
            AND {end_pos} BETWEEN alt_start_i AND alt_end_i;
            """  # noqa: S608
        results = await self.execute_query(query)
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
                pos_cond = f"""
                    AND {start_pos} + T.cds_start_i
                        BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                    AND {end_pos} + T.cds_start_i
                        BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                    """
            else:
                pos_cond = f"""
                    AND {start_pos} BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                    AND {end_pos} BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                    """

        order_by_cond = """
        ORDER BY SUBSTR(ALIGN.alt_ac, 0, position('.' in ALIGN.alt_ac)),
        CAST(SUBSTR(ALIGN.alt_ac, position('.' in ALIGN.alt_ac) + 1,
            LENGTH(ALIGN.alt_ac)) AS INT) DESC,
        ALIGN.tx_end_i - ALIGN.tx_start_i DESC;
        """
        if alt_ac:
            alt_ac_cond = f"AND ALIGN.alt_ac = '{alt_ac}'"
            if alt_ac.startswith("EN"):
                order_by_cond = "ORDER BY ALIGN.alt_ac;"
        else:
            alt_ac_cond = "AND ALIGN.alt_ac LIKE 'NC_00%'"

        gene_cond = f"AND T.hgnc = '{gene}'" if gene else ""

        query = f"""
            SELECT AA.pro_ac, AA.tx_ac, ALIGN.alt_ac, T.cds_start_i
            FROM {self.schema}.associated_accessions as AA
            JOIN {self.schema}.transcript as T ON T.ac = AA.tx_ac
            JOIN {self.schema}.tx_exon_aln_v as ALIGN ON T.ac = ALIGN.tx_ac
            WHERE ALIGN.alt_aln_method = 'splign'
            {gene_cond}
            {alt_ac_cond}
            {pos_cond}
            {order_by_cond}
            """  # noqa: S608
        results = await self.execute_query(query)
        results = [
            (r["pro_ac"], r["tx_ac"], r["alt_ac"], r["cds_start_i"]) for r in results
        ]
        results_df = pl.DataFrame(results, schema=schema, orient="row")
        if results:
            results_df = results_df.unique()
        return results_df

    async def get_chr_assembly(self, ac: str) -> tuple[str, Assembly] | None:
        """Get chromosome and assembly for NC accession if not in GRCh38.

        >>> import asyncio
        >>> from cool_seq_tool.sources.uta_database import UtaDatabase
        >>> uta_db = asyncio.run(UtaDatabase.create())
        >>> result = asyncio.run(uta_db.get_chr_assembly("NC_000007.13"))
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
        except ValueError as e:
            _logger.error(e)
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
            FROM {self.schema}.associated_accessions
            WHERE pro_ac = '{p_ac}'
            {order_by_cond}
            """  # noqa: S608
        result = await self.execute_query(query)
        if result:
            result = [r["tx_ac"] for r in result]
        return result

    async def get_transcripts_from_genomic_pos(
        self, alt_ac: str, g_pos: int
    ) -> list[str]:
        """Get transcripts associated to a genomic ac and position.

        :param alt_ac: Genomic accession
        :param g_pos: Genomic position
        :return: RefSeq transcripts on c. coordinate
        """
        query = f"""
               SELECT distinct tx_ac
               FROM {self.schema}.tx_exon_aln_v
               WHERE alt_ac = '{alt_ac}'
               AND {g_pos} BETWEEN alt_start_i AND alt_end_i
               AND tx_ac LIKE 'NM_%';
               """  # noqa: S608
        results = await self.execute_query(query)
        if not results:
            return []
        return [item for sublist in results for item in sublist]

    @staticmethod
    def get_secret() -> str:
        """Get secrets for UTA DB instances. Used for deployment on AWS.

        :raises ClientError: If unable to retrieve secret value due to decryption
            decryption failure, internal service error, invalid parameter, invalid
            request, or resource not found.
        """
        secret_name = environ["UTA_DB_SECRET"]
        region_name = "us-east-2"

        # Create a Secrets Manager client
        session = boto3.session.Session()
        client = session.client(service_name="secretsmanager", region_name=region_name)

        try:
            get_secret_value_response = client.get_secret_value(SecretId=secret_name)
        except ClientError as e:
            # For a list of exceptions thrown, see
            # https://docs.aws.amazon.com/secretsmanager/latest/apireference/API_GetSecretValue.html
            _logger.error(e)
            raise e
        else:
            return get_secret_value_response["SecretString"]


class ParseResult(UrlLibParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.
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
        return path_elems[2] if len(path_elems) > 2 else None
