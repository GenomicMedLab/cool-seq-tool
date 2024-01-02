"""Provide transcript lookup and metadata tools via the UTA database."""
import ast
import base64
import logging
from os import environ
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar, Union
from urllib.parse import ParseResult as UrlLibParseResult
from urllib.parse import quote, unquote, urlparse

import asyncpg
import boto3
import polars as pl
from asyncpg.exceptions import InterfaceError, InvalidAuthorizationSpecificationError
from botocore.exceptions import ClientError
from pyliftover import LiftOver

from cool_seq_tool.schemas import AnnotationLayer, Assembly, Strand

# use `bound` to upper-bound UtaDatabase or child classes
UTADatabaseType = TypeVar("UTADatabaseType", bound="UtaDatabase")

# Environment variables for paths to chain files for pyliftover
LIFTOVER_CHAIN_37_TO_38 = environ.get("LIFTOVER_CHAIN_37_TO_38")
LIFTOVER_CHAIN_38_TO_37 = environ.get("LIFTOVER_CHAIN_38_TO_37")

UTA_DB_URL = environ.get(
    "UTA_DB_URL", "postgresql://uta_admin:uta@localhost:5433/uta/uta_20210129b"
)

logger = logging.getLogger(__name__)


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

    def __init__(
        self,
        db_url: str = UTA_DB_URL,
        chain_file_37_to_38: Optional[str] = None,
        chain_file_38_to_37: Optional[str] = None,
    ) -> None:
        """Initialize DB class. Should only be used by ``create()`` method, and not
        be called directly by a user.

        :param db_url: PostgreSQL connection URL
            Format: ``driver://user:password@host/database/schema``
        :param chain_file_37_to_38: Optional path to chain file for 37 to 38 assembly.
            This is used for ``pyliftover``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_37_TO_38`` env var is set. If neither is provided, will
            allow ``pyliftover`` to download a chain file from UCSC
        :param chain_file_38_to_37: Optional path to chain file for 38 to 37 assembly.
            This is used for ``pyliftover``. If this is not provided, will check to see
            if ``LIFTOVER_CHAIN_38_TO_37`` env var is set. If neither is provided, will
            allow ``pyliftover`` to download a chain file from UCSC
        """
        self.schema = None
        self._connection_pool = None
        original_pwd = db_url.split("//")[-1].split("@")[0].split(":")[-1]
        self.db_url = db_url.replace(original_pwd, quote(original_pwd))
        self.args = self._get_conn_args()

        chain_file_37_to_38 = chain_file_37_to_38 or LIFTOVER_CHAIN_37_TO_38
        if chain_file_37_to_38:
            self.liftover_37_to_38 = LiftOver(chain_file_37_to_38)
        else:
            self.liftover_37_to_38 = LiftOver("hg19", "hg38")

        chain_file_38_to_37 = chain_file_38_to_37 or LIFTOVER_CHAIN_38_TO_37
        if chain_file_38_to_37:
            self.liftover_38_to_37 = LiftOver(chain_file_38_to_37)
        else:
            self.liftover_38_to_37 = LiftOver("hg38", "hg19")

    def _get_conn_args(self) -> Dict:
        """Return connection arguments.

        :param db_url: raw connection URL
        :return: Database credentials
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
            environ[
                "UTA_DB_URL"
            ] = f"postgresql://{username}@{host}:{port}/{database}/{schema}"
            return dict(
                host=host,
                port=int(port),
                database=database,
                user=username,
                password=password,
            )
        else:
            url = ParseResult(urlparse(self.db_url))
            self.schema = url.schema
            password = unquote(url.password) if url.password else ""
            return dict(
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
                    host=self.args["host"],
                    port=self.args["port"],
                    user=self.args["user"],
                    password=self.args["password"],
                    database=self.args["database"],
                )
            except InterfaceError as e:
                logger.error(
                    f"While creating connection pool, " f"encountered exception {e}"
                )
                raise Exception("Could not create connection pool")

    @classmethod
    async def create(
        cls: Type[UTADatabaseType], db_url: str = UTA_DB_URL
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
            async with self._connection_pool.acquire() as connection:
                async with connection.transaction():
                    r = await connection.fetch(q)
                    return r

        if not self._connection_pool:
            await self.create_pool()
        try:
            result = await _execute_query(query)
            return result
        except InvalidAuthorizationSpecificationError:
            self._connection_pool = None
            await self.create_pool()
            result = await _execute_query(query)
            return result

    async def _create_genomic_table(self) -> None:
        """Create table containing genomic accession information."""
        check_table_exists = f"""
            SELECT EXISTS (
               SELECT FROM information_schema.tables
               WHERE table_schema = '{self.schema}'
               AND table_name = 'genomic'
            );
            """
        genomic_table_exists = await self.execute_query(check_table_exists)
        genomic_table_exists = genomic_table_exists[0].get("exists")
        if genomic_table_exists is None:
            logger.critical(
                "SELECT EXISTS query in UtaDatabase._create_genomic_table "
                "returned invalid response"
            )
            raise ValueError("SELECT EXISTS query returned invalid response")
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
                """
            await self.execute_query(create_genomic_table)

            indexes = [
                f"""CREATE INDEX alt_pos_index ON {self.schema}.genomic (alt_ac, alt_start_i, alt_end_i);""",
                f"""CREATE INDEX gene_alt_index ON {self.schema}.genomic (hgnc, alt_ac);""",
                f"""CREATE INDEX alt_ac_index ON {self.schema}.genomic (alt_ac);""",
            ]
            for create_index in indexes:
                await self.execute_query(create_index)

    @staticmethod
    def _transform_list(li: List) -> List[List[Any]]:
        """Transform list to only contain field values

        :param li: List of asyncpg.Record objects
        :return: List of list of objects
        """
        results = list()
        for item in li:
            results.append([field for field in item])
        return results

    async def get_genes_and_alt_acs(
        self,
        pos: int,
        strand: Optional[Strand] = None,
        chromosome: Optional[int] = None,
        alt_ac: Optional[str] = None,
        gene: Optional[str] = None,
    ) -> Tuple[Optional[Dict], Optional[str]]:
        """Return genes and genomic accessions for a position on a chromosome or alt_ac

        :param pos: Genomic position
        :param strand: Strand
        :param chromosome: Chromosome. Must give chromosome without a prefix
            (i.e. ``1`` or ``X``). If not provided, must provide ``alt_ac``.
            If ``alt_ac`` is also provided, ``alt_ac`` will be used.
        :param alt_ac: Genomic accession (i.e. ``NC_000001.11``). If not provided,
            must provide ``chromosome. If ``chromosome`` is also provided, ``alt_ac``
            will be used.
        :param gene: Gene symbol
        :return: Dictionary containing genes and genomic accessions and warnings if found
        """
        alt_ac_cond = (
            f"WHERE alt_ac = '{alt_ac}'"
            if alt_ac
            else f"WHERE alt_ac ~ '^NC_[0-9]+0{chromosome}.[0-9]+$'"
        )
        strand_cond = f"AND alt_strand = '{strand.value}'" if strand else ""
        gene_cond = f"AND hgnc = '{gene}'" if gene else ""

        query = f"""
            SELECT hgnc, alt_ac
            FROM {self.schema}.tx_exon_aln_v
            {alt_ac_cond}
            AND alt_aln_method = 'splign'
            AND {pos} BETWEEN alt_start_i AND alt_end_i
            {strand_cond}
            {gene_cond};
            """

        results = await self.execute_query(query)
        if not results:
            msg = (
                f"Unable to find a result for chromosome "
                f"{alt_ac or chromosome} where genomic coordinate {pos}"
                f" is mapped between an exon's start and end coordinates"
            )
            if strand:
                msg += (
                    f" on the "
                    f"{'positive' if strand == Strand.POSITIVE else 'negative'} strand"
                )
            if gene:
                msg += f" and on gene {gene}"
            return None, msg

        results = self._transform_list(results)
        genes = set()
        alt_acs = set()
        for r in results:
            genes.add(r[0])
            alt_acs.add(r[1])
        return dict(genes=genes, alt_acs=alt_acs), None

    async def get_tx_exons(
        self, tx_ac: str, alt_ac: Optional[str] = None
    ) -> Tuple[Optional[List[Tuple[int, int]]], Optional[str]]:
        """Get list of transcript exons start/end coordinates.

        :param tx_ac: Transcript accession
        :param alt_ac: Genomic accession
        :return: List of a transcript's accessions and warnings if found
        """
        if alt_ac:
            # We know what assembly we're looking for since we have the
            # genomic accession
            query = f"""
                SELECT DISTINCT tx_start_i, tx_end_i
                FROM {self.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_aln_method = 'splign'
                AND alt_ac = '{alt_ac}'
                """
        else:
            # Use GRCh38 by default if no genomic accession is provided
            query = f"""
                SELECT DISTINCT tx_start_i, tx_end_i
                FROM {self.schema}.tx_exon_aln_v as t
                INNER JOIN {self.schema}._seq_anno_most_recent as s
                ON t.alt_ac = s.ac
                WHERE s.descr = ''
                AND t.tx_ac = '{tx_ac}'
                AND t.alt_aln_method = 'splign'
                AND t.alt_ac like 'NC_000%'
                """
        result = await self.execute_query(query)

        if not result:
            msg = f"Unable to get exons for {tx_ac}"
            logger.warning(msg)
            return None, msg
        else:
            tx_exons = [(r["tx_start_i"], r["tx_end_i"]) for r in result]
            return tx_exons, None

    async def get_alt_ac_start_or_end(
        self, tx_ac: str, tx_exon_start: int, tx_exon_end: int, gene: Optional[str]
    ) -> Tuple[Optional[Tuple[str, str, int, int, int]], Optional[str]]:
        """Get genomic data for related transcript exon start or end.

        :param tx_ac: Transcript accession
        :param tx_exon_start: Transcript's exon start coordinate
        :param tx_exon_end: Transcript's exon end coordinate
        :param gene: HGNC gene symbol
        :return: [hgnc symbol, genomic accession for chromosome,
            aligned genomic start coordinate, aligned genomic end coordinate, strand],
            and warnings if found
        """
        if gene:
            gene_query = f"AND T.hgnc = '{gene}'"
        else:
            gene_query = ""

        query = f"""
            SELECT T.hgnc, T.alt_ac, T.alt_start_i, T.alt_end_i, T.alt_strand
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
            """
        result = await self.execute_query(query)
        if not result:
            msg = (
                f"Unable to find a result where {tx_ac} has transcript "
                f"coordinates {tx_exon_start} and {tx_exon_end} between "
                f"an exon's start and end coordinates"
            )
            if gene_query:
                msg += f" on gene {gene}"
            logger.warning(msg)
            return None, msg
        else:
            result = result[0]
        return (result[0], result[1], result[2], result[3], result[4]), None

    async def get_cds_start_end(self, tx_ac: str) -> Optional[Tuple[int, int]]:
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
            """
        cds_start_end = await self.execute_query(query)
        if cds_start_end:
            cds_start_end = cds_start_end[0]
            if cds_start_end[0] is not None and cds_start_end[1] is not None:
                return cds_start_end[0], cds_start_end[1]
        else:
            logger.warning(
                f"Unable to get coding start/end site for " f"accession: {tx_ac}"
            )
            return None

    async def get_newest_assembly_ac(self, ac: str) -> List[str]:
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
            """
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
            """
        result = await self.execute_query(query)
        return result[0][0]

    async def get_ac_descr(self, ac: str) -> Optional[str]:
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
            """
        result = await self.execute_query(query)
        if not result:
            logger.warning(f"Accession {ac} does not have a description")
            return None
        else:
            result = result[0][0]
            if result == "":
                result = None
            return result

    async def get_tx_exon_aln_v_data(
        self,
        tx_ac: str,
        start_pos: int,
        end_pos: int,
        alt_ac: Optional[str] = None,
        use_tx_pos: bool = True,
        like_tx_ac: bool = False,
    ) -> List:
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
        :return: List of tx_exon_aln_v data
        """
        if end_pos is None:
            end_pos = start_pos

        if tx_ac.startswith("EN"):
            temp_ac = tx_ac.split(".")[0]
            aln_method = f"AND alt_aln_method='genebuild'"  # noqa: F541
        else:
            temp_ac = tx_ac
            aln_method = f"AND alt_aln_method='splign'"  # noqa: F541

        if like_tx_ac:
            tx_q = f"WHERE tx_ac LIKE '{temp_ac}%'"  # noqa: F541
        else:
            tx_q = f"WHERE tx_ac='{temp_ac}'"  # noqa: F541

        order_by_cond = "ORDER BY CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1, LENGTH(alt_ac)) AS INT)"  # noqa: E501
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
                alt_end_i, alt_strand, alt_aln_method, tx_exon_id, alt_exon_id
            FROM {self.schema}.tx_exon_aln_v
            {tx_q}
            {alt_ac_q}
            {aln_method}
            AND {start_pos} BETWEEN {pos_q}
            AND {end_pos} BETWEEN {pos_q}
            {order_by_cond}
            """
        result = await self.execute_query(query)
        if not result:
            logger.warning(
                f"Unable to find transcript alignment for query: " f"{query}"
            )
            return []
        if alt_ac and not use_tx_pos:
            if len(result) > 1:
                logger.debug(
                    f"Found more than one match for tx_ac {temp_ac} "
                    f"and alt_ac = {alt_ac}"
                )
        results = list()
        for r in result:
            results.append([field for field in r])
        return results

    @staticmethod
    def data_from_result(result: List) -> Optional[Dict]:
        """Return data found from result.

        :param result: Data from tx_exon_aln_v table
        :return: Gene, strand, and position ranges for tx and alt_ac
        """
        gene = result[0]
        tx_pos_range = result[2], result[3]
        alt_pos_range = result[5], result[6]
        strand = Strand(result[7])
        alt_aln_method = result[8]
        tx_exon_id = result[9]
        alt_exon_id = result[10]

        if (tx_pos_range[1] - tx_pos_range[0]) != (alt_pos_range[1] - alt_pos_range[0]):
            logger.warning(
                f"tx_pos_range {tx_pos_range} "
                f"is not the same length as alt_pos_range "
                f"{alt_pos_range}."
            )
            return None

        return dict(
            gene=gene,
            strand=strand,
            tx_pos_range=tx_pos_range,
            alt_pos_range=alt_pos_range,
            alt_aln_method=alt_aln_method,
            tx_exon_id=tx_exon_id,
            alt_exon_id=alt_exon_id,
        )

    async def get_mane_c_genomic_data(
        self, ac: str, alt_ac: Optional[str], start_pos: int, end_pos: int
    ) -> Optional[Dict]:
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
        :return: MANE transcript results if successful
        """
        results = await self.get_tx_exon_aln_v_data(
            ac, start_pos, end_pos, alt_ac=alt_ac, use_tx_pos=False
        )
        if not results:
            return None
        result = results[0]

        data = self.data_from_result(result)
        if not data:
            return None

        coding_start_site = await self.get_cds_start_end(ac)
        if coding_start_site is None:
            logger.warning(f"Accession {ac} not found in UTA")
            return None

        data["tx_ac"] = result[1]
        data["alt_ac"] = result[4]
        data["coding_start_site"] = coding_start_site[0]
        data["coding_end_site"] = coding_start_site[1]

        if data["strand"] == Strand.NEGATIVE:
            data["alt_pos_change_range"] = (end_pos, start_pos)
            data["alt_pos_change"] = (
                data["alt_pos_range"][1] - data["alt_pos_change_range"][0],
                data["alt_pos_change_range"][1] - data["alt_pos_range"][0],
            )
        else:
            data["alt_pos_change_range"] = (start_pos, end_pos)
            data["alt_pos_change"] = (
                data["alt_pos_change_range"][0] - data["alt_pos_range"][0],
                data["alt_pos_range"][1] - data["alt_pos_change_range"][1],
            )

        return data

    async def get_genomic_tx_data(
        self,
        tx_ac: str,
        pos: Tuple[int, int],
        annotation_layer: Union[
            AnnotationLayer.CDNA, AnnotationLayer.GENOMIC
        ] = AnnotationLayer.CDNA,
        alt_ac: Optional[str] = None,
        target_genome_assembly: Assembly = Assembly.GRCH38,
    ) -> Optional[Dict]:
        """Get transcript mapping to genomic data.

        :param tx_ac: Accession on c. coordinate
        :param pos: (start pos, end pos)
        :param annotation_layer: Annotation layer for ``ac`` and ``pos``
        :param alt_ac: Accession on g. coordinate
        :param target_genome_assembly: Genome assembly to get genomic data for.
            If ``alt_ac`` is provided, it will return the associated assembly.
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
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

        data = self.data_from_result(result)
        if not data:
            return None
        data["tx_ac"] = result[1]
        data["alt_ac"] = result[4]

        data["pos_change"] = (
            pos[0] - data["tx_pos_range"][0],
            data["tx_pos_range"][1] - pos[1],
        )

        if annotation_layer == AnnotationLayer.CDNA:
            if data["strand"] == Strand.NEGATIVE:
                data["alt_pos_change_range"] = (
                    data["alt_pos_range"][1] - data["pos_change"][0],
                    data["alt_pos_range"][0] + data["pos_change"][1],
                )
            else:
                data["alt_pos_change_range"] = (
                    data["alt_pos_range"][0] + data["pos_change"][0],
                    data["alt_pos_range"][1] - data["pos_change"][1],
                )
        else:
            if data["strand"] == Strand.NEGATIVE:
                data["alt_pos_change_range"] = (pos[1], pos[0])
            else:
                data["alt_pos_change_range"] = pos

        return data

    async def get_ac_from_gene(self, gene: str) -> List[str]:
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
            """

        records = await self.execute_query(query)
        if not records:
            return []

        alt_acs = [r["alt_ac"] for r in records]
        alt_acs.sort(key=lambda x: int(x.split(".")[-1]), reverse=True)
        return alt_acs

    async def get_gene_from_ac(
        self, ac: str, start_pos: int, end_pos: int
    ) -> Optional[List[str]]:
        """Get gene(s) within the provided coordinate range

        >>> import asyncio
        >>> from cool_seq_tool.sources import UtaDatabase
        >>> async def get_gene():
        ...     uta_db = await UtaDatabase.create()
        ...     result = await uta_db.get_gene_from_ac("NC_000017.11", 43044296, 43045802)
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
            """
        results = await self.execute_query(query)
        if not results:
            logger.warning(
                f"Unable to find gene between {start_pos} and" f" {end_pos} on {ac}"
            )
            return None
        else:
            if len(results) > 1:
                logger.info(
                    f"Found more than one gene between "
                    f"{start_pos} and {end_pos} on {ac}"
                )

        return [r[0] for r in results]

    async def get_transcripts(
        self,
        start_pos: Optional[int] = None,
        end_pos: Optional[int] = None,
        gene: Optional[str] = None,
        use_tx_pos: bool = True,
        alt_ac: Optional[str] = None,
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
            """
        results = await self.execute_query(query)
        results = [
            (r["pro_ac"], r["tx_ac"], r["alt_ac"], r["cds_start_i"]) for r in results
        ]
        results_df = pl.DataFrame(results, schema=schema)
        if results:
            results_df = results_df.unique()
        return results_df

    async def get_chr_assembly(self, ac: str) -> Optional[Tuple[str, str]]:
        """Get chromosome and assembly for NC accession if not in GRCh38.

        :param ac: NC accession
        :return: Chromosome and Assembly accession is on
        """
        descr = await self.get_ac_descr(ac)
        if not descr:
            # Already GRCh38 Assembly
            return None
        descr = descr.split(",")
        chromosome = f"chr{descr[0].split()[-1]}"
        assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

        if assembly not in ["GRCh37", "GRCh38"]:
            logger.warning(
                f"Assembly not supported: {assembly}. "
                f"Only GRCh37 and GRCh38 are supported."
            )
            return None

        return chromosome, assembly

    async def liftover_to_38(self, genomic_tx_data: Dict) -> None:
        """Liftover genomic_tx_data to hg38 assembly.

        :param genomic_tx_data: Dictionary containing gene, nc_accession, alt_pos, and
            strand
        """
        descr = await self.get_chr_assembly(genomic_tx_data["alt_ac"])
        if descr is None:
            # already grch38
            return None
        chromosome, _ = descr

        query = f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac = '{genomic_tx_data['tx_ac']}';
            """
        nc_acs = await self.execute_query(query)
        nc_acs = [nc_ac[0] for nc_ac in nc_acs]
        if nc_acs == [genomic_tx_data["alt_ac"]]:
            logger.warning(
                f"UTA does not have GRCh38 assembly for "
                f"{genomic_tx_data['alt_ac'].split('.')[0]}"
            )
            return None

        # Get most recent assembly version position
        # Liftover range
        self._set_liftover(
            genomic_tx_data, "alt_pos_range", chromosome, Assembly.GRCH38
        )

        # Liftover changes range
        self._set_liftover(
            genomic_tx_data, "alt_pos_change_range", chromosome, Assembly.GRCH38
        )

        # Change alt_ac to most recent
        if genomic_tx_data["alt_ac"].startswith("EN"):
            order_by_cond = "ORDER BY alt_ac DESC;"
        else:
            order_by_cond = """
            ORDER BY CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1,
            LENGTH(alt_ac)) AS INT) DESC;
            """
        query = f"""
            SELECT alt_ac
            FROM {self.schema}.genomic
            WHERE alt_ac LIKE '{genomic_tx_data['alt_ac'].split('.')[0]}%'
            {order_by_cond}
            """
        nc_acs = await self.execute_query(query)
        genomic_tx_data["alt_ac"] = nc_acs[0][0]

    def get_liftover(
        self, chromosome: str, pos: int, liftover_to_assembly: Assembly
    ) -> Optional[Tuple]:
        """Get new genome assembly data for a position on a chromosome.

        :param chromosome: The chromosome number. Must be prefixed with ``chr``
        :param pos: Position on the chromosome
        :param liftover_to_assembly: Assembly to liftover to
        :return: [Target chromosome, target position, target strand,
            conversion_chain_score] for assembly
        """
        if not chromosome.startswith("chr"):
            logger.warning("`chromosome` must be prefixed with chr")
            return None

        if liftover_to_assembly == Assembly.GRCH38:
            liftover = self.liftover_37_to_38.convert_coordinate(chromosome, pos)
        elif liftover_to_assembly == Assembly.GRCH37:
            liftover = self.liftover_38_to_37.convert_coordinate(chromosome, pos)
        else:
            logger.warning(f"{liftover_to_assembly} assembly not supported")
            liftover = None

        if liftover is None or len(liftover) == 0:
            logger.warning(f"{pos} does not exist on {chromosome}")
            return None
        else:
            return liftover[0]

    def _set_liftover(
        self,
        genomic_tx_data: Dict,
        key: str,
        chromosome: str,
        liftover_to_assembly: Assembly,
    ) -> None:
        """Update genomic_tx_data to have coordinates for given assembly.

        :param genomic_tx_data: Dictionary containing gene, nc_accession, alt_pos, and
            strand
        :param key: Key to access coordinate positions
        :param chromosome: Chromosome, must be prefixed with ``chr``
        :param liftover_to_assembly: Assembly to liftover to
        """
        liftover_start_i = self.get_liftover(
            chromosome, genomic_tx_data[key][0], liftover_to_assembly
        )
        if liftover_start_i is None:
            logger.warning(
                f"Unable to liftover position "
                f"{genomic_tx_data[key][0]} on {chromosome}"
            )
            return None

        liftover_end_i = self.get_liftover(
            chromosome, genomic_tx_data[key][1], liftover_to_assembly
        )
        if liftover_end_i is None:
            logger.warning(
                f"Unable to liftover position "
                f"{genomic_tx_data[key][1]} on {chromosome}"
            )
            return None

        genomic_tx_data[key] = liftover_start_i[1], liftover_end_i[1]

    async def p_to_c_ac(self, p_ac: str) -> List[str]:
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
            """
        result = await self.execute_query(query)
        if result:
            result = [r["tx_ac"] for r in result]
        return result

    async def get_transcripts_from_genomic_pos(
        self, alt_ac: str, g_pos: int
    ) -> List[str]:
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
               """
        results = await self.execute_query(query)
        if not results:
            return []
        return [item for sublist in results for item in sublist]

    @staticmethod
    def get_secret() -> str:
        """Get secrets for UTA DB instances. Used for deployment on AWS."""
        secret_name = environ["UTA_DB_SECRET"]
        region_name = "us-east-2"

        # Create a Secrets Manager client
        session = boto3.session.Session()
        client = session.client(service_name="secretsmanager", region_name=region_name)

        try:
            get_secret_value_response = client.get_secret_value(SecretId=secret_name)
        except ClientError as e:
            logger.warning(e)
            if e.response["Error"]["Code"] == "DecryptionFailureException":
                # Secrets Manager can"t decrypt the protected
                # secret text using the provided KMS key.
                raise e
            elif e.response["Error"]["Code"] == "InternalServiceErrorException":
                # An error occurred on the server side.
                raise e
            elif e.response["Error"]["Code"] == "InvalidParameterException":
                # You provided an invalid value for a parameter.
                raise e
            elif e.response["Error"]["Code"] == "InvalidRequestException":
                # You provided a parameter value that is not valid for
                # the current state of the resource.
                raise e
            elif e.response["Error"]["Code"] == "ResourceNotFoundException":
                # We can"t find the resource that you asked for.
                raise e
        else:
            # Decrypts secret using the associated KMS CMK.
            # Depending on whether the secret is a string or binary,
            # one of these fields will be populated.
            if "SecretString" in get_secret_value_response:
                secret = get_secret_value_response["SecretString"]
                return secret
            else:
                decoded_binary_secret = base64.b64decode(
                    get_secret_value_response["SecretBinary"]
                )
                return decoded_binary_secret


class ParseResult(UrlLibParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.
    Source: https://github.com/biocommons/hgvs
    """

    def __new__(cls, pr):  # noqa
        """Create new instance."""
        return super(ParseResult, cls).__new__(cls, *pr)

    @property
    def database(self) -> Optional[str]:
        """Create database property."""
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self) -> Optional[str]:
        """Create schema property."""
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None
