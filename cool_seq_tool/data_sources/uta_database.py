"""Module for UTA queries."""
import ast
import base64
from typing import Dict, List, Optional, Tuple, Any, TypeVar, Type, Union
from os import environ
from urllib.parse import quote, unquote

from six.moves.urllib import parse as urlparse
import pandas as pd
import asyncpg
import boto3
from pyliftover import LiftOver
from asyncpg.exceptions import InvalidAuthorizationSpecificationError, \
    InterfaceError
from botocore.exceptions import ClientError

from cool_seq_tool import UTA_DB_URL, logger
from cool_seq_tool.schemas import AnnotationLayer, Assembly


# use `bound` to upper-bound UTADatabase or child classes
UTADatabaseType = TypeVar("UTADatabaseType", bound="UTADatabase")


class UTADatabase:
    """Class for connecting and querying UTA database."""

    def __init__(self, db_url: str = UTA_DB_URL, db_pwd: str = "") -> None:
        """Initialize DB class. Downstream libraries should use the create()
        method to construct a new instance: await UTADatabase.create()

        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        """
        self.schema = None
        self.db_url = db_url
        self.db_pwd = db_pwd
        self._connection_pool = None
        self.args = self._get_conn_args()
        self.liftover_37_to_38 = LiftOver("hg19", "hg38")
        self.liftover_38_to_37 = LiftOver("hg38", "hg19")

    @staticmethod
    def _update_db_url(db_pwd: str, db_url: str) -> str:
        """Return new db_url containing password.

        :param str db_pwd: User's password for uta database
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :return: PostgreSQL connection URL
        """
        if "UTA_DB_URL" in environ:
            return environ["UTA_DB_URL"]
        if not db_pwd and "UTA_PASSWORD" not in environ:
            raise Exception("Environment variable UTA_PASSWORD "
                            "or `db_pwd` param must be set")
        else:
            uta_password_in_environ = "UTA_PASSWORD" in environ
            if uta_password_in_environ and db_pwd:
                if db_pwd != environ["UTA_PASSWORD"]:
                    raise Exception("If both environment variable UTA_PASSWORD"
                                    " and param db_pwd is set, they must "
                                    "both be the same")
            else:
                if uta_password_in_environ and not db_pwd:
                    db_pwd = environ["UTA_PASSWORD"]
            db_url = db_url.split("@")
            db_url_with_pass = f"{db_url[0]}:{db_pwd}@{db_url[1]}"
            environ["UTA_DB_URL"] = db_url_with_pass
            return db_url_with_pass

    def _get_conn_args(self) -> Dict:
        """Return connection arguments.

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
            environ["UTA_DB_URL"] = f"postgresql://{username}@{host}:{port}/{database}/{schema}"  # noqa: E501
            return dict(host=host, port=int(port), database=database, user=username,
                        password=password)
        else:
            db_url = self._update_db_url(self.db_pwd, self.db_url)
            original_pwd = db_url.split("//")[-1].split("@")[0].split(":")[-1]
            db_url = db_url.replace(original_pwd, quote(original_pwd))

            url = ParseResult(urlparse.urlparse(db_url))
            self.schema = url.schema
            password = unquote(url.password) if url.password else ""
            return dict(host=url.hostname, port=url.port,
                        database=url.database, user=url.username,
                        password=password)

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
                logger.error(f"While creating connection pool, "
                             f"encountered exception {e}")
                raise Exception("Could not create connection pool")

    @classmethod
    async def create(
            cls: Type[UTADatabaseType], db_url: str = UTA_DB_URL,
            db_pwd: str = "") -> UTADatabaseType:
        """Provide fully-initialized class instance (a la factory pattern)
        :param UTADatabaseType cls: supplied implicitly
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :return: UTA DB access class instance
        """
        self = cls(db_url, db_pwd)
        await self._create_genomic_table()
        await self.create_pool()
        return self

    async def execute_query(self, query: str) -> Any:  # noqa: ANN401
        """Execute a query and return its result.

        :param str query: Query to make on database
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
        check_table_exists = (
            f"""
            SELECT EXISTS (
               SELECT FROM information_schema.tables
               WHERE table_schema = '{self.schema}'
               AND table_name = 'genomic'
            );
            """
        )
        genomic_table_exists = await self.execute_query(check_table_exists)
        genomic_table_exists = genomic_table_exists[0].get("exists")
        if genomic_table_exists is None:
            logger.critical(
                "SELECT EXISTS query in UTADatabase._create_genomic_table "
                "returned invalid response"
            )
            raise ValueError("SELECT EXISTS query returned invalid response")
        if not genomic_table_exists:
            create_genomic_table = (
                f"""
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
            )
            await self.execute_query(create_genomic_table)

            indexes = [
                f"""CREATE INDEX alt_pos_index ON {self.schema}.genomic (alt_ac, alt_start_i, alt_end_i);""",  # noqa: E501
                f"""CREATE INDEX gene_alt_index ON {self.schema}.genomic (hgnc, alt_ac);""",  # noqa: E501
                f"""CREATE INDEX alt_ac_index ON {self.schema}.genomic (alt_ac);"""  # noqa: E501
            ]
            for create_index in indexes:
                await self.execute_query(create_index)

    @staticmethod
    def _transform_list(li: List) -> List[List[Any]]:
        """Transform list to only contain field values

        :param List li: List of asyncpg.Record objects
        :return: List of list of objects
        """
        results = list()
        for item in li:
            results.append([field for field in item])
        return results

    async def chr_to_gene_and_accessions(
            self, chromosome: int, pos: int, strand: Optional[int] = None,
            alt_ac: Optional[str] = None,
            gene: Optional[str] = None) -> Tuple[Optional[Dict], Optional[str]]:
        """Return genes and genomic accessions related to a position on a chr.

        :param int chromosome: Chromosome number
        :param int pos: Genomic position
        :param Optional[int] strand: Strand. Must be either `-1` or `1`
        :param Optional[str] alt_ac: Genomic accession
        :param Optional[str] gene: Gene symbol
        :return: Dictionary containing genes and genomic accessions and
            warnings if found
        """
        alt_ac_cond = f"WHERE alt_ac = '{alt_ac}'" if alt_ac else f"WHERE alt_ac ~ '^NC_[0-9]+0{chromosome}.[0-9]+$'"  # noqa: E501
        strand_cond = f"AND alt_strand = '{strand}'" if strand else ""
        gene_cond = f"AND hgnc = '{gene}'" if gene else ""

        query = (
            f"""
            SELECT hgnc, alt_ac
            FROM {self.schema}.tx_exon_aln_v
            {alt_ac_cond}
            AND alt_aln_method = 'splign'
            AND {pos} BETWEEN alt_start_i AND alt_end_i
            {strand_cond}
            {gene_cond};
            """
        )

        results = await self.execute_query(query)
        if not results:
            msg = f"Unable to find a result for chromosome " \
                  f"{alt_ac or chromosome} where genomic coordinate {pos}" \
                  f" is mapped between an exon's start and end coordinates"
            if strand:
                msg += f" on the " \
                       f"{'positive' if strand == 1 else 'negative'} strand"
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
    ) -> Tuple[Optional[List[Tuple[int, int]]], Optional[str]]:  # noqa: E501
        """Get list of transcript exons start/end coordinates.

        :param str tx_ac: Transcript accession
        :param Optional[str] alt_ac: Genomic accession
        :return: List of a transcript's accessions and warnings if found
        """
        if alt_ac:
            # We know what asesmbly we're looking for since we have the
            # genomic accession
            query = (
                f"""
                SELECT DISTINCT tx_start_i, tx_end_i
                FROM {self.schema}.tx_exon_aln_v
                WHERE tx_ac = '{tx_ac}'
                AND alt_aln_method = 'splign'
                AND alt_ac = '{alt_ac}'
                """
            )
        else:
            # Use GRCh38 by default if no genomic accession is provided
            query = (
                f"""
                SELECT DISTINCT tx_start_i, tx_end_i
                FROM {self.schema}.tx_exon_aln_v as t
                INNER JOIN {self.schema}._seq_anno_most_recent as s
                ON t.alt_ac = s.ac
                WHERE s.descr = ''
                AND t.tx_ac = '{tx_ac}'
                AND t.alt_aln_method = 'splign'
                AND t.alt_ac like 'NC_000%'
                """
            )
        result = await self.execute_query(query)

        if not result:
            msg = f"Unable to get exons for {tx_ac}"
            logger.warning(msg)
            return None, msg
        else:
            tx_exons = [(r["tx_start_i"], r["tx_end_i"]) for r in result]
            return tx_exons, None

    @staticmethod
    def _validate_exon(
            transcript: str, tx_exons: List[str],
            exon_number: Optional[int] = None) -> Tuple[Optional[List], Optional[str]]:  # noqa: E501
        """Validate that exon number is valid

        :param str transcript: Transcript accession
        :param List tx_exons: List of transcript's exons
        :param Optional[int] exon_number: Exon number to validate
        :return: Transcript coordinates and warnings if found
        """
        msg = f"Exon {exon_number} does not exist on {transcript}"
        try:
            if exon_number < 1:
                return None, msg
            exon = tx_exons[exon_number - 1]
        except IndexError:
            return None, msg
        return exon, None

    def get_tx_exon_coords(
            self, transcript: str, tx_exons: List[str],
            exon_start: Optional[int] = None,
            exon_end: Optional[int] = None) -> Tuple[Optional[Tuple[List, List]], Optional[str]]:  # noqa: E501
        """Get transcript exon coordinates

        :param str transcript: Transcript accession
        :param List tx_exons: List of transcript exons
        :param Optional[int] exon_start: Start exon number
        :param Optional[int] exon_end: End exon number
        :return: [Transcript start exon coords, Transcript end exon coords],
            and warnings if found
        """
        if exon_start is not None:
            tx_exon_start, warning = self._validate_exon(
                transcript, tx_exons, exon_start)
            if not tx_exon_start:
                return None, warning
        else:
            tx_exon_start = None

        if exon_end is not None:
            tx_exon_end, warning = self._validate_exon(
                transcript, tx_exons, exon_end)
            if not tx_exon_end:
                return None, warning
        else:
            tx_exon_end = None
        return (tx_exon_start, tx_exon_end), None

    async def get_alt_ac_start_and_end(
            self, tx_ac: str, tx_exon_start: Optional[List[str]] = None,
            tx_exon_end: Optional[List[str]] = None,
            gene: str = None) -> Tuple[Optional[Tuple[Tuple, Tuple]], Optional[str]]:  # noqa: E501
        """Get genomic coordinates for related transcript exon start and end.

        :param str tx_ac: Transcript accession
        :param Optional[List[str]] tx_exon_start: Transcript's exon start
            coordinates
        :param Optional[List[str]] tx_exon_end: Transcript's exon end
            coordinates
        :param str gene: Gene symbol
        :return: Alt ac start and end data, and warnings if found
        """
        if tx_exon_start:
            alt_ac_start, warning = await self.get_alt_ac_start_or_end(
                tx_ac, int(tx_exon_start[0]), int(tx_exon_start[1]), gene=gene)
            if not alt_ac_start:
                return None, warning
        else:
            alt_ac_start = None

        if tx_exon_end:
            alt_ac_end, warning = await self.get_alt_ac_start_or_end(
                tx_ac, int(tx_exon_end[0]), int(tx_exon_end[1]), gene=gene)
            if not alt_ac_end:
                return None, warning
        else:
            alt_ac_end = None

        if alt_ac_start is None and alt_ac_end is None:
            msg = "Unable to find `alt_ac_start` or `alt_ac_end`"
            logger.warning(msg)
            return None, msg

        # validate
        if alt_ac_start and alt_ac_end:
            for i in (0, 1, 4):
                if alt_ac_start[i] != alt_ac_end[i]:
                    if i == 0:
                        error = "Gene symbol does not match"
                    elif i == 1:
                        error = "Chromosome does not match"
                    else:
                        error = "Strand does not match"
                    logger.warning(f"{error}: "
                                   f"{alt_ac_start[i]} != {alt_ac_end[i]}")
        return (alt_ac_start, alt_ac_end), None

    async def get_alt_ac_start_or_end(self, tx_ac: str, tx_exon_start: int,
                                      tx_exon_end: int, gene: Optional[str])\
            -> Tuple[Optional[Tuple[str, str, int, int, int]], Optional[str]]:
        """Get genomic data for related transcript exon start or end.

        :param str tx_ac: Transcript accession
        :param int tx_exon_start: Transcript's exon start coordinate
        :param int tx_exon_end: Transcript's exon end coordinate
        :param Optional[str] gene: Gene symbol
        :return: [hgnc symbol, genomic accession for chromosome,
            start exon's end coordinate, end exon's start coordinate, strand],
            and warnings if found
        """
        if gene:
            gene_query = f"AND T.hgnc = '{gene}'"
        else:
            gene_query = ""

        query = (
            f"""
            SELECT T.hgnc, T.alt_ac, T.alt_start_i, T.alt_end_i, T.alt_strand
            FROM {self.schema}._cds_exons_fp_v as C
            JOIN {self.schema}.tx_exon_aln_v as T ON T.tx_ac = C.tx_ac
            WHERE T.tx_ac = '{tx_ac}'
            {gene_query}
            AND {tx_exon_start} BETWEEN T.tx_start_i AND T.tx_end_i
            AND {tx_exon_end} BETWEEN T.tx_start_i AND T.tx_end_i
            AND T.alt_aln_method = 'splign'
            AND T.alt_ac LIKE 'NC_00%'
            ORDER BY (CAST(SUBSTR(T.alt_ac, position('.' in T.alt_ac) + 1,
                LENGTH(T.alt_ac)) AS INT)) DESC;
            """
        )
        result = await self.execute_query(query)
        if not result:
            msg = f"Unable to find a result where {tx_ac} has transcript " \
                  f"coordinates {tx_exon_start} and {tx_exon_end} between " \
                  f"an exon's start and end coordinates"
            if gene_query:
                msg += f" on gene {gene}"
            logger.warning(msg)
            return None, msg
        else:
            result = result[0]
        return (result[0], result[1], result[2], result[3], result[4]), None

    async def get_cds_start_end(self, tx_ac: str) -> Optional[Tuple[int, int]]:
        """Get coding start and end site

        :param str tx_ac: Transcript accession
        :return: [Coding start site, Coding end site]
        """
        if tx_ac.startswith("ENS"):
            tx_ac = tx_ac.split(".")[0]
        query = (
            f"""
            SELECT cds_start_i, cds_end_i
            FROM {self.schema}.transcript
            WHERE ac='{tx_ac}';
            """
        )
        cds_start_end = await self.execute_query(query)
        if cds_start_end:
            cds_start_end = cds_start_end[0]
            if cds_start_end[0] is not None \
                    and cds_start_end[1] is not None:
                return cds_start_end[0], cds_start_end[1]
        else:
            logger.warning(f"Unable to get coding start/end site for "
                           f"accession: {tx_ac}")
            return None

    async def get_newest_assembly_ac(self, ac: str) -> List:
        """Find most recent genomic accession version.

        :param str ac: Genomic accession
        :return: List of most recent genomic accession version
        """
        query = (
            f"""
            SELECT ac
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac LIKE '{ac.split('.')[0]}%'
            AND ((descr IS NULL) OR (descr = ''))
            ORDER BY ac;
            """
        )
        result = await self.execute_query(query)
        if result:
            ret = list()
            for r in result:
                ret.append([field for field in r])
            return ret
        return []

    async def validate_genomic_ac(self, ac: str) -> bool:
        """Return whether or not genomic accession exists.

        :param str ac: Genomic accession
        :return: `True` if genomic accession exists. `False` otherwise.
        """
        query = (
            f"""
            SELECT EXISTS(
                SELECT ac
                FROM {self.schema}._seq_anno_most_recent
                WHERE ac = '{ac}'
            );
            """
        )
        result = await self.execute_query(query)
        return result[0][0]

    async def get_ac_descr(self, ac: str) -> Optional[str]:
        """Return accession description.
        Typically description exists if not GRCh38 assembly.

        :param str ac: Accession
        :return: Description containing assembly and chromosome
        """
        query = (
            f"""
            SELECT descr
            FROM {self.schema}._seq_anno_most_recent
            WHERE ac = '{ac}';
            """
        )
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
        self, tx_ac: str, start_pos: int, end_pos: int, alt_ac: str = None,
        use_tx_pos: bool = True, like_tx_ac: bool = False
    ) -> List:
        """Return queried data from tx_exon_aln_v table.

        :param str tx_ac: accession on c. coordinate
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param str alt_ac: accession on g. coordinate
        :param bool use_tx_pos: `True` if querying on transcript position. This means
            `start_pos` and `end_pos` are on the c. coordinate
            `False` if querying on genomic position. This means `start_pos` and
            `end_pos` are on the g. coordinate
        :param bool like_tx_ac: `True` if tx_ac condition should be a like statement.
            This is used when you want to query an accession regardless of its version
            `False` if tx_condition will be exact match
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

        if alt_ac:
            alt_ac_q = f"AND alt_ac = '{alt_ac}'"
        else:
            alt_ac_q = f"AND alt_ac LIKE 'NC_00%'"  # noqa: F541

        if use_tx_pos:
            pos_q = f"""tx_start_i AND tx_end_i"""  # noqa: F541
        else:
            pos_q = f"""alt_start_i AND alt_end_i"""  # noqa: F541

        query = (
            f"""
            SELECT hgnc, tx_ac, tx_start_i, tx_end_i, alt_ac, alt_start_i,
                alt_end_i, alt_strand, alt_aln_method, tx_exon_id, alt_exon_id
            FROM {self.schema}.tx_exon_aln_v
            {tx_q}
            {alt_ac_q}
            {aln_method}
            AND {start_pos} BETWEEN {pos_q}
            AND {end_pos} BETWEEN {pos_q}
            ORDER BY (CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1,
                LENGTH(alt_ac)) AS INT))
            """
        )
        result = await self.execute_query(query)
        if not result:
            logger.warning(f"Unable to find transcript alignment for query: "
                           f"{query}")
            return []
        if alt_ac and not use_tx_pos:
            if len(result) > 1:
                logger.debug(f"Found more than one match for tx_ac {temp_ac} "
                             f"and alt_ac = {alt_ac}")
        results = list()
        for r in result:
            results.append([field for field in r])
        return results

    @staticmethod
    def data_from_result(result: List) -> Optional[Dict]:
        """Return data found from result.

        :param List result: Data from tx_exon_aln_v table
        :return: Gene, strand, and position ranges for tx and alt_ac
        """
        gene = result[0]
        if result[7] == -1:
            strand = "-"
        else:
            strand = "+"
        tx_pos_range = result[2], result[3]
        alt_pos_range = result[5], result[6]
        alt_aln_method = result[8]
        tx_exon_id = result[9]
        alt_exon_id = result[10]

        if (tx_pos_range[1] - tx_pos_range[0]) != \
                (alt_pos_range[1] - alt_pos_range[0]):
            logger.warning(f"tx_pos_range {tx_pos_range} "
                           f"is not the same length as alt_pos_range "
                           f"{alt_pos_range}.")
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

    async def get_mane_c_genomic_data(self, ac: str, alt_ac: Optional[str],
                                      start_pos: int,
                                      end_pos: int) -> Optional[Dict]:
        """Get MANE Transcript and genomic data.

        Used when going from g -> MANE c
        :param str ac: MANE Transcript accession
        :param str alt_ac: NC Accession
        :param int start_pos: Genomic start position change
        :param int end_pos: Genomic end position change
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

        if data["strand"] == "-":
            end_pos += 1
            start_pos += 1
            data["alt_pos_change_range"] = (end_pos, start_pos)
            data["alt_pos_change"] = (
                data["alt_pos_range"][1] - data["alt_pos_change_range"][0],
                data["alt_pos_change_range"][1] - data["alt_pos_range"][0]
            )
        else:
            data["alt_pos_change_range"] = (start_pos, end_pos)
            data["alt_pos_change"] = (
                data["alt_pos_change_range"][0] - data["alt_pos_range"][0],
                data["alt_pos_range"][1] - data["alt_pos_change_range"][1]
            )

        return data

    async def get_genomic_tx_data(
        self, tx_ac: str, pos: Tuple[int, int],
        annotation_layer: Union[AnnotationLayer.CDNA, AnnotationLayer.GENOMIC] = AnnotationLayer.CDNA,  # noqa: E501
        alt_ac: Optional[str] = None,
        target_genome_assembly: Assembly = Assembly.GRCH38
    ) -> Optional[Dict]:
        """Get transcript mapping to genomic data.

        :param str tx_ac: Accession on c. coordinate
        :param Tuple pos: (start pos, end pos)
        :param Union[AnnotationLayer.CDNA, AnnotationLayer.GENOMIC] annotation_layer:
            Annotation layer for `ac` and `pos`
        :param Optional[str] alt_ac: Accession on g. coordinate
        :param Assembly target_genome_assembly: Genome assembly to get genomic data for.
            If `alt_ac` is provided, it will return the associated assembly.
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        results = await self.get_tx_exon_aln_v_data(
            tx_ac, pos[0], pos[1], use_tx_pos=annotation_layer == AnnotationLayer.CDNA,
            alt_ac=alt_ac)
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
            data["tx_pos_range"][1] - pos[1]
        )

        if annotation_layer == AnnotationLayer.CDNA:
            if data["strand"] == "-":
                data["alt_pos_change_range"] = (
                    data["alt_pos_range"][1] - data["pos_change"][0],
                    data["alt_pos_range"][0] + data["pos_change"][1]
                )
            else:
                data["alt_pos_change_range"] = (
                    data["alt_pos_range"][0] + data["pos_change"][0],
                    data["alt_pos_range"][1] - data["pos_change"][1]
                )
        else:
            if data["strand"] == "-":
                data["alt_pos_change_range"] = (pos[1] + 1, pos[0] + 1)
            else:
                data["alt_pos_change_range"] = pos

        return data

    async def get_ac_from_gene(self, gene: str) -> Optional[List[str]]:
        """Return genomic accession for a gene.

        :param str gene: Gene symbol
        :return: List of genomic accessions
        """
        query = (
            f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.genomic
            WHERE hgnc = '{gene}'
            AND alt_ac LIKE 'NC_00%'
            ORDER BY alt_ac DESC;
            """
        )
        results = await self.execute_query(query)
        if not results:
            return []
        return [item for sublist in results for item in sublist]

    async def get_gene_from_ac(self, ac: str, start_pos: int,
                               end_pos: int) -> Optional[List[str]]:
        """Get transcripts from NC accession and positions.

        :param str ac: NC Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :return: List of HGNC gene symbols
        """
        if end_pos is None:
            end_pos = start_pos
        query = (
            f"""
            SELECT DISTINCT hgnc
            FROM {self.schema}.genomic
            WHERE alt_ac = '{ac}'
            AND {start_pos} BETWEEN alt_start_i AND alt_end_i
            AND {end_pos} BETWEEN alt_start_i AND alt_end_i;
            """
        )
        results = await self.execute_query(query)
        if not results:
            logger.warning(f"Unable to find gene between {start_pos} and"
                           f" {end_pos} on {ac}")
            return None
        else:
            if len(results) > 1:
                logger.info(f"Found more than one gene between "
                            f"{start_pos} and {end_pos} on {ac}")

        return [r[0] for r in results]

    async def get_transcripts_from_gene(
        self, gene: str, start_pos: int, end_pos: int, use_tx_pos: bool = True,
        alt_ac: Optional[str] = None
    ) -> pd.core.frame.DataFrame:
        """Get transcripts associated to a gene.

        :param str gene: Gene symbol
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param bool use_tx_pos: `True` if querying on transcript position.
            This means `start_pos` and `end_pos` are c. coordinate positions
            `False` if querying on genomic position. This means `start_pos`
            and `end_pos` are g. coordinate positions
        :param Optional[str] alt_ac: Genomic accession
        :return: Data Frame containing transcripts associated with a gene.
            Transcripts are ordered by most recent NC accession, then by
            descending transcript length.
        """
        if use_tx_pos:
            pos_cond = (
                f"""
                AND {start_pos} + T.cds_start_i
                    BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                AND {end_pos} + T.cds_start_i
                    BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
                """
            )
        else:
            pos_cond = (
                f"""
                AND {start_pos} BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                AND {end_pos} BETWEEN ALIGN.alt_start_i AND ALIGN.alt_end_i
                """
            )

        if alt_ac:
            alt_ac_cond = f"AND ALIGN.alt_ac = '{alt_ac}'"
        else:
            alt_ac_cond = "AND ALIGN.alt_ac LIKE 'NC_00%'"

        query = (
            f"""
            SELECT AA.pro_ac, AA.tx_ac, ALIGN.alt_ac, T.cds_start_i
            FROM {self.schema}.associated_accessions as AA
            JOIN {self.schema}.transcript as T ON T.ac = AA.tx_ac
            JOIN {self.schema}.tx_exon_aln_v as ALIGN ON T.ac = ALIGN.tx_ac
            WHERE T.hgnc = '{gene}'
            {alt_ac_cond}
            AND ALIGN.alt_aln_method = 'splign'
            {pos_cond}
            ORDER BY ALIGN.alt_ac, ALIGN.tx_end_i - ALIGN.tx_start_i DESC;
            """
        )
        results = await self.execute_query(query)
        return pd.DataFrame(
            results,
            columns=["pro_ac", "tx_ac", "alt_ac", "cds_start_i"]
        ).drop_duplicates()

    async def get_chr_assembly(self, ac: str) -> Optional[Tuple[str, str]]:
        """Get chromosome and assembly for NC accession if not in GRCh38.

        :param str ac: NC accession
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
            logger.warning(f"Assembly not supported: {assembly}. "
                           f"Only GRCh37 and GRCh38 are supported.")
            return None

        return chromosome, assembly

    async def liftover_to_38(self, genomic_tx_data: Dict) -> None:
        """Liftover genomic_tx_data to hg38 assembly.

        :param Dict genomic_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        """
        descr = await self.get_chr_assembly(genomic_tx_data["alt_ac"])
        if descr is None:
            # already grch38
            return None
        chromosome, _ = descr

        query = (
            f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac = '{genomic_tx_data['tx_ac']}';
            """
        )
        nc_acs = await self.execute_query(query)
        nc_acs = [nc_ac[0] for nc_ac in nc_acs]
        if nc_acs == [genomic_tx_data["alt_ac"]]:
            logger.warning(f"UTA does not have GRCh38 assembly for "
                           f"{genomic_tx_data['alt_ac'].split('.')[0]}")
            return None

        # Get most recent assembly version position
        # Liftover range
        self._set_liftover(genomic_tx_data, "alt_pos_range", chromosome,
                           Assembly.GRCH38)

        # Liftover changes range
        self._set_liftover(genomic_tx_data, "alt_pos_change_range", chromosome,
                           Assembly.GRCH38)

        # Change alt_ac to most recent
        query = (
            f"""
            SELECT alt_ac
            FROM {self.schema}.genomic
            WHERE alt_ac LIKE '{genomic_tx_data['alt_ac'].split('.')[0]}%'
            ORDER BY (CAST(SUBSTR(alt_ac, position('.' in alt_ac) + 1,
                LENGTH(alt_ac)) AS INT)) DESC;
            """
        )
        nc_acs = await self.execute_query(query)
        genomic_tx_data["alt_ac"] = nc_acs[0][0]

    def get_liftover(self, chromosome: str, pos: int,
                     liftover_to_assembly: Assembly) -> Optional[Tuple]:
        """Get new genome assembly data for a position on a chromosome.

        :param str chromosome: The chromosome number. Must be prefixed with `chr`
        :param int pos: Position on the chromosome
        :param Assembly liftover_to_assembly: Assembly to liftover to
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

    def _set_liftover(self, genomic_tx_data: Dict, key: str, chromosome: str,
                      liftover_to_assembly: Assembly) -> None:
        """Update genomic_tx_data to have coordinates for given assembly.

        :param Dict genomic_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        :param str key: Key to access coordinate positions
        :param str chromosome: Chromosome, must be prefixed with `chr`
        :param Assembly liftover_to_assembly: Assembly to liftover to
        """
        liftover_start_i = self.get_liftover(chromosome, genomic_tx_data[key][0],
                                             liftover_to_assembly)
        if liftover_start_i is None:
            logger.warning(f"Unable to liftover position "
                           f"{genomic_tx_data[key][0]} on {chromosome}")
            return None

        liftover_end_i = self.get_liftover(chromosome, genomic_tx_data[key][1],
                                           liftover_to_assembly)
        if liftover_end_i is None:
            logger.warning(f"Unable to liftover position "
                           f"{genomic_tx_data[key][1]} on {chromosome}")
            return None

        genomic_tx_data[key] = liftover_start_i[1], liftover_end_i[1]

    async def p_to_c_ac(self, p_ac: str) -> List[str]:
        """Return c. accession from p. accession.

        :param str p_ac: Protein accession
        :return: List of rows containing c. accessions that are associated
            with the given p. accession.
        """
        query = (
            f"""
            SELECT tx_ac
            FROM {self.schema}.associated_accessions
            WHERE pro_ac = '{p_ac}'
            ORDER BY tx_ac;
            """
        )
        result = await self.execute_query(query)
        if result:
            result = [r[0] for r in result]
        return result

    async def get_transcripts_from_genomic_pos(
            self, alt_ac: str, g_pos: int) -> List[str]:
        """Get transcripts associated to a genomic ac and position.

        :param str alt_ac: Genomic accession
        :param int g_pos: Genomic position
        :return: RefSeq transcripts on c. coordinate
        """
        query = (
            f"""
               SELECT distinct tx_ac
               FROM {self.schema}.tx_exon_aln_v
               WHERE alt_ac = '{alt_ac}'
               AND {g_pos} BETWEEN alt_start_i AND alt_end_i
               AND tx_ac LIKE 'NM_%';
               """
        )
        results = await self.execute_query(query)
        if not results:
            return []
        return [item for sublist in results for item in sublist]

    @staticmethod
    def get_secret() -> str:
        """Get secrets for UTA DB instances."""
        secret_name = environ["UTA_DB_SECRET"]
        region_name = "us-east-2"

        # Create a Secrets Manager client
        session = boto3.session.Session()
        client = session.client(
            service_name="secretsmanager",
            region_name=region_name
        )

        try:
            get_secret_value_response = client.get_secret_value(
                SecretId=secret_name
            )
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
                    get_secret_value_response["SecretBinary"])
                return decoded_binary_secret


class ParseResult(urlparse.ParseResult):
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
