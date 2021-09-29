"""Module for UTA queries."""
import asyncpg
from uta_tools import UTA_DB_URL, logger
from six.moves.urllib import parse as urlparse
from asyncpg.exceptions import InterfaceError
from typing import Dict, List, Optional, Tuple
from os import environ
from urllib.parse import quote, unquote
import boto3


class UTADatabase:
    """Class for connecting and querying UTA database."""

    def __init__(self, db_url: str = UTA_DB_URL, db_pwd: str = '',
                 is_prod_env: bool = False) -> None:
        """Initialize DB class.

        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :param bool is_prod_env: `True` if working in prod environment.
            `False` otherwise.
        """
        self.schema = None
        self.args = self._get_conn_args(is_prod_env, db_url, db_pwd)
        self._connection_pool = None

    @staticmethod
    def _update_db_url(db_pwd: str, db_url: str) -> Optional[str]:
        """Return new db_url containing password.

        :param str db_pwd: User's password for uta database
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :return: PostgreSQL connection URL
        """
        if 'UTA_DB_URL' in environ:
            return environ["UTA_DB_URL"]
        if not db_pwd and 'UTA_PASSWORD' not in environ:
            raise Exception('Environment variable UTA_PASSWORD '
                            'or `db_pwd` param must be set')
        else:
            uta_password_in_environ = 'UTA_PASSWORD' in environ
            db_url = db_url.split('@')
            if uta_password_in_environ and db_pwd:
                if db_pwd != environ['UTA_PASSWORD']:
                    raise Exception('If both environment variable UTA_PASSWORD'
                                    ' and param db_pwd is set, they must '
                                    'both be the same')
            else:
                if uta_password_in_environ and not db_pwd:
                    db_pwd = environ['UTA_PASSWORD']
            return f"{db_url[0]}:{db_pwd}@{db_url[1]}"

    def _get_args(self, host: str, port: int, database: str, user: str,
                  password: str) -> Dict:
        """Return db credentials.

        :param str host: Database host name
        :param int port: Database port number
        :param str database: Database name
        :param str user: Database username
        :param str password: Database password for user
        :return: Database credentials
        """
        return dict(
            host=host,
            port=port,
            database=database,
            user=user,
            password=password,
            application_name='uta_tools'
        )

    def _get_conn_args(self, is_prod_env: bool, db_url: str,
                       db_pwd: str = '') -> Dict:
        """Return connection arguments.

        :param bool is_prod_env: `True` if production environment.
            `False` otherwise.
        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :return: Database credentials
        """
        if is_prod_env:
            self.schema = environ['UTA_SCHEMA']
            region = 'us-east-2'
            client = boto3.client('rds', region_name=region)
            token = client.generate_db_auth_token(
                DBHostname=environ['UTA_HOST'], Port=environ['UTA_PORT'],
                DBUsername=environ['UTA_USER'], Region=region
            )
            return self._get_args(environ['UTA_HOST'],
                                  int(environ['UTA_PORT']),
                                  environ['UTA_DATABASE'], environ['UTA_USER'],
                                  token)
        else:
            db_url = self._update_db_url(db_pwd, db_url)
            original_pwd = db_url.split('//')[-1].split('@')[0].split(':')[-1]
            db_url = db_url.replace(original_pwd, quote(original_pwd))

            url = ParseResult(urlparse.urlparse(db_url))
            self.schema = url.schema
            password = unquote(url.password) if url.password else ''
            return self._get_args(url.hostname, url.port, url.database,
                                  url.username, password)

    async def create_pool(self) -> None:
        """Create connection pool if not already created."""
        if not self._connection_pool:
            try:
                self._connection_pool = await asyncpg.create_pool(
                    min_size=1,
                    max_size=10,
                    max_inactive_connection_lifetime=3,
                    command_timeout=60,
                    host=self.args['host'],
                    port=self.args['port'],
                    user=self.args['user'],
                    password=self.args['password'],
                    database=self.args['database'],
                )
            except InterfaceError as e:
                logger.error(f'While creating connection pool, '
                             f'encountered exception {e}')

    async def execute_query(self, query: str) -> any:
        """Execute a query and return its result.

        :param str query: Query to make on database
        :return: Query's result
        """
        if not self._connection_pool:
            await self.create_pool()
        async with self._connection_pool.acquire() as connection:
            async with connection.transaction():
                result = await connection.fetch(query)
                return result

    async def transcript_to_genomic(
            self, tx_ac: str, exon_start: int, exon_end: int,
            exon_start_offset: int = 0, exon_end_offset: int = 0,
            gene: str = None) -> Optional[Dict]:
        """Get genomic data given transcript data.

        :param str tx_ac: Transcript accession
        :param int exon_start: Starting exon number
        :param int exon_end: Ending exon number
        :param int exon_start_offset: Starting exon offset
        :param int exon_end_offset: Ending exon offset
        :param str gene: Gene symbol
        :return: Dictionary containing transcript exon data and genomic
            start/end coordinates, or None if lookup fails
        """
        if gene:
            gene = gene.upper().strip()

        if tx_ac:
            tx_ac = tx_ac.strip()

        tx_exon_start_end = await self.get_tx_exon_start_end(
            tx_ac, exon_start, exon_end)
        if not tx_exon_start_end:
            return None
        (tx_exons, exon_start, exon_end) = tx_exon_start_end

        tx_exon_coords = self.get_tx_exon_coords(
            tx_exons, exon_start, exon_end)
        if not tx_exon_coords:
            return None
        tx_exon_start, tx_exon_end = tx_exon_coords

        alt_ac_start_end = await self.get_alt_ac_start_and_end(
            tx_ac, tx_exon_start, tx_exon_end, gene=gene)
        if not alt_ac_start_end:
            return None
        alt_ac_start, alt_ac_end = alt_ac_start_end

        start = alt_ac_start[3]
        end = alt_ac_end[2]
        strand = alt_ac_start[4]
        if strand == -1:
            start_offset = exon_start_offset * -1
            end_offset = exon_end_offset * -1
        else:
            start_offset = exon_start_offset
            end_offset = exon_end_offset
        start += start_offset
        end += end_offset

        gene = alt_ac_start[0]
        chr = alt_ac_start[1]
        if gene is None or chr is None:
            return None
        return {
            "gene": gene,
            "chr": chr,
            "start": start,
            "end": end,
            "exon_start": exon_start,
            "exon_start_offset": exon_start_offset,
            "exon_end": exon_end,
            "exon_end_offset": exon_end_offset,
        }

    async def get_tx_exons(self, tx_ac: str) -> Optional[List[str]]:
        """Get list of transcript exons start/end coordinates.

        :param str tx_ac: Transcript accession
        :return: List of a transcript's accessions, or None if lookup fails
        """
        query = (
            f"""
            SELECT cds_se_i
            FROM {self.schema}._cds_exons_fp_v
            WHERE tx_ac = '{tx_ac}';
            """
        )
        cds_se_i = await self.execute_query(query)
        if not cds_se_i:
            logger.warning(f"Unable to get exons for {tx_ac}")
            return None
        return cds_se_i[0][0].split(';')

    async def get_tx_exon_start_end(self, tx_ac: str, exon_start: int,
                                    exon_end: int)\
            -> Optional[Tuple[List[str], int, int]]:
        """Get exon start/end coordinates given accession and gene.

        :param str tx_ac: Transcript accession
        :param int exon_start: Starting exon number
        :param int exon_end: Ending exon number
        :return: Transcript's exons and start/end exon coordinates, or
            None if lookup fails
        """
        if exon_start and exon_end:
            if exon_start > exon_end:
                logger.warning(f"start exon, {exon_start},"
                               f"is greater than end exon, {exon_end}")
                return None
            elif exon_end < exon_start:
                logger.warning(f"end exon, {exon_end}, "
                               f"is less than start exon, {exon_start}")
                return None

        tx_exons = await self.get_tx_exons(tx_ac)
        if not tx_exons:
            return None

        if exon_start == 0:
            exon_start = 1

        if exon_end == 0:
            exon_end = len(tx_exons)

        return tx_exons, exon_start, exon_end

    @staticmethod
    def get_tx_exon_coords(tx_exons: List[str], exon_start: int,
                           exon_end: int) -> Optional[Tuple[List, List]]:
        """Get transcript exon coordinates.

        :param list tx_exons: List of transcript exons
        :param int exon_start: Start exon number
        :param int exon_end: End exon number
        :return: [Transcript start exon coords, Transcript end exon coords], or
            None if exon start/end is not valid
        """
        try:
            tx_exon_start = tx_exons[exon_start - 1].split(',')
            tx_exon_end = tx_exons[exon_end - 1].split(',')
        except IndexError as e:
            logger.warning(e)
            return None
        return tx_exon_start, tx_exon_end

    async def get_alt_ac_start_and_end(self, tx_ac: str,
                                       tx_exon_start: List[str],
                                       tx_exon_end: List[str],
                                       gene: str = None)\
            -> Optional[Tuple[Tuple, Tuple]]:
        """Get genomic coordinates for related transcript exon start and end.

        :param str tx_ac: Transcript accession
        :param List[str] tx_exon_start: Transcript's exon start coordinates
        :param List[str] tx_exon_end: Transcript's exon end coordinates
        :param str gene: Gene symbol
        :return: Alt ac start and end data
        """
        alt_ac_start = await self.get_alt_ac_start_or_end(
            tx_ac, int(tx_exon_start[0]), int(tx_exon_start[1]), gene=gene)
        if not alt_ac_start:
            return None

        alt_ac_end = await self.get_alt_ac_start_or_end(
            tx_ac, int(tx_exon_end[0]), int(tx_exon_end[1]), gene=gene)
        if not alt_ac_end:
            return None

        # validate
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
        return alt_ac_start, alt_ac_end

    async def get_alt_ac_start_or_end(self, tx_ac: str, tx_exon_start: int,
                                      tx_exon_end: int, gene: Optional[str])\
            -> Optional[Tuple[str, str, int, int, int]]:
        """Get genomic data for related transcript exon start or end.

        :param str tx_ac: Transcript accession
        :param int tx_exon_start: Transcript's exon start coordinate
        :param int tx_exon_end: Transcript's exon end coordinate
        :param Optional[str] gene: Gene symbol
        :return: hgnc symbol, genomic accession for chromosome,
            start exon's end coordinate, end exon's start coordinate, strand
        """
        if gene:
            gene_query = f"AND T.hgnc = '{gene}'"
        else:
            gene_query = ''

        query = (
            f"""
            SELECT T.hgnc, T.alt_ac, T.alt_start_i, T.alt_end_i, T.alt_strand
            FROM uta_20210129._cds_exons_fp_v as C
            JOIN uta_20210129.tx_exon_aln_v as T ON T.tx_ac = C.tx_ac
            WHERE T.tx_ac = '{tx_ac}'
            {gene_query}
            AND {tx_exon_start} BETWEEN T.tx_start_i AND T.tx_end_i
            AND {tx_exon_end} BETWEEN T.tx_start_i AND T.tx_end_i
            AND T.alt_aln_method = 'splign'
            AND T.alt_ac LIKE 'NC_00%'
            ORDER BY T.alt_ac DESC;
            """
        )
        results = await self.execute_query(query)
        if not results:
            logger.warning(f"Unable to get genomic data for {tx_ac}"
                           f" on start exon {tx_exon_start} and "
                           f"end exon {tx_exon_end}")
            return None
        result = results[0]
        return result[0], result[1], result[2], result[3], result[4]


class ParseResult(urlparse.ParseResult):
    """Subclass of url.ParseResult that adds database and schema methods,
    and provides stringification.
    Source: https://github.com/biocommons/hgvs
    """

    def __new__(cls, pr):
        """Create new instance."""
        return super(ParseResult, cls).__new__(cls, *pr)

    @property
    def database(self):
        """Create database property."""
        path_elems = self.path.split("/")
        return path_elems[1] if len(path_elems) > 1 else None

    @property
    def schema(self):
        """Create schema property."""
        path_elems = self.path.split("/")
        return path_elems[2] if len(path_elems) > 2 else None
