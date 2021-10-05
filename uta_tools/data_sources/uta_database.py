"""Module for UTA queries."""
import asyncpg
import boto3
import pandas as pd
from uta_tools import UTA_DB_URL, logger, IS_PROD_ENV
from six.moves.urllib import parse as urlparse
from asyncpg.exceptions import InterfaceError
from typing import Dict, List, Optional, Tuple
from os import environ
from urllib.parse import quote, unquote
from pyliftover import LiftOver


class UTADatabase:
    """Class for connecting and querying UTA database."""

    def __init__(self, db_url: str = UTA_DB_URL, db_pwd: str = '',
                 liftover_from: str = 'hg19',
                 liftover_to: str = 'hg38') -> None:
        """Initialize DB class.
        After initializing, you must run `_create_genomic_table()`

        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :param str liftover_from: Assembly to liftover from
        :param str liftover_to: Assembly to liftover to
        """
        self.schema = None
        self.args = self._get_conn_args(db_url, db_pwd)
        self._connection_pool = None
        self.liftover = LiftOver(liftover_from, liftover_to)

    @staticmethod
    def _update_db_url(db_pwd: str, db_url: str) -> str:
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
            if uta_password_in_environ and db_pwd:
                if db_pwd != environ['UTA_PASSWORD']:
                    raise Exception('If both environment variable UTA_PASSWORD'
                                    ' and param db_pwd is set, they must '
                                    'both be the same')
            else:
                if uta_password_in_environ and not db_pwd:
                    db_pwd = environ['UTA_PASSWORD']
            db_url = db_url.split('@')
            return f"{db_url[0]}:{db_pwd}@{db_url[1]}"

    @staticmethod
    def _get_args(host: str, port: int, database: str, user: str,
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

    def _get_conn_args(self, db_url: str, db_pwd: str = '') -> Dict:
        """Return connection arguments.

        :param str db_url: PostgreSQL connection URL
            Format: `driver://user:pass@host/database/schema`
        :param str db_pwd: User's password for uta database
        :return: Database credentials
        """
        if IS_PROD_ENV:
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
                raise Exception("Could not create connection pool")

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

    async def _create_genomic_table(self):
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
        genomic_table_exists = genomic_table_exists[0]
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
                f"""CREATE INDEX gene_alt_index ON {self.schema}.genomic (hgnc, alt_ac);"""  # noqa: E501
                f"""CREATE INDEX alt_ac_index ON {self.schema}.genomic (alt_ac);"""  # noqa: E501
            ]
            for create_index in indexes:
                await self.execute_query(create_index)

    def _transform_list(self, li):
        """Transform list to only contain field values"""
        results = list()
        for item in li:
            results.append([field for field in item])
        return results

    async def chr_to_accession(self, chromosome, pos, strand=None,
                               alt_ac=None):
        """Chromosome number to accession.

        :return: List of transcripts associated with input arguments
        """
        if alt_ac:
            alt_ac_cond = f"WHERE alt_ac = '{alt_ac}'"
        else:
            alt_ac_cond = f"WHERE alt_ac ~ '^NC_[0-9]+0{chromosome}.[0-9]+$'"

        if strand:
            strand_cond = f"AND alt_strand = '{strand}'"
        else:
            strand_cond = ""
        query = (
            f"""
            SELECT hgnc, alt_ac
            FROM {self.schema}.tx_exon_aln_v
            {alt_ac_cond}
            AND alt_aln_method = 'splign'
            AND {pos} BETWEEN alt_start_i AND alt_end_i
            {strand_cond}
            """
        )
        results = await self.execute_query(query)
        results = self._transform_list(results)
        genes = set()
        alt_acs = set()
        for r in results:
            genes.add(r[0])
            alt_acs.add(r[1])
        return {
            "genes": genes,
            "alt_acs": alt_acs
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
        # if exon_start or exon_end is 0, we will default to the transcript's
        # first and last transcript
        if exon_start and exon_end:
            if exon_start > exon_end:
                logger.warning(f"start exon, {exon_start},"
                               f"is greater than end exon, {exon_end}")
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
        result = result[0]
        if not result:
            logger.warning(f"Unable to get genomic data for {tx_ac}"
                           f" on start exon {tx_exon_start} and "
                           f"end exon {tx_exon_end}")
            return None
        return result[0], result[1], result[2], result[3], result[4]

    async def get_cds_start_end(self, tx_ac: str) -> Optional[tuple[int, int]]:
        """Get coding start and end site

        :param str tx_ac: Transcript accession
        :return: [Coding start site, Coding end site]
        """
        if tx_ac.startswith('ENS'):
            tx_ac = tx_ac.split('.')[0]
        query = (
            f"""
            SELECT cds_start_i, cds_end_i
            FROM {self.schema}.transcript
            WHERE ac='{tx_ac}'
            """
        )
        cds_start_end = await self.execute_query(query)
        if cds_start_end:
            cds_start_end = cds_start_end[0]
            if cds_start_end[0] is not None \
                    and cds_start_end[1] is not None:
                return cds_start_end
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
            ORDER BY ac
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
            )
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
            WHERE ac = '{ac}'
            """
        )
        result = await self.execute_query(query)
        if not result:
            logger.warning(f"Accession {ac} does not have a description")
            return None
        else:
            result = result[0][0]
            if result == '':
                result = None
            return result

    async def get_tx_exon_aln_v_data(self, ac: str, start_pos: int,
                                     end_pos: int, alt_ac: str = None,
                                     use_tx_pos: bool = True) -> Optional[List]:  # noqa: E501
        """Return queried data from tx_exon_aln_v table.

        :param str ac: Accession
        :param int start_pos: Start position change
        :param int end_pos: End position change
        :param str alt_ac: NC accession
        :param bool use_tx_pos: `True` if querying on transcript position.
            `False` if querying on genomic position.
        :return: tx_exon_aln_v data
        """
        if end_pos is None:
            end_pos = start_pos

        if ac.startswith('ENST'):
            temp_ac = ac.split('.')[0]
            aln_method = f"AND alt_aln_method='genebuild'"  # noqa: F541
        else:
            temp_ac = ac
            aln_method = f"AND alt_aln_method='splign'"  # noqa: F541

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
            WHERE tx_ac='{temp_ac}'
            {alt_ac_q}
            {aln_method}
            AND {start_pos} BETWEEN {pos_q}
            AND {end_pos} BETWEEN {pos_q}
            ORDER BY alt_ac;
            """
        )
        result = await self.execute_query(query)
        if not result:
            logger.warning(f"Unable to find transcript alignment for query: "
                           f"{query}")
            return None
        if alt_ac and not use_tx_pos:
            if len(result) > 1:
                logger.debug(f"Found more than one match for tx_ac {temp_ac} "
                             f"and alt_ac = {alt_ac}")
        results = list()
        for r in result:
            results.append([field for field in r])
        return results

    @staticmethod
    def data_from_result(result: list) -> Optional[Dict]:
        """Return data found from result.

        :param list result: Data from tx_exon_aln_v table
        :return: Gene, strand, and position ranges for tx and alt_ac
        """
        gene = result[0]
        if result[7] == -1:
            strand = '-'
        else:
            strand = '+'
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

    async def get_mane_c_genomic_data(self, ac: str, alt_ac: str,
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

        data['tx_ac'] = result[1]
        data['alt_ac'] = result[4]
        data['coding_start_site'] = coding_start_site[0]
        data['coding_end_site'] = coding_start_site[1]
        data['alt_pos_change'] = (
            start_pos - data['alt_pos_range'][0],
            data['alt_pos_range'][1] - end_pos
        )
        data['alt_pos_change_range'] = (
            data['alt_pos_range'][0] + data['alt_pos_change'][0],
            data['alt_pos_range'][1] - data['alt_pos_change'][1]
        )
        return data

    async def get_genomic_tx_data(self, ac: str,
                                  pos: tuple[int, int]) -> Optional[Dict]:
        """Get transcript mapping to genomic data.

        Used when going from c -> g
        :param str ac: cDNA transcript
        :param tuple pos: [cDNA pos start, cDNA pos end]
        :return: Gene, Transcript accession and position change,
            Altered transcript accession and position change, Strand
        """
        results = await self.get_tx_exon_aln_v_data(ac, pos[0], pos[1])
        if not results:
            return None
        result = results[-1]

        data = self.data_from_result(result)
        if not data:
            return None
        data['tx_ac'] = ac
        data['alt_ac'] = result[4]
        data['pos_change'] = (
            pos[0] - data['tx_pos_range'][0],
            data['tx_pos_range'][1] - pos[1]
        )
        data['alt_pos_change_range'] = (
            data['alt_pos_range'][0] + data['pos_change'][0],
            data['alt_pos_range'][1] - data['pos_change'][1]
        )
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
            ORDER BY alt_ac DESC
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
            AND {end_pos} BETWEEN alt_start_i AND alt_end_i
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

    async def get_transcripts_from_gene(self, gene: str, start_pos: int,
                                        end_pos: int) -> pd.core.frame.DataFrame:  # noqa: E501
        """Get transcripts on c coordinate associated to a gene.

        :param str gene: Gene symbol
        :param int start_pos: Start position change on c. coordinate
        :param int end_pos: End position change on c. coordinate
        :return: Data Frame containing transcripts associated with a gene.
            Transcripts are ordered by most recent NC accession, then by
            descending transcript length.
        """
        query = (
            f"""
            SELECT AA.pro_ac, AA.tx_ac, ALIGN.alt_ac, T.cds_start_i
            FROM {self.schema}.associated_accessions as AA
            JOIN {self.schema}.transcript as T ON T.ac = AA.tx_ac
            JOIN {self.schema}.tx_exon_aln_v as ALIGN ON T.ac = ALIGN.tx_ac
            WHERE T.hgnc = '{gene}'
            AND ALIGN.alt_ac LIKE 'NC_00%'
            AND ALIGN.alt_aln_method = 'splign'
            AND {start_pos} + T.cds_start_i
                BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
            AND {end_pos} + T.cds_start_i
                BETWEEN ALIGN.tx_start_i AND ALIGN.tx_end_i
            ORDER BY ALIGN.alt_ac, ALIGN.tx_end_i - ALIGN.tx_start_i DESC;
            """
        )
        results = await self.execute_query(query)
        return pd.DataFrame(
            results, columns=["pro_ac", "tx_ac", "alt_ac", "cds_start_i"])

    async def get_chr_assembly(self, ac: str) -> Optional[Tuple[str, str]]:
        """Get chromosome and assembly for NC accession if not in GRCh38.

        :param str ac: NC accession
        :return: Chromosome and Assembly accession is on
        """
        descr = await self.get_ac_descr(ac)
        if not descr:
            # Already GRCh38 Assembly
            return None
        descr = descr.split(',')
        chromosome = f"chr{descr[0].split()[-1]}"
        assembly = f"GRCh{descr[1].split('.')[0].split('GRCh')[-1]}"

        if assembly not in ['GRCh37', 'GRCh38']:
            logger.warning(f"Assembly not supported: {assembly}. "
                           f"Only GRCh37 and GRCh38 are supported.")
            return None

        return chromosome, assembly

    async def liftover_to_38(self, genomic_tx_data: dict) -> None:
        """Liftover genomic_tx_data to hg38 assembly.

        :param dict genomic_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        """
        descr = await self.get_chr_assembly(genomic_tx_data['alt_ac'])
        if descr is None:
            return None
        chromosome, assembly = descr

        query = (
            f"""
            SELECT DISTINCT alt_ac
            FROM {self.schema}.tx_exon_aln_v
            WHERE tx_ac = '{genomic_tx_data['tx_ac']}'
            """
        )
        nc_acs = await self.execute_query(query)
        if len(nc_acs) == 1:
            logger.warning(f"UTA does not have GRCh38 assembly for "
                           f"{genomic_tx_data['alt_ac'].split('.')[0]}")
            return None

        # Get most recent assembly version position
        # Liftover range
        self._set_liftover(
            genomic_tx_data, 'alt_pos_range', chromosome
        )

        # Liftover changes range
        self._set_liftover(
            genomic_tx_data, 'alt_pos_change_range', chromosome
        )

        # Change alt_ac to most recent
        query = (
            f"""
            SELECT alt_ac
            FROM {self.schema}.genomic
            WHERE alt_ac LIKE '{genomic_tx_data['alt_ac'].split('.')[0]}%'
            ORDER BY alt_ac;
            """
        )
        nc_acs = await self.execute_query(query)
        genomic_tx_data['alt_ac'] = nc_acs[-1][0]

    def get_liftover(self, chromosome: str, pos: int) -> Optional[Tuple]:
        """Get new genome assembly data for a position on a chromosome.

        :param str chromosome: The chromosome number
        :param int pos: Position on the chromosome
        :return: [Target chromosome, target position, target strand,
            conversion_chain_score] for hg38 assembly
        """
        liftover = self.liftover.convert_coordinate(chromosome, pos)
        if liftover is None or len(liftover) == 0:
            logger.warning(f"{pos} does not exist on {chromosome}")
            return None
        else:
            return liftover[0]

    def _set_liftover(self, genomic_tx_data: dict, key: str,
                      chromosome: str) -> None:
        """Update genomic_tx_data to have hg38 coordinates.

        :param dict genomic_tx_data: Dictionary containing gene, nc_accession,
            alt_pos, and strand
        :param str key: Key to access coordinate positions
        :param str chromosome: Chromosome
        """
        liftover_start_i = self.get_liftover(chromosome,
                                             genomic_tx_data[key][0])
        if liftover_start_i is None:
            logger.warning(f"Unable to liftover position "
                           f"{genomic_tx_data[key][0]} on {chromosome}")
            return None

        liftover_end_i = self.get_liftover(chromosome,
                                           genomic_tx_data[key][1])
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
