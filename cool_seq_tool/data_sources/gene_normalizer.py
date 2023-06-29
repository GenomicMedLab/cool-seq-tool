"""Module for accessing Gene Normalizer"""
import logging
from typing import Dict, Optional

from gene.database.dynamodb import DynamoDbDatabase
from gene.query import QueryHandler
from gene.schemas import SourceName


logger = logging.getLogger("cool_seq_tool")


class GeneNormalizer:
    """Gene Normalizer class for getting gene data"""

    def __init__(
        self, query_handler: Optional[QueryHandler] = None, db_url: str = "",
        db_region: str = "us-east-2"
    ) -> None:
        """Initialize gene normalizer class

        :param query_handler: Gene normalizer query handler instance that will be
            reused by cool_seq_tool. If left empty or given ``None``, will create a
            new instance instead.
        :param db_url: URL to gene normalizer dynamodb. Ignored unless ``query_handler``
            is ``None``.
        :param db_region: AWS region for gene normalizer db. Ignored unless
            ``query_handler`` is ``None``.
        """
        if query_handler:
            self.query_handler = query_handler
        else:
            ddb = DynamoDbDatabase(db_url=db_url, region_name=db_region)
            self.query_handler = QueryHandler(ddb)

    def get_hgnc_data(self, gene: str) -> Optional[Dict]:
        """Return HGNC data for a given gene

        :param gene: Gene query
        :return: HGNC data corresponding to normalized concept matching given gene term,
            if available
        """
        hgnc_data = dict()
        gene_resp = self.query_handler.normalize_unmerged(gene)
        hgnc_matches = gene_resp.source_matches.get(SourceName.HGNC)
        if hgnc_matches and hgnc_matches.records:
            hgnc_data = hgnc_matches.records[0].dict()
        else:
            logger.warning(f"Unable to get HGNC symbol for {gene}")
        return hgnc_data
