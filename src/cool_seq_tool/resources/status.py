"""Enable quick status check of Cool-Seq-Tool resources."""

import logging
from collections import namedtuple
from pathlib import Path

from agct._core import ChainfileError
from asyncpg import InvalidCatalogNameError, UndefinedTableError
from biocommons.seqrepo import SeqRepo

from cool_seq_tool.handlers.seqrepo_access import SEQREPO_ROOT_DIR, SeqRepoAccess
from cool_seq_tool.mappers.liftover import LiftOver
from cool_seq_tool.resources.data_files import DataFile, get_data_file
from cool_seq_tool.sources.uta_database import UTA_DB_URL, UtaDatabase

_logger = logging.getLogger(__name__)


ResourceStatus = namedtuple(
    "ResourceStatus",
    (
        "uta",
        "seqrepo",
        DataFile.TRANSCRIPT_MAPPINGS.lower(),
        DataFile.MANE_SUMMARY.lower(),
        DataFile.LRG_REFSEQGENE.lower(),
        "liftover",
    ),
)


async def check_status(
    transcript_file_path: Path | None = None,
    lrg_refseqgene_path: Path | None = None,
    mane_data_path: Path | None = None,
    db_url: str = UTA_DB_URL,
    sr: SeqRepo | None = None,
    chain_file_37_to_38: str | None = None,
    chain_file_38_to_37: str | None = None,
) -> ResourceStatus:
    """Perform basic status checks on availability of required data resources.

    Arguments are intended to mirror arguments to :py:meth:`cool_seq_tool.app.CoolSeqTool.__init__`.

    Additional arguments are available for testing paths to specific chainfiles (same
    signature as :py:meth:`cool_seq_tool.mappers.liftover.LiftOver.__init__`).

    >>> from cool_seq_tool.resources.status import check_status
    >>> await check_status()
    ResourceStatus(uta=True, seqrepo=True, transcript_mappings=True, mane_summary=True, lrg_refseqgene=True, liftover=True)

    :param transcript_file_path: The path to ``transcript_mapping.tsv``
    :param lrg_refseqgene_path: The path to the LRG_RefSeqGene file
    :param mane_data_path: Path to RefSeq MANE summary data
    :param db_url: PostgreSQL connection URL
        Format: ``driver://user:password@host/database/schema``
    :param chain_file_37_to_38: Optional path to chain file for 37 to 38 assembly. This
        is used for ``agct``. If this is not provided, will check to see if
        ``LIFTOVER_CHAIN_37_TO_38`` env var is set. If neither is provided, will allow
        ``agct`` to download a chain file from UCSC
    :param chain_file_38_to_37: Optional path to chain file for 38 to 37 assembly. This
        is used for ``agct``. If this is not provided, will check to see if
        ``LIFTOVER_CHAIN_38_TO_37`` env var is set. If neither is provided, will allow
        ``agct`` to download a chain file from UCSC
    :return: boolean description of availability of each resource, given current
        environment configurations
    """
    file_path_params = {
        DataFile.TRANSCRIPT_MAPPINGS.lower(): transcript_file_path,
        DataFile.LRG_REFSEQGENE.lower(): lrg_refseqgene_path,
        DataFile.MANE_SUMMARY.lower(): mane_data_path,
    }

    status = {
        DataFile.TRANSCRIPT_MAPPINGS.lower(): False,
        DataFile.LRG_REFSEQGENE.lower(): False,
        DataFile.MANE_SUMMARY.lower(): False,
        "liftover": False,
        "uta": False,
        "seqrepo": False,
    }
    for r in list(DataFile):
        name_lower = r.lower()
        declared_path = file_path_params[name_lower]
        if declared_path and declared_path.exists() and declared_path.is_file():
            status[name_lower] = True
            continue
        try:
            get_data_file(r)
        except FileNotFoundError:
            _logger.error(
                "%s does not exist at configured location %s", name_lower, declared_path
            )
        except ValueError:
            _logger.error(
                "%s configured at %s is not a valid file.", name_lower, declared_path
            )
        except Exception as e:
            _logger.critical(
                "Encountered unexpected error fetching %s: %s", name_lower, e
            )
        else:
            status[name_lower] = True

    try:
        LiftOver(
            chain_file_37_to_38=chain_file_37_to_38,
            chain_file_38_to_37=chain_file_38_to_37,
        )
    except (FileNotFoundError, ChainfileError) as e:
        _logger.error("agct converter setup failed: %s", e)
    except Exception as e:
        _logger.critical("Encountered unexpected error setting up agct: %s", e)
    else:
        status["liftover"] = True

    try:
        await UtaDatabase.create(db_url)
    except (OSError, InvalidCatalogNameError, UndefinedTableError) as e:
        _logger.error(
            "Encountered error instantiating UTA at URI %s: %s", UTA_DB_URL, e
        )
    except Exception as e:
        _logger.critical(
            "Encountered unexpected error instantiating UTA from URI %s: %s",
            UTA_DB_URL,
            e,
        )
    else:
        status["uta"] = True

    try:
        if not sr:
            sr = SeqRepo(root_dir=SEQREPO_ROOT_DIR)
        sra = SeqRepoAccess(sr)
        sra.sr["NC_000001.11"][1000:1001]
    except OSError as e:
        _logger.error("Encountered error while instantiating SeqRepo: %s", e)
    except KeyError:
        _logger.error("SeqRepo data fetch test failed -- is it populated?")
    except Exception as e:
        _logger.critical("Encountered unexpected error setting up SeqRepo: %s", e)
    else:
        status["seqrepo"] = True

    structured_status = ResourceStatus(**status)
    if all(status.values()):
        _logger.info("Cool-Seq-Tool resource status passed")
    else:
        _logger.error(
            "Cool-Seq-Tool resource check failed. Result: %s", structured_status
        )
    return structured_status
