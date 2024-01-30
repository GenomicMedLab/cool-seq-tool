"""Module containing default routes"""
import logging
import os
import tempfile
from pathlib import Path

from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from starlette.background import BackgroundTasks

from cool_seq_tool.routers import (
    RESP_DESCR,
    SERVICE_NAME,
    UNHANDLED_EXCEPTION_MSG,
    cool_seq_tool,
)
from cool_seq_tool.schemas import (
    GenomicDataResponse,
    GenomicRequestBody,
    TranscriptRequestBody,
)
from cool_seq_tool.utils import service_meta

logger = logging.getLogger("cool_seq_tool")

router = APIRouter(prefix=f"/{SERVICE_NAME}")


@router.post(
    "/genomic_to_transcript_exon_coordinates",
    summary="Get transcript exon data given genomic coordinate data",
    response_description=RESP_DESCR,
    description="Return transcript exon data",
    response_model=GenomicDataResponse,
)
async def genomic_to_transcript_exon_coordinates(
    request_body: GenomicRequestBody,
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param GenomicRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.model_dump()

    response = GenomicDataResponse(
        genomic_data=None, warnings=[], service_meta=service_meta()
    )

    try:
        response = await cool_seq_tool.ex_g_coords_mapper.genomic_to_transcript_exon_coordinates(
            **request_body
        )
    except Exception as e:
        logger.error("genomic_to_transcript_exon_coordinates unhandled exception %s", e)
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response


@router.post(
    "/transcript_to_genomic_coordinates",
    summary="Get genomic coordinate data given transcript exon data",
    response_description=RESP_DESCR,
    description="Return genomic coordinate data",
    response_model=GenomicDataResponse,
)
async def transcript_to_genomic_coordinates(
    request_body: TranscriptRequestBody,
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param TranscriptRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.model_dump()

    response = GenomicDataResponse(
        genomic_data=None, warnings=[], service_meta=service_meta()
    )

    try:
        response = (
            await cool_seq_tool.ex_g_coords_mapper.transcript_to_genomic_coordinates(
                **request_body
            )
        )
    except Exception as e:
        logger.error("transcript_to_genomic_coordinates unhandled exception %s", e)
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response


@router.get(
    "/download_sequence",
    summary="Get sequence for ID",
    response_description=RESP_DESCR,
    description="Given a known accession identifier, retrieve sequence data and return"
    "as a FASTA file",
    response_class=FileResponse,
)
async def get_sequence(
    background_tasks: BackgroundTasks,
    sequence_id: str = Query(
        ..., description="ID of sequence to retrieve, sans namespace"
    ),
) -> FileResponse:
    """Get sequence for requested sequence ID.
    :param sequence_id: accession ID, sans namespace, eg `NM_152263.3`
    :param background_tasks: Starlette background tasks object. Use to clean up
        tempfile after get method returns.
    :return: FASTA file if successful, or 404 if unable to find matching resource
    """
    _, path = tempfile.mkstemp(suffix=".fasta")
    try:
        cool_seq_tool.seqrepo_access.get_fasta_file(sequence_id, Path(path))
    except KeyError as e:
        raise HTTPException(
            status_code=404, detail="No sequence available for requested identifier"
        ) from e
    background_tasks.add_task(lambda p: os.unlink(p), path)  # noqa: PTH108
    return FileResponse(path)
