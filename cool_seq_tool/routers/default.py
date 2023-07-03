"""Module containing default routes"""
import logging
import os
import tempfile
from pathlib import Path

from fastapi import APIRouter
from fastapi import Query, HTTPException
from fastapi.responses import FileResponse
from starlette.background import BackgroundTasks


from cool_seq_tool.routers import cool_seq_tool, SERVICE_NAME, RESP_DESCR, \
    UNHANDLED_EXCEPTION_MSG
from cool_seq_tool.schemas import GenomicDataResponse, GenomicRequestBody, \
    TranscriptRequestBody


logger = logging.getLogger("cool_seq_tool")

router = APIRouter(prefix=f"/{SERVICE_NAME}")


@router.post(
    "/genomic_to_transcript_exon_coordinates",
    summary="Get transcript exon data given genomic coordinate data",
    response_description=RESP_DESCR,
    description="Return transcript exon data",
    response_model=GenomicDataResponse
)
async def genomic_to_transcript_exon_coordinates(
    request_body: GenomicRequestBody
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param GenomicRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.dict()

    response = GenomicDataResponse(
        genomic_data=None, warnings=list(), service_meta=cool_seq_tool.service_meta())

    try:
        response = \
            await cool_seq_tool.genomic_to_transcript_exon_coordinates(**request_body)
    except Exception as e:
        logger.error(f"genomic_to_transcript_exon_coordinates unhandled exception {str(e)}")  # noqa: E501
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response


@router.post(
    "/transcript_to_genomic_coordinates",
    summary="Get genomic coordinate data given transcript exon data",
    response_description=RESP_DESCR,
    description="Return genomic coordinate data",
    response_model=GenomicDataResponse
)
async def transcript_to_genomic_coordinates(
    request_body: TranscriptRequestBody
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param TranscriptRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.dict()

    response = GenomicDataResponse(
        genomic_data=None, warnings=list(), service_meta=cool_seq_tool.service_meta())

    try:
        response = await cool_seq_tool.transcript_to_genomic_coordinates(**request_body)
    except Exception as e:
        logger.error(f"transcript_to_genomic_coordinates unhandled exception {str(e)}")
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response


@router.get(
    "/download_sequence",
    summary="Get sequence for ID",
    response_description=RESP_DESCR,
    description="Given a known accession identifier, retrieve sequence data and return"
                "as a FASTA file",
    response_class=FileResponse
)
async def get_sequence(
    background_tasks: BackgroundTasks,
    sequence_id: str = Query(
        ...,
        description="ID of sequence to retrieve, sans namespace"
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
        cool_seq_tool.get_fasta_file(sequence_id, Path(path))
    except KeyError:
        raise HTTPException(
            status_code=404,
            detail="No sequence available for requested identifier"
        )
    background_tasks.add_task(
        lambda p: os.unlink(p),
        path
    )
    return FileResponse(path)
