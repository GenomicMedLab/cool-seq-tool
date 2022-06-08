"""Main application for FastAPI"""
from enum import Enum
from typing import Dict, List, Optional

from fastapi import FastAPI, Query
from fastapi.openapi.utils import get_openapi

from uta_tools import UTATools, logger
from uta_tools.data_sources.mane_transcript import MANETranscriptError
from uta_tools.schemas import AnnotationLayer, Assembly, GenomicDataResponse,\
    GenomicRequestBody, ManeDataService, MappedManeDataService, ResidueMode,\
    TranscriptRequestBody
from uta_tools.version import __version__


SERVICE_NAME = "uta_tools"


class Tags(str, Enum):
    """Define tags for endpoints"""

    MANE_TRANSCRIPT = "MANE Transcript"


app = FastAPI(
    docs_url=f"/{SERVICE_NAME}",
    openapi_url=f"/{SERVICE_NAME}/openapi.json",
    swagger_ui_parameters={"tryItOutEnabled": True}
)


def custom_openapi() -> Dict:
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The GenomicMedLab UTA Tools",
        version=__version__,
        description="Service for querying data in the biocommons UTA database and retrieving MANE data.",  # noqa: E501
        routes=app.routes
    )

    openapi_schema["info"]["contact"] = {
        "name": "Alex H. Wagner",
        "email": "Alex.Wagner@nationwidechildrens.org",
        "url": "https://www.nationwidechildrens.org/specialties/institute-for-genomic-medicine/research-labs/wagner-lab"  # noqa: E501
    }
    app.openapi_schema = openapi_schema
    return app.openapi_schema


app.openapi = custom_openapi

uta_tools = UTATools()

RESP_DESCR = "A response to a validly-formed query."
UNHANDLED_EXCEPTION_MSG = "Unhandled exception occurred. Check logs for more details."


@app.post(f"/{SERVICE_NAME}/genomic_to_transcript_exon_coordinates",
          summary="Get transcript exon data given genomic coordinate data",
          response_description=RESP_DESCR,
          description="Return transcript exon data",
          response_model=GenomicDataResponse)
async def genomic_to_transcript_exon_coordinates(
    request_body: GenomicRequestBody
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param GenomicRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.dict()

    response = GenomicDataResponse(
        genomic_data=None, warnings=list(), service_meta=uta_tools.service_meta())

    try:
        response = \
            await uta_tools.genomic_to_transcript_exon_coordinates(**request_body)
    except Exception as e:
        logger.error(f"genomic_to_transcript_exon_coordinates unhandled exception {str(e)}")  # noqa: E501
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response


@app.post(f"/{SERVICE_NAME}/transcript_to_genomic_coordinates",
          summary="Get genomic coordinate data given transcript exon data",
          response_description=RESP_DESCR,
          description="Return genomic coordinate data",
          response_model=GenomicDataResponse)
async def transcript_to_genomic_coordinates(
    request_body: TranscriptRequestBody
) -> GenomicDataResponse:
    """Get transcript exon data given genomic coordinate data

    :param TranscriptRequestBody request_body: Request body

    Returns: GenomicDataResponse with data and warnings
    """
    request_body = request_body.dict()

    response = GenomicDataResponse(
        genomic_data=None, warnings=list(), service_meta=uta_tools.service_meta())

    try:
        response = await uta_tools.transcript_to_genomic_coordinates(**request_body)
    except Exception as e:
        logger.error(f"transcript_to_genomic_coordinates unhandled exception {str(e)}")
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response

ref_descr = "Reference at position given during input. When this is set, it will "\
            "ensure that the reference sequences match for the final result."
try_longest_compatible_descr = "`True` if should try longest compatible remaining if"\
                               " mane transcript was not compatible. `False` otherwise."


@app.get(f"/{SERVICE_NAME}/get_mane_data",
         summary="Retrieve MANE data in inter-residue coordinates",
         response_description=RESP_DESCR,
         description="Return MANE Select, MANE Plus Clinical, or Longest Remaining "
                     "Transcript data in inter-residue coordinates. See our docs for "
                     "more information on transcript priority.",
         response_model=ManeDataService,
         tags=[Tags.MANE_TRANSCRIPT])
async def get_mane_data(
    ac: str = Query(..., description="Accession"),
    start_pos: int = Query(..., description="Start position"),
    start_annotation_layer: AnnotationLayer = Query(..., description="Starting annotation layer for query"),  # noqa: E501
    end_pos: Optional[int] = Query(None, description="End position. If not set, will set to `start_pos`."),  # noqa: #501
    gene: Optional[str] = Query(None, description="HGNC gene symbol"),
    ref: Optional[str] = Query(None, description=ref_descr),
    try_longest_compatible: bool = Query(True, description=try_longest_compatible_descr),  # noqa: E501
    residue_mode: ResidueMode = Query(ResidueMode.RESIDUE, description="Residue mode for position(s)")  # noqa: E501
) -> ManeDataService:
    """Return MANE or Longest Compatible Remaining Transcript data on inter-residue
    coordinates

    :param str ac: Accession
    :param int start_pos: Start position
    :param AnnotationLayer start_annotation_layer: Starting annotation layer for query
    :param Optional[int] end_pos: End position. If `None` assumes
        both  `start_pos` and `end_pos` have same values.
    :param Optional[str] gene: Gene symbol
    :param Optional[str] ref: Reference at position given during input
    :param bool try_longest_compatible: `True` if should try longest
        compatible remaining if mane transcript was not compatible.
        `False` otherwise.
    :param ResidueMode residue_mode: Starting residue mode for `start_pos`
        and `end_pos`. Will always return coordinates in inter-residue
    """
    warnings = list()
    mane_data = None
    try:
        mane_data = await uta_tools.mane_transcript.get_mane_transcript(
            ac=ac, start_pos=start_pos, start_annotation_layer=start_annotation_layer,
            end_pos=end_pos, gene=gene, ref=ref,
            try_longest_compatible=try_longest_compatible, residue_mode=residue_mode)

        if not mane_data:
            warnings.append("Unable to retrieve MANE data")
    except Exception as e:
        logger.exception(f"get_mane_data unhandled exception {e}")
        warnings.append(UNHANDLED_EXCEPTION_MSG)

    return ManeDataService(
        mane_data=mane_data,
        warnings=warnings,
        service_meta=uta_tools.service_meta()
    )


@app.get(f"/{SERVICE_NAME}/get_mapped_mane_data",
         summary="Retrieve MANE Transcript mapped to a given assembly",
         response_description=RESP_DESCR,
         description="Return mapped MANE Transcript data to a given assembly",
         response_model=MappedManeDataService,
         tags=[Tags.MANE_TRANSCRIPT])
async def get_mapped_mane_data(
    gene: str = Query(..., description="HGNC Symbol or Identifier"),
    assembly: Assembly = Query(..., description="Genomic assembly to use"),
    genomic_position: int = Query(..., description="Genomic position associated to the given gene and assembly"),  # noqa: E501
    residue_mode: ResidueMode = Query(ResidueMode.INTER_RESIDUE,
                                      description="Residue mode for `genomic_position`")
) -> MappedManeDataService:
    """Get MANE data for gene, assembly, and position. If GRCh37 assembly is given,
    will return mapped MANE data.

    :param str gene: HGNC symbol or identifier
    :param Assembly assembly: Assembly for the provided genomic position
    :param int genomic_position: Position on the genomic reference sequence to find
        MANE data for
    :param ResidueMode residue_mode: Starting residue mode for `start_pos`
        and `end_pos`. Will always return coordinates in inter-residue
    :return: Mapped MANE or Longest Compatible Remaining data
    """
    warnings: List = list()
    mapped_mane_data = None
    try:
        mapped_mane_data = await uta_tools.mane_transcript.get_mapped_mane_data(
            gene, assembly, genomic_position, residue_mode)
        if not mapped_mane_data:
            warnings.append(f"Unable to find mapped data for gene {gene} at position "
                            f"{genomic_position} ({residue_mode} coordinates) on "
                            f"assembly {assembly}")
    except MANETranscriptError as e:
        e = str(e)
        logger.exception(e)
        warnings.append(e)
    except Exception as e:
        logger.exception(f"get_mapped_mane_data unhandled exception {e}")
        warnings.append(UNHANDLED_EXCEPTION_MSG)

    return MappedManeDataService(
        mapped_mane_data=mapped_mane_data,
        warnings=warnings,
        service_meta=uta_tools.service_meta()
    )
