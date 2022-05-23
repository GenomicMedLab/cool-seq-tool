"""Main application for FastAPI"""
from typing import Dict

from fastapi import FastAPI
from fastapi.openapi.utils import get_openapi

from uta_tools import UTATools, logger
from uta_tools.schemas import GenomicDataResponse, GenomicRequestBody, \
    TranscriptRequestBody
from uta_tools.version import __version__


SERVICE_NAME = "uta_tools"


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
    except TypeError as e:
        error = str(e)
        logger.error(f"genomic_to_transcript_exon_coordinates TypeError: {error}")
        response.warnings.append(error)
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
    except TypeError as e:
        error = str(e)
        logger.error(f"transcript_to_genomic_coordinates TypeError: {error}")
        response.warnings.append(error)
    except Exception as e:
        logger.error(f"transcript_to_genomic_coordinates unhandled exception {str(e)}")  # noqa: E501
        response.warnings.append(UNHANDLED_EXCEPTION_MSG)

    return response
