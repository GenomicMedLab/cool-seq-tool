"""Main application for FastAPI"""
from typing import Dict

from fastapi import FastAPI
from fastapi.openapi.utils import get_openapi


from cool_seq_tool.routers import default, mane, mappings, SERVICE_NAME
from cool_seq_tool.version import __version__


app = FastAPI(
    docs_url=f"/{SERVICE_NAME}",
    openapi_url=f"/{SERVICE_NAME}/openapi.json",
    swagger_ui_parameters={"tryItOutEnabled": True}
)


app.include_router(default.router)
app.include_router(mane.router)
app.include_router(mappings.router)


def custom_openapi() -> Dict:
    """Generate custom fields for OpenAPI response."""
    if app.openapi_schema:
        return app.openapi_schema
    openapi_schema = get_openapi(
        title="The GenomicMedLab Cool Seq Tool",
        version=__version__,
        description="Common Operations On Lots-of Sequences Tool.",
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
