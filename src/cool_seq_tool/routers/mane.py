"""Module containing routes related to MANE data"""

import logging

from fastapi import APIRouter, Query

from cool_seq_tool.routers import (
    RESP_DESCR,
    SERVICE_NAME,
    UNHANDLED_EXCEPTION_MSG,
    Tags,
    cool_seq_tool,
)
from cool_seq_tool.schemas import AnnotationLayer, ManeDataService, ResidueMode
from cool_seq_tool.utils import service_meta

logger = logging.getLogger("cool_seq_tool")

router = APIRouter(prefix=f"/{SERVICE_NAME}/mane")


ref_descr = (
    "Reference at position given during input. When this is set, it will "
    "ensure that the reference sequences match for the final result."
)
try_longest_compatible_descr = (
    "`True` if should try longest compatible remaining if"
    " mane transcript was not compatible. `False` otherwise."
)


@router.get(
    "/get_mane_data",
    summary="Retrieve MANE data in inter-residue coordinates",
    response_description=RESP_DESCR,
    description="Return MANE Select, MANE Plus Clinical, or Longest Remaining "
    "Transcript data in inter-residue coordinates. See our docs for "
    "more information on transcript priority.",
    response_model=ManeDataService,
    tags=[Tags.MANE_TRANSCRIPT],
)
async def get_mane_data(
    ac: str = Query(..., description="Accession"),
    start_pos: int = Query(..., description="Start position"),
    start_annotation_layer: AnnotationLayer = Query(
        ..., description="Starting annotation layer for query"
    ),
    end_pos: int | None = Query(
        None, description="End position. If not set, will set to `start_pos`."
    ),
    gene: str | None = Query(None, description="HGNC gene symbol"),
    ref: str | None = Query(None, description=ref_descr),
    try_longest_compatible: bool = Query(
        True, description=try_longest_compatible_descr
    ),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE, description="Residue mode for position(s)"
    ),
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
    warnings = []
    mane_data = None
    try:
        mane_data = await cool_seq_tool.mane_transcript.get_mane_transcript(
            ac=ac,
            start_pos=start_pos,
            end_pos=end_pos,
            start_annotation_layer=start_annotation_layer,
            gene=gene,
            ref=ref,
            try_longest_compatible=try_longest_compatible,
            residue_mode=residue_mode,
        )

        if not mane_data:
            warnings.append("Unable to retrieve MANE data")
    except Exception as e:
        logger.exception("get_mane_data unhandled exception %s", e)
        warnings.append(UNHANDLED_EXCEPTION_MSG)

    return ManeDataService(
        mane_data=mane_data, warnings=warnings, service_meta=service_meta()
    )
