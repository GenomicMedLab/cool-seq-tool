"""Module containing routes related to alignment mapping"""
import logging
from typing import Optional

from fastapi import APIRouter
from fastapi import Query


from cool_seq_tool.routers import cool_seq_tool, SERVICE_NAME, RESP_DESCR, Tags
from cool_seq_tool.schemas import Assembly, ToGenomicService, ToCdnaService, \
    ResidueMode

logger = logging.getLogger("cool_seq_tool")

router = APIRouter(prefix=f"/{SERVICE_NAME}/alignment_mapper")


@router.get(
    "/p_to_c",
    summary="Translate protein representation to cDNA representation",
    response_description=RESP_DESCR,
    description="Given protein accession and positions, return associated cDNA "
                "accession and positions to codon(s)",
    response_model=ToCdnaService,
    tags=[Tags.ALIGNMENT_MAPPER]
)
async def p_to_c(
    p_ac: str = Query(..., description="Protein RefSeq accession"),
    p_start_pos: int = Query(..., description="Protein start position"),
    p_end_pos: int = Query(..., description="Protein end position"),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `p_start_pos` and `p_end_pos`")
) -> ToCdnaService:
    """Translate protein representation to cDNA representation

    :param str p_ac: Protein RefSeq accession
    :param int p_start_pos: Protein start position
    :param int p_end_pos: Protein end position
    :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`.
    :return: ToCdnaService containing cDNA representation, warnings, and
        service meta
    """
    try:
        c_data, w = await cool_seq_tool.alignment_mapper.p_to_c(
            p_ac, p_start_pos, p_end_pos, residue_mode)
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        c_data = None
    return ToCdnaService(
        c_data=c_data,
        warnings=[w] if w else [],
        service_meta=cool_seq_tool.service_meta()
    )


@router.get(
    "/c_to_g",
    summary="Translate cDNA representation to genomic representation",
    response_description=RESP_DESCR,
    description="Given cDNA accession and positions for codon(s), return associated genomic"  # noqa: E501
                " accession and positions for a given target genome assembly",
    response_model=ToGenomicService,
    tags=[Tags.ALIGNMENT_MAPPER]
)
async def c_to_g(
    c_ac: str = Query(..., description="cDNA RefSeq accession"),
    c_start_pos: int = Query(..., description="cDNA start position for codon"),
    c_end_pos: int = Query(..., description="cDNA end position for codon"),
    cds_start: Optional[int] = Query(
        None, description="CDS start site. If not provided, this will be computed."),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `c_start_pos` and `c_end_pos`"),
    target_genome_assembly: Assembly = Query(Assembly.GRCH38,
                                             description="Genomic assembly to map to")
) -> ToGenomicService:
    """Translate cDNA representation to genomic representation

    :param str c_ac: cDNA RefSeq accession
    :param int c_start_pos: cDNA start position for codon
    :param int c_end_pos: cDNA end position for codon
    :param Optional[int] cds_start: CDS start site. If not provided, this will be
        computed.
    :param ResidueMode residue_mode: Residue mode for `c_start_pos` and `c_end_pos`.
    :param Assembly target_genome_assembly: Genome assembly to get genomic data for
    :return: ToGenomicService containing genomic representation, warnings, and
        service meta
    """
    try:
        g_data, w = await cool_seq_tool.alignment_mapper.c_to_g(
            c_ac, c_start_pos, c_end_pos, cds_start=cds_start,
            residue_mode=residue_mode,
            target_genome_assembly=target_genome_assembly)
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        g_data = None
    return ToGenomicService(
        g_data=g_data,
        warnings=[w] if w else [],
        service_meta=cool_seq_tool.service_meta()
    )


@router.get(
    "/p_to_g",
    summary="Translate protein representation to genomic representation",
    response_description=RESP_DESCR,
    description="Given protein accession and positions, return associated genomic "
                "accession and positions for a given target genome assembly",
    response_model=ToGenomicService,
    tags=[Tags.ALIGNMENT_MAPPER]
)
async def p_to_g(
    p_ac: str = Query(..., description="Protein RefSeq accession"),
    p_start_pos: int = Query(..., description="Protein start position"),
    p_end_pos: int = Query(..., description="Protein end position"),
    residue_mode: ResidueMode = Query(
        ResidueMode.RESIDUE,
        description="Residue mode for `p_start_pos` and `p_end_pos`"),
    target_genome_assembly: Assembly = Query(Assembly.GRCH38,
                                             description="Genomic assembly to map to")
) -> ToGenomicService:
    """Translate protein representation to genomic representation

    :param str p_ac: Protein RefSeq accession
    :param int p_start_pos: Protein start position
    :param int p_end_pos: Protein end position
    :param ResidueMode residue_mode: Residue mode for `p_start_pos` and `p_end_pos`.
    :param Assembly target_genome_assembly: Genome assembly to get genomic data for
    :return: ToGenomicService containing genomic representation, warnings, and
        service meta
    """
    try:
        g_data, w = await cool_seq_tool.alignment_mapper.p_to_g(
            p_ac, p_start_pos, p_end_pos, residue_mode=residue_mode,
            target_genome_assembly=target_genome_assembly)
    except Exception as e:
        logger.error("Unhandled exception: %s", str(e))
        w = "Unhandled exception. See logs for more information."
        g_data = None
    return ToGenomicService(
        g_data=g_data,
        warnings=[w] if w else [],
        service_meta=cool_seq_tool.service_meta()
    )
