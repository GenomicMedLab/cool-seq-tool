"""Module containing validation functions"""
from typing import Dict, Tuple, Optional

from cool_seq_tool.schemas import AnnotationLayer, ResidueMode
from cool_seq_tool.data_sources import SeqRepoAccess


class ValidationError(Exception):
    """Custom exception for validation errors"""


# TODO: Figure out what residue_mode is being used
def validate_reading_frames(
    ac: str, start_pos: int, end_pos: int, transcript_data: Dict
) -> None:
    """Return whether reading frames are the same after translation.

    :param str ac: Query accession
    :param int start_pos: Original start cDNA position change
    :param int end_pos: Original end cDNA position change
    :param Dict transcript_data: Ensembl and RefSeq transcripts with
        corresponding position change
    :raises ValidationError: If reading frame validation checks do not pass
    """
    for pos, pos_index in [(start_pos, 0), (end_pos, 1)]:
        if pos is not None:
            og_rf = pos % 3 or 3
            new_rf = pos % 3 or 3

            if og_rf != new_rf:
                raise ValidationError(
                    f"{ac} original reading frame ({og_rf}) does not match new "
                    f"{transcript_data['ensembl']}, {transcript_data['refseq']} "
                    f"reading frame ({new_rf})"
                )
        else:
            if pos_index == 0:
                raise ValidationError(f"{ac} must have start position")

def validate_references(
    seqrepo_access: SeqRepoAccess, ac: str, coding_start_site: int, start_pos: int,
    end_pos: int, mane_transcript: Dict, expected_ref: str, anno: AnnotationLayer,
    residue_mode: str
) -> None:
    """Return whether or not reference changes are the same.

    :param SeqRepoAccess: SeqRepoAccess client
    :param str ac: Query accession
    :param int coding_start_site: ac's coding start site
    :param int start_pos: Original start position change
    :param int end_pos: Origin end position change
    :param Dict mane_transcript: Ensembl and RefSeq transcripts with
        corresponding position change
    :param str expected_ref: Reference at position given during input
    :param AnnotationLayer anno: Annotation layer we are starting from
    :param ResidueMode residue_mode: Residue mode
    :raise ValidationError: If reference sequence checks do not pass
    """
    if anno == AnnotationLayer.CDNA:
        start_pos += coding_start_site
        end_pos += coding_start_site

    ref, _ = seqrepo_access.get_reference_sequence(
        ac, start_pos, end=end_pos, residue_mode=residue_mode
    )
    if ref is None:
        return False

    if mane_transcript:
        mane_start_pos = mane_transcript["pos"][0]
        mane_end_pos = mane_transcript["pos"][1]
        if anno == AnnotationLayer.CDNA:
            mane_cds = mane_transcript["coding_start_site"]
            mane_start_pos += mane_cds
            mane_end_pos += mane_cds
        mane_ref, _ = seqrepo_access.get_reference_sequence(
            mane_transcript["refseq"],
            mane_start_pos,
            end=mane_end_pos if mane_start_pos != mane_end_pos else None,
            residue_mode=residue_mode
        )
        if not mane_ref:
            raise ValidationError("Unable to validate reference for MANE Transcript")

        if expected_ref != mane_ref:
            raise ValidationError(f"Expected ref, {expected_ref}, but got {mane_ref}"
                                  f" on MANE accession, {mane_transcript['refseq']}")

    if expected_ref != ref:
        raise ValidationError(f"Expected ref, {expected_ref}, but got {ref} on "
                              f"accession, {ac}")


def validate_index(
    seqrepo_access: SeqRepoAccess, ac: str, pos: Tuple[int, int],
    coding_start_site: int
) -> None:
    """Validate that positions actually exist on accession

    :param SeqRepoAccess: SeqRepoAccess client
    :param str ac: Accession
    :param Tuple[int, int] pos: Start position change, End position change
    :param int coding_start_site: coding start site for accession
    :return: `True` if positions exist on accession. `False` otherwise
    :raise ValidationError: If index checks do not pass
    """
    start_pos = pos[0] + coding_start_site
    end_pos = pos[1] + coding_start_site
    if not seqrepo_access.get_reference_sequence(
        ac, start_pos, end_pos, residue_mode=ResidueMode.INTER_RESIDUE)[0]:
        raise ValidationError("Invalid index")
