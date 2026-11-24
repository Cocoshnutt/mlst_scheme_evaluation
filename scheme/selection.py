"""Deterministic modal-length-proximal allele selection."""

from collections import Counter

from .fasta import FastaRecord
from .sequence import evaluate_allele_sequence, is_eligible_iupac_dna_sequence

__all__ = [
    "select_allele_closest_to_modal_length",
    "select_iupac_allele_closest_to_modal_length",
    "select_valid_cds_closest_to_modal_length",
]


def select_allele_closest_to_modal_length(
    locus_records: tuple[FastaRecord, ...],
    *,
    require_valid_coding_sequence: bool,
) -> FastaRecord | None:
    """Select the first eligible record closest to the shortest tied mode.

    The mode uses every non-empty record, including invalid CDS records.
    Original FASTA order resolves equal distances.
    """
    non_empty_allele_lengths = [
        len(record.sequence) for record in locus_records if record.sequence
    ]
    if not non_empty_allele_lengths:
        return None
    length_counts = Counter(non_empty_allele_lengths)
    highest_length_frequency = max(length_counts.values())
    modal_allele_length = min(
        allele_length
        for allele_length, frequency in length_counts.items()
        if frequency == highest_length_frequency
    )
    for _, fasta_record in sorted(
        enumerate(locus_records),
        key=lambda indexed_record: (
            abs(len(indexed_record[1].sequence) - modal_allele_length),
            indexed_record[0],
        ),
    ):
        is_eligible = (
            evaluate_allele_sequence(fasta_record.sequence).is_valid
            if require_valid_coding_sequence
            else is_eligible_iupac_dna_sequence(fasta_record.sequence)
        )
        if is_eligible:
            return fasta_record
    return None


def select_valid_cds_closest_to_modal_length(
    locus_records: tuple[FastaRecord, ...],
) -> FastaRecord | None:
    """Select a valid CDS using the established modal-length rule."""
    return select_allele_closest_to_modal_length(
        locus_records,
        require_valid_coding_sequence=True,
    )


def select_iupac_allele_closest_to_modal_length(
    locus_records: tuple[FastaRecord, ...],
) -> FastaRecord | None:
    """Select an IUPAC DNA record using the established modal-length rule."""
    return select_allele_closest_to_modal_length(
        locus_records,
        require_valid_coding_sequence=False,
    )
