"""Shared nucleotide composition, CDS validation, and translation rules."""

from dataclasses import dataclass

VALID_DNA_BASES = frozenset("ACGT")
IUPAC_DNA_BASES = frozenset("ACGTRYSWKMBDHVN")
BACTERIAL_START_CODONS = frozenset({"ATG", "GTG", "TTG"})
STOP_CODONS = frozenset({"TAA", "TAG", "TGA"})
CODON_TO_AMINO_ACID = {
    "AAA": "K", "AAC": "N", "AAG": "K", "AAT": "N",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
    "AGA": "R", "AGC": "S", "AGG": "R", "AGT": "S",
    "ATA": "I", "ATC": "I", "ATG": "M", "ATT": "I",
    "CAA": "Q", "CAC": "H", "CAG": "Q", "CAT": "H",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
    "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
    "GAA": "E", "GAC": "D", "GAG": "E", "GAT": "D",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
    "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
    "TAA": "*", "TAC": "Y", "TAG": "*", "TAT": "Y",
    "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
    "TGA": "*", "TGC": "C", "TGG": "W", "TGT": "C",
    "TTA": "L", "TTC": "F", "TTG": "L", "TTT": "F",
}

__all__ = [
    "BACTERIAL_START_CODONS",
    "CODON_TO_AMINO_ACID",
    "IUPAC_DNA_BASES",
    "STOP_CODONS",
    "VALID_DNA_BASES",
    "AlleleSequenceValidation",
    "NucleotideBaseComposition",
    "contains_in_frame_stop_codon",
    "count_nucleotide_base_composition",
    "evaluate_allele_sequence",
    "is_eligible_iupac_dna_sequence",
    "is_valid_allele_sequence",
    "normalize_nucleotide_sequence",
    "translate_valid_allele",
]


@dataclass(frozen=True)
class NucleotideBaseComposition:
    """Canonical and excluded nucleotide-base counts."""

    guanine_cytosine_base_count: int
    canonical_nucleotide_base_count: int
    excluded_nucleotide_base_count: int

    @property
    def guanine_cytosine_content_percentage(self) -> float | None:
        """Return GC percentage, or ``None`` without canonical bases."""
        if self.canonical_nucleotide_base_count == 0:
            return None
        return (
            self.guanine_cytosine_base_count
            / self.canonical_nucleotide_base_count
            * 100
        )


@dataclass(frozen=True)
class AlleleSequenceValidation:
    """Independent coding-sequence validation flags."""

    has_incomplete_orf: bool
    has_ambiguous_bases: bool
    is_missing_start_or_stop_codon: bool
    has_in_frame_stop_codon: bool

    @property
    def is_valid(self) -> bool:
        """Return whether every coding-sequence rule passes."""
        return not any(
            (
                self.has_incomplete_orf,
                self.has_ambiguous_bases,
                self.is_missing_start_or_stop_codon,
                self.has_in_frame_stop_codon,
            )
        )


def normalize_nucleotide_sequence(nucleotide_sequence: str) -> str:
    """Remove whitespace and normalize nucleotide letters to uppercase."""
    return "".join(nucleotide_sequence.split()).upper()


def count_nucleotide_base_composition(
    nucleotide_sequence: str,
) -> NucleotideBaseComposition:
    """Count canonical GC/ACGT bases and excluded non-whitespace characters."""
    normalized_nucleotide_sequence = normalize_nucleotide_sequence(
        nucleotide_sequence
    )
    canonical_nucleotide_base_count = sum(
        nucleotide_base in VALID_DNA_BASES
        for nucleotide_base in normalized_nucleotide_sequence
    )
    return NucleotideBaseComposition(
        guanine_cytosine_base_count=sum(
            nucleotide_base in {"G", "C"}
            for nucleotide_base in normalized_nucleotide_sequence
        ),
        canonical_nucleotide_base_count=canonical_nucleotide_base_count,
        excluded_nucleotide_base_count=(
            len(normalized_nucleotide_sequence) - canonical_nucleotide_base_count
        ),
    )


def contains_in_frame_stop_codon(normalized_allele_sequence: str) -> bool:
    """Return whether a stop codon occurs before a valid terminal stop."""
    final_complete_codon_end_position = len(normalized_allele_sequence) - (
        len(normalized_allele_sequence) % 3
    )
    in_frame_codon_end_position = final_complete_codon_end_position
    if (
        len(normalized_allele_sequence) % 3 == 0
        and normalized_allele_sequence[-3:] in STOP_CODONS
    ):
        in_frame_codon_end_position -= 3
    return any(
        normalized_allele_sequence[
            codon_start_position : codon_start_position + 3
        ]
        in STOP_CODONS
        for codon_start_position in range(3, in_frame_codon_end_position, 3)
    )


def evaluate_allele_sequence(
    normalized_allele_sequence: str,
) -> AlleleSequenceValidation:
    """Evaluate the four CDS rules for an already normalized sequence."""
    return AlleleSequenceValidation(
        has_incomplete_orf=len(normalized_allele_sequence) % 3 != 0,
        has_ambiguous_bases=not set(normalized_allele_sequence).issubset(
            VALID_DNA_BASES
        ),
        is_missing_start_or_stop_codon=(
            normalized_allele_sequence[:3] not in BACTERIAL_START_CODONS
            or normalized_allele_sequence[-3:] not in STOP_CODONS
        ),
        has_in_frame_stop_codon=contains_in_frame_stop_codon(
            normalized_allele_sequence
        ),
    )


def is_valid_allele_sequence(allele_sequence: str) -> bool:
    """Return whether an allele passes every coding-sequence rule."""
    return evaluate_allele_sequence(
        normalize_nucleotide_sequence(allele_sequence)
    ).is_valid


def is_eligible_iupac_dna_sequence(
    normalized_nucleotide_sequence: str,
) -> bool:
    """Return whether a non-empty normalized sequence contains IUPAC DNA."""
    return bool(normalized_nucleotide_sequence) and set(
        normalized_nucleotide_sequence
    ).issubset(IUPAC_DNA_BASES)


def translate_valid_allele(allele_sequence: str) -> str:
    """Translate a valid CDS, normalizing its bacterial start codon to M.

    The caller must supply a sequence that passes :func:`evaluate_allele_sequence`.
    The terminal stop codon is deliberately omitted from the returned protein.
    """
    normalized_allele_sequence = normalize_nucleotide_sequence(allele_sequence)
    translated_amino_acids = ["M"]
    translated_amino_acids.extend(
        CODON_TO_AMINO_ACID[
            normalized_allele_sequence[
                codon_start_position : codon_start_position + 3
            ]
        ]
        for codon_start_position in range(
            3,
            len(normalized_allele_sequence) - 3,
            3,
        )
    )
    return "".join(translated_amino_acids)
