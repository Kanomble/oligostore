from dataclasses import dataclass
from typing import List

from Bio.Seq import Seq

from core.services.sequence_loader import load_sequences


@dataclass
class PrimerBindingHit:
    record_id: str
    start: int
    end: int
    strand: str
    mismatches: int


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def count_mismatches(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y)


def scan_sequence(
    sequence: str,
    primer: str,
    strand: str,
    max_mismatches: int,
    block_3prime_mismatch: bool = True,
) -> List[PrimerBindingHit]:
    hits = []
    primer_len = len(primer)

    for i in range(len(sequence) - primer_len + 1):
        window = sequence[i : i + primer_len]

        # Enforce perfect 3' base
        if block_3prime_mismatch and window[-1] != primer[-1]:
            continue

        mismatches = count_mismatches(window, primer)
        if mismatches <= max_mismatches:
            hits.append(
                PrimerBindingHit(
                    record_id="",
                    start=i,
                    end=i + primer_len,
                    strand=strand,
                    mismatches=mismatches,
                )
            )

    return hits


def analyze_primer_binding(
    primer_sequence: str,
    sequence_file,
    max_mismatches: int = 2,
    block_3prime_mismatch: bool = True,
) -> List[PrimerBindingHit]:
    """
    Analyze primer binding against a SequenceFile (FASTA or GenBank).
    """
    primer = primer_sequence.upper()
    primer_rc = reverse_complement(primer)

    results: List[PrimerBindingHit] = []

    for record in load_sequences(
        sequence_file.file.path,
        sequence_file.file_type,
    ):
        seq = str(record.seq).upper()

        fwd_hits = scan_sequence(
            seq,
            primer,
            strand="+",
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )

        rev_hits = scan_sequence(
            seq,
            primer_rc,
            strand="-",
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )

        for hit in fwd_hits + rev_hits:
            hit.record_id = record.id
            results.append(hit)

    return results
