from dataclasses import dataclass
from typing import List
from core.services.sequence_loader import load_sequences
from core.services.sequence_utils import reverse_complement

@dataclass
class PrimerBindingHit:
    record_id: str
    start: int
    end: int
    strand: str
    mismatches: int

def iter_mismatch_counts(sequence: str, primer: str):
    primer_len = len(primer)
    window_count = len(sequence) - primer_len + 1
    if window_count <= 0:
        return

    mismatch_counts = [0] * window_count
    for offset, primer_base in enumerate(primer):
        for index, seq_base in enumerate(
            sequence[offset : offset + window_count]
        ):
            if seq_base != primer_base:
                mismatch_counts[index] += 1

    yield from mismatch_counts

def scan_sequence(
    sequence: str,
    primer: str,
    strand: str,
    max_mismatches: int,
    block_3prime_mismatch: bool = True,
) -> List[PrimerBindingHit]:
    hits = []
    primer_len = len(primer)
    mismatch_counts = list(iter_mismatch_counts(sequence, primer))

    for i, mismatches in enumerate(mismatch_counts):
        window_end = i + primer_len

        # Enforce perfect 3' base
        if block_3prime_mismatch and sequence[window_end - 1] != primer[-1]:
            continue

        if mismatches <= max_mismatches:
            hits.append(
                PrimerBindingHit(
                    record_id="",
                    start=i,
                    end=window_end,
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
