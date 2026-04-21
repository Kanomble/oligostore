from dataclasses import dataclass
from typing import List

from core.services.sequence_loader import load_sequences_from_sequence_file
from core.services.sequence_utils import reverse_complement

@dataclass
class PrimerBindingHit:
    record_id: str
    start: int
    end: int
    strand: str
    mismatches: int


@dataclass
class PCRProductCandidate:
    record_id: str
    record_name: str
    is_circular_record: bool
    wraps_origin: bool
    forward_start: int
    forward_end: int
    reverse_start: int
    reverse_end: int
    product_start: int
    product_end: int
    product_length: int
    product_sequence: str
    forward_overhang_sequence: str
    reverse_overhang_sequence: str
    forward_mismatches: int
    reverse_mismatches: int


def _record_is_circular(record) -> bool:
    topology = str(getattr(record, "annotations", {}).get("topology", "") or "").strip().lower()
    return topology == "circular"


def _extract_product_sequence(sequence: str, start: int, end: int, *, is_circular: bool) -> str:
    if not sequence:
        return ""
    if start <= end:
        return sequence[start:end]
    if not is_circular:
        return ""
    return sequence[start:] + sequence[:end]


def _build_pcr_product_sequence(
    template_product_sequence: str,
    *,
    forward_overhang_sequence: str,
    reverse_overhang_sequence: str,
) -> str:
    forward_overhang = (forward_overhang_sequence or "").upper().strip()
    reverse_overhang = (reverse_overhang_sequence or "").upper().strip()
    reverse_overhang_on_product = reverse_complement(reverse_overhang) if reverse_overhang else ""
    return f"{forward_overhang}{template_product_sequence}{reverse_overhang_on_product}"

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
    Analyze primer binding against a SequenceFile.
    """
    primer = primer_sequence.upper()
    primer_rc = reverse_complement(primer)

    results: List[PrimerBindingHit] = []

    for record in load_sequences_from_sequence_file(sequence_file):
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


def analyze_primerpair_products(
    *,
    forward_primer_sequence: str,
    reverse_primer_sequence: str,
    forward_overhang_sequence: str = "",
    reverse_overhang_sequence: str = "",
    sequence_file,
    max_mismatches: int = 0,
    block_3prime_mismatch: bool = True,
) -> List[PCRProductCandidate]:
    """
    Find candidate PCR products for a forward/reverse primer pair.

    Coordinates returned here are 1-based and inclusive.
    """
    forward_primer = (forward_primer_sequence or "").upper().strip()
    reverse_primer = (reverse_primer_sequence or "").upper().strip()
    reverse_binding_sequence = reverse_complement(reverse_primer)

    if not forward_primer or not reverse_primer:
        return []

    products: List[PCRProductCandidate] = []

    for record in load_sequences_from_sequence_file(sequence_file):
        seq = str(record.seq).upper()
        record_name = getattr(record, "name", "") or str(record.id)
        is_circular_record = _record_is_circular(record)

        forward_hits = scan_sequence(
            seq,
            forward_primer,
            strand="+",
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )
        reverse_hits = scan_sequence(
            seq,
            reverse_binding_sequence,
            strand="-",
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )

        for forward_hit in forward_hits:
            for reverse_hit in reverse_hits:
                wraps_origin = forward_hit.start > reverse_hit.start
                if wraps_origin and not is_circular_record:
                    continue

                template_product_sequence = _extract_product_sequence(
                    seq,
                    forward_hit.start,
                    reverse_hit.end,
                    is_circular=is_circular_record,
                )
                if not template_product_sequence:
                    continue
                product_sequence = _build_pcr_product_sequence(
                    template_product_sequence,
                    forward_overhang_sequence=forward_overhang_sequence,
                    reverse_overhang_sequence=reverse_overhang_sequence,
                )

                products.append(
                    PCRProductCandidate(
                        record_id=str(record.id),
                        record_name=str(record_name),
                        is_circular_record=is_circular_record,
                        wraps_origin=wraps_origin,
                        forward_start=forward_hit.start + 1,
                        forward_end=forward_hit.end,
                        reverse_start=reverse_hit.start + 1,
                        reverse_end=reverse_hit.end,
                        product_start=forward_hit.start + 1,
                        product_end=reverse_hit.end,
                        product_length=len(product_sequence),
                        product_sequence=product_sequence,
                        forward_overhang_sequence=(forward_overhang_sequence or "").upper().strip(),
                        reverse_overhang_sequence=(reverse_overhang_sequence or "").upper().strip(),
                        forward_mismatches=forward_hit.mismatches,
                        reverse_mismatches=reverse_hit.mismatches,
                    )
                )

    return sorted(
        products,
        key=lambda item: (
            item.record_id,
            item.product_start,
            item.product_end,
            item.forward_mismatches + item.reverse_mismatches,
        ),
    )
