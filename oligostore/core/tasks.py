from celery import shared_task
from django.core.exceptions import ObjectDoesNotExist

from core.models import Primer, PrimerPair, SequenceFile
from core.services.async_jobs import mark_job_failure, mark_job_success
from core.services.primer_binding import (
    analyze_primer_binding,
    analyze_primerpair_products,
)


def _serialize_binding_hits(hits):
    return [
        {
            "record_id": hit.record_id,
            "start": hit.start,
            "end": hit.end,
            "strand": hit.strand,
            "mismatches": hit.mismatches,
        }
        for hit in hits
    ]


def _serialize_pcr_products(products):
    return [
        {
            "record_id": product.record_id,
            "record_name": product.record_name,
            "forward_start": product.forward_start,
            "forward_end": product.forward_end,
            "reverse_start": product.reverse_start,
            "reverse_end": product.reverse_end,
            "product_start": product.product_start,
            "product_end": product.product_end,
            "product_length": product.product_length,
            "product_sequence": product.product_sequence,
            "forward_mismatches": product.forward_mismatches,
            "reverse_mismatches": product.reverse_mismatches,
        }
        for product in products
    ]


@shared_task(bind=True)
def analyze_primer_binding_task(
    self,
    analysis_job_id: int,
    primer_id: int,
    sequence_file_id: int,
    max_mismatches: int = 2,
    block_3prime_mismatch: bool = True,
):
    try:
        primer = Primer.objects.get(id=primer_id)
        sequence_file = SequenceFile.objects.get(id=sequence_file_id)
    except ObjectDoesNotExist as exc:
        message = "Primer or sequence file not found."
        mark_job_failure(analysis_job_id, message)
        raise ValueError(message) from exc

    try:
        hits = analyze_primer_binding(
            primer_sequence=primer.sequence,
            sequence_file=sequence_file,
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )
        payload = _serialize_binding_hits(hits)
        mark_job_success(analysis_job_id, payload)
        return payload
    except Exception as exc:
        mark_job_failure(analysis_job_id, str(exc))
        raise


@shared_task(bind=True)
def analyze_primerpair_products_task(
    self,
    analysis_job_id: int,
    primer_pair_id: int,
    sequence_file_id: int,
    max_mismatches: int = 0,
    block_3prime_mismatch: bool = True,
):
    try:
        primer_pair = PrimerPair.objects.select_related(
            "forward_primer",
            "reverse_primer",
        ).get(id=primer_pair_id)
        sequence_file = SequenceFile.objects.get(id=sequence_file_id)
    except ObjectDoesNotExist as exc:
        message = "Primer pair or sequence file not found."
        mark_job_failure(analysis_job_id, message)
        raise ValueError(message) from exc

    try:
        products = analyze_primerpair_products(
            forward_primer_sequence=primer_pair.forward_primer.sequence,
            reverse_primer_sequence=primer_pair.reverse_primer.sequence,
            sequence_file=sequence_file,
            max_mismatches=max_mismatches,
            block_3prime_mismatch=block_3prime_mismatch,
        )
        payload = {
            "primer_pair": {
                "id": primer_pair.id,
                "name": primer_pair.name,
                "forward_name": primer_pair.forward_primer.primer_name,
                "reverse_name": primer_pair.reverse_primer.primer_name,
            },
            "sequence_file": {
                "id": sequence_file.id,
                "name": sequence_file.name,
            },
            "products": _serialize_pcr_products(products),
        }
        mark_job_success(analysis_job_id, payload)
        return payload
    except Exception as exc:
        mark_job_failure(analysis_job_id, str(exc))
        raise
