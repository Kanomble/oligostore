from celery import shared_task
from django.core.exceptions import ObjectDoesNotExist

from core.models import Primer, SequenceFile
from core.services.primer_binding import analyze_primer_binding


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


@shared_task(bind=True)
def analyze_primer_binding_task(
    self,
    primer_id: int,
    sequence_file_id: int,
    max_mismatches: int = 2,
    block_3prime_mismatch: bool = True,
):
    try:
        primer = Primer.objects.get(id=primer_id)
        sequence_file = SequenceFile.objects.get(id=sequence_file_id)
    except ObjectDoesNotExist as exc:
        raise ValueError("Primer or sequence file not found.") from exc

    hits = analyze_primer_binding(
        primer_sequence=primer.sequence,
        sequence_file=sequence_file,
        max_mismatches=max_mismatches,
        block_3prime_mismatch=block_3prime_mismatch,
    )
    return _serialize_binding_hits(hits)