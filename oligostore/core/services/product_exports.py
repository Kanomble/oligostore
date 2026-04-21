from copy import deepcopy
from io import StringIO

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ..models import SequenceFeature
from .sequence_loader import load_sequences
from .sequence_records import get_sequence_records
from .sequence_utils import reverse_complement


def _get_product_segments(product_start, product_end, record_length, wraps_origin):
    start_index = int(product_start) - 1
    end_index = int(product_end)
    if wraps_origin:
        return [
            (start_index, record_length, 0),
            (0, end_index, record_length - start_index),
        ]
    return [(start_index, end_index, 0)]


def _extract_product_sequence(record_sequence, product_start, product_end, wraps_origin):
    segments = _get_product_segments(
        product_start=product_start,
        product_end=product_end,
        record_length=len(record_sequence),
        wraps_origin=wraps_origin,
    )
    return "".join(record_sequence[start:end] for start, end, _ in segments)


def _normalize_overhang(sequence):
    return str(sequence or "").strip().upper()


def _map_location_parts(parts, segments):
    mapped_parts = []
    for part in parts:
        part_start = int(part.start)
        part_end = int(part.end)
        part_strand = getattr(part, "strand", None)
        for segment_start, segment_end, segment_offset in segments:
            overlap_start = max(part_start, segment_start)
            overlap_end = min(part_end, segment_end)
            if overlap_start >= overlap_end:
                continue
            mapped_parts.append(
                FeatureLocation(
                    segment_offset + (overlap_start - segment_start),
                    segment_offset + (overlap_end - segment_start),
                    strand=part_strand,
                )
            )
    return mapped_parts


def _build_feature_location(feature, segments):
    feature_location = getattr(feature, "location", None)
    if feature_location is None:
        return None
    parts = list(feature_location.parts) if isinstance(feature_location, CompoundLocation) else [feature_location]
    mapped_parts = _map_location_parts(parts, segments)
    if not mapped_parts:
        return None
    if len(mapped_parts) == 1:
        return mapped_parts[0]
    return CompoundLocation(mapped_parts)


def _iter_user_features(sequence_file, record_id):
    for feature in (
        SequenceFeature.objects.filter(sequence_file=sequence_file, record_id=str(record_id))
        .select_related("primer")
        .order_by("start", "end", "label")
    ):
        qualifiers = {
            "label": [str(feature.label)],
            "note": [str(feature.get_feature_type_display())],
        }
        if feature.primer_id:
            qualifiers["primer_id"] = [str(feature.primer_id)]
        yield SeqFeature(
            location=FeatureLocation(int(feature.start) - 1, int(feature.end), strand=int(feature.strand)),
            type=str(feature.feature_type),
            qualifiers=qualifiers,
        )


def build_product_record(
    *,
    sequence_file,
    record_id,
    product_start,
    product_end,
    wraps_origin=False,
    forward_overhang_sequence="",
    reverse_overhang_sequence="",
    exported_name=None,
):
    records = get_sequence_records(sequence_file, load_sequences)
    source_record = next((record for record in records if str(record.id) == str(record_id)), None)
    if source_record is None:
        raise ValueError("Selected source record was not found.")

    record_sequence = str(source_record.seq).upper()
    template_product_sequence = _extract_product_sequence(
        record_sequence,
        product_start=product_start,
        product_end=product_end,
        wraps_origin=wraps_origin,
    )
    if not template_product_sequence:
        raise ValueError("Could not build PCR product sequence.")
    forward_overhang = _normalize_overhang(forward_overhang_sequence)
    reverse_overhang = _normalize_overhang(reverse_overhang_sequence)
    product_sequence = (
        f"{forward_overhang}"
        f"{template_product_sequence}"
        f"{reverse_complement(reverse_overhang) if reverse_overhang else ''}"
    )

    product_name = exported_name or f"{source_record.id}_{product_start}_{product_end}"
    product_record = SeqRecord(
        Seq(product_sequence),
        id=str(product_name),
        name=str(product_name),
        description=f"PCR product from {source_record.id}",
    )
    product_record.annotations["molecule_type"] = "DNA"
    product_record.annotations["topology"] = "linear"
    product_record.annotations["source"] = str(source_record.annotations.get("source", "") or "")
    product_record.annotations["organism"] = str(source_record.annotations.get("organism", ".") or ".")
    if forward_overhang or reverse_overhang:
        product_record.annotations["comment"] = (
            f"Forward primer overhang included: {forward_overhang}; "
            f"reverse primer overhang included on product strand as "
            f"{reverse_complement(reverse_overhang) if reverse_overhang else ''}."
        )
    segments = _get_product_segments(
        product_start=product_start,
        product_end=product_end,
        record_length=len(record_sequence),
        wraps_origin=wraps_origin,
    )
    template_offset = len(forward_overhang)

    source_features = list(getattr(source_record, "features", [])) + list(
        _iter_user_features(sequence_file, record_id)
    )
    for source_feature in source_features:
        mapped_location = _build_feature_location(source_feature, segments)
        if mapped_location is None:
            continue
        if isinstance(mapped_location, CompoundLocation):
            shifted_location = CompoundLocation(
                [
                    FeatureLocation(
                        int(part.start) + template_offset,
                        int(part.end) + template_offset,
                        strand=getattr(part, "strand", None),
                    )
                    for part in mapped_location.parts
                ]
            )
        else:
            shifted_location = FeatureLocation(
                int(mapped_location.start) + template_offset,
                int(mapped_location.end) + template_offset,
                strand=getattr(mapped_location, "strand", None),
            )
        product_record.features.append(
            SeqFeature(
                location=shifted_location,
                type=str(getattr(source_feature, "type", "misc_feature")),
                id=getattr(source_feature, "id", "<unknown id>"),
                qualifiers=deepcopy(getattr(source_feature, "qualifiers", {}) or {}),
            )
        )

    return product_record


def export_product_fasta(product_record):
    return f">{product_record.id}\n{str(product_record.seq)}\n"


def export_product_genbank(product_record):
    handle = StringIO()
    SeqIO.write(product_record, handle, "genbank")
    return handle.getvalue()
