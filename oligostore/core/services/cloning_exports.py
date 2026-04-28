from copy import deepcopy
from io import StringIO
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ..models import CloningConstruct, PCRProduct, SequenceFeature
from .cloning import _resolve_construct_asset, _validate_restriction_ligation
from .sequence_loader import load_sequences
from .sequence_records import get_sequence_records


def _build_genbank_locus_name(value, fallback="construct"):
    normalized = re.sub(r"\s+", "_", str(value or "").strip())
    normalized = re.sub(r"[^A-Za-z0-9_.-]", "_", normalized)
    normalized = normalized.strip("._-") or fallback
    return normalized[:16]


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


def _clone_feature(feature, mapped_location):
    return SeqFeature(
        location=mapped_location,
        type=str(getattr(feature, "type", "misc_feature")),
        id=getattr(feature, "id", "<unknown id>"),
        qualifiers=deepcopy(getattr(feature, "qualifiers", {}) or {}),
    )


def _iter_user_features(sequence_file, record_id):
    for feature in (
        SequenceFeature.objects.filter(sequence_file=sequence_file, record_id=str(record_id))
        .select_related("primer")
        .order_by("start", "end", "label")
    ):
        qualifiers = {
            "label": [str(feature.label)],
            "note": [str(feature.get_feature_type_display())],
            "source": ["user"],
        }
        if feature.primer_id:
            qualifiers["primer_id"] = [str(feature.primer_id)]
        yield SeqFeature(
            location=FeatureLocation(int(feature.start) - 1, int(feature.end), strand=int(feature.strand)),
            type=str(feature.feature_type),
            qualifiers=qualifiers,
        )


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
    return "".join(
        record_sequence[start:end]
        for start, end, _ in _get_product_segments(product_start, product_end, len(record_sequence), wraps_origin)
    )


def _resolve_sequence_record(sequence_file, record_id):
    records = get_sequence_records(sequence_file, load_sequences)
    record = next((record for record in records if str(record.id) == str(record_id)), None)
    if record is None:
        raise ValueError(f"Record '{record_id}' is not available in sequence file '{sequence_file.name}'.")
    return record


def _build_sequence_file_asset_bundle(asset):
    source_record = _resolve_sequence_record(asset.sequence_file, asset.record_id)
    source_features = list(getattr(source_record, "features", [])) + list(
        _iter_user_features(asset.sequence_file, asset.record_id)
    )
    return {
        "sequence": asset.sequence,
        "features": source_features,
        "annotations": dict(getattr(source_record, "annotations", {}) or {}),
    }


def _build_pcr_product_asset_bundle(asset):
    pcr_product = asset.pcr_product
    source_features = []
    annotations = {}
    if pcr_product is None:
        return {
            "sequence": asset.sequence,
            "features": source_features,
            "annotations": annotations,
        }

    if pcr_product.forward_primer_label:
        source_features.append(
            SeqFeature(
                location=FeatureLocation(0, min(len(asset.sequence), 1)),
                type="primer_bind",
                qualifiers={"label": [str(pcr_product.forward_primer_label)]},
            )
        )
    if pcr_product.reverse_primer_label:
        reverse_start = max(0, len(asset.sequence) - 1)
        source_features.append(
            SeqFeature(
                location=FeatureLocation(reverse_start, len(asset.sequence), strand=-1),
                type="primer_bind",
                qualifiers={"label": [str(pcr_product.reverse_primer_label)]},
            )
        )

    sequence_file = getattr(pcr_product, "sequence_file", None)
    if sequence_file is None or not pcr_product.record_id:
        return {
            "sequence": asset.sequence,
            "features": source_features,
            "annotations": annotations,
        }

    source_record = _resolve_sequence_record(sequence_file, pcr_product.record_id)
    annotations = dict(getattr(source_record, "annotations", {}) or {})
    wraps_origin = int(pcr_product.start) > int(pcr_product.end)
    template_sequence = _extract_product_sequence(
        str(source_record.seq).upper(),
        pcr_product.start,
        pcr_product.end,
        wraps_origin,
    )
    sequence = str(asset.sequence or "").upper()
    template_offset = sequence.find(template_sequence) if template_sequence else -1
    if template_offset == -1:
        return {
            "sequence": sequence,
            "features": source_features,
            "annotations": annotations,
        }

    segments = _get_product_segments(
        product_start=pcr_product.start,
        product_end=pcr_product.end,
        record_length=len(source_record.seq),
        wraps_origin=wraps_origin,
    )
    source_record_features = list(getattr(source_record, "features", [])) + list(
        _iter_user_features(sequence_file, pcr_product.record_id)
    )
    for source_feature in source_record_features:
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
        source_features.append(_clone_feature(source_feature, shifted_location))

    return {
        "sequence": sequence,
        "features": source_features,
        "annotations": annotations,
    }


def _build_asset_bundle(asset):
    if asset.source_type == CloningConstruct.SOURCE_SEQUENCE_FILE:
        return _build_sequence_file_asset_bundle(asset)
    if asset.source_type == CloningConstruct.SOURCE_PCR_PRODUCT:
        return _build_pcr_product_asset_bundle(asset)
    raise ValueError("Unknown cloning asset type.")


def _append_features(record, source_features, segments):
    for source_feature in source_features:
        mapped_location = _build_feature_location(source_feature, segments)
        if mapped_location is None:
            continue
        record.features.append(_clone_feature(source_feature, mapped_location))


def build_cloning_construct_record(construct):
    vector_asset = _resolve_construct_asset(
        source_type=construct.vector_source_type,
        sequence_file=construct.vector_sequence_file,
        pcr_product=construct.vector_pcr_product,
        record_id=construct.vector_record_id,
        label="Vector",
    )
    insert_asset = _resolve_construct_asset(
        source_type=construct.insert_source_type,
        sequence_file=construct.insert_sequence_file,
        pcr_product=construct.insert_pcr_product,
        record_id=construct.insert_record_id,
        label="Insert",
    )
    validation_result = _validate_restriction_ligation(
        vector_asset.sequence,
        insert_asset.sequence,
        construct.left_enzyme,
        construct.right_enzyme,
    )
    if not validation_result.is_valid or validation_result.assembled_sequence != construct.assembled_sequence:
        raise ValueError("Construct cannot be exported because the assembled sequence could not be reproduced.")

    vector_bundle = _build_asset_bundle(vector_asset)
    insert_bundle = _build_asset_bundle(insert_asset)

    record = SeqRecord(
        Seq(construct.assembled_sequence),
        id=_build_genbank_locus_name(construct.name),
        name=_build_genbank_locus_name(construct.name),
        description=str(construct.description or construct.name),
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    source_annotations = vector_bundle.get("annotations") or insert_bundle.get("annotations") or {}
    record.annotations["source"] = str(source_annotations.get("source", "") or "")
    record.annotations["organism"] = str(source_annotations.get("organism", ".") or ".")
    record.annotations["comment"] = (
        f"Restriction ligation construct assembled from {construct.vector_asset_label} and "
        f"{construct.insert_asset_label} using {construct.left_enzyme}/{construct.right_enzyme}."
    )

    vector_features = vector_bundle["features"]
    insert_features = insert_bundle["features"]

    if construct.left_enzyme == construct.right_enzyme:
        from .cloning import _find_cut_positions, _get_enzyme_by_name

        enzyme = _get_enzyme_by_name(construct.left_enzyme)
        if enzyme is None:
            raise ValueError("Construct cannot be exported because the selected enzyme is not available.")
        source_vector_cuts = sorted(_find_cut_positions(vector_asset.sequence, enzyme))
        source_insert_cuts = sorted(_find_cut_positions(insert_asset.sequence, enzyme))

        if len(source_vector_cuts) == 1 and len(source_insert_cuts) == 0:
            left_cut = source_vector_cuts[0]
            _append_features(
                record,
                vector_features,
                [
                    (0, left_cut, 0),
                    (left_cut, len(vector_asset.sequence), left_cut + len(insert_asset.sequence)),
                ],
            )
            _append_features(
                record,
                insert_features,
                [(0, len(insert_asset.sequence), left_cut)],
            )
        elif len(source_vector_cuts) == 2 and len(source_insert_cuts) == 2:
            vector_left_cut, vector_right_cut = source_vector_cuts
            insert_left_cut, insert_right_cut = source_insert_cuts
            insert_fragment_length = insert_right_cut - insert_left_cut
            _append_features(
                record,
                vector_features,
                [
                    (0, vector_left_cut, 0),
                    (vector_right_cut, len(vector_asset.sequence), vector_left_cut + insert_fragment_length),
                ],
            )
            _append_features(
                record,
                insert_features,
                [(insert_left_cut, insert_right_cut, vector_left_cut)],
            )
        else:
            raise ValueError("Construct cannot be exported because the same-enzyme topology is not supported.")
    else:
        from .cloning import _find_cut_positions, _get_enzyme_by_name

        left_enzyme = _get_enzyme_by_name(construct.left_enzyme)
        right_enzyme = _get_enzyme_by_name(construct.right_enzyme)
        if left_enzyme is None or right_enzyme is None:
            raise ValueError("Construct cannot be exported because the selected enzymes are not available.")
        left_hits = sorted(_find_cut_positions(vector_asset.sequence, left_enzyme))
        right_hits = sorted(_find_cut_positions(vector_asset.sequence, right_enzyme))
        if len(left_hits) != 1 or len(right_hits) != 1:
            raise ValueError("Construct cannot be exported because vector cut sites could not be reproduced.")
        left_cut = left_hits[0]
        right_cut = right_hits[0]
        _append_features(
            record,
            vector_features,
            [
                (0, left_cut, 0),
                (right_cut, len(vector_asset.sequence), left_cut + len(insert_asset.sequence)),
            ],
        )
        _append_features(
            record,
            insert_features,
            [(0, len(insert_asset.sequence), left_cut)],
        )

    record.features.append(
        SeqFeature(
            location=FeatureLocation(0, len(record.seq)),
            type="source",
            qualifiers={"label": [construct.name]},
        )
    )
    return record


def export_cloning_construct_genbank(construct):
    handle = StringIO()
    SeqIO.write(build_cloning_construct_record(construct), handle, "genbank")
    return handle.getvalue()
