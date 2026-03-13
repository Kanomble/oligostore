import json
import os
from functools import lru_cache

from Bio.Restriction import CommOnly
from Bio.Seq import Seq

from ..models import SequenceFeature


@lru_cache(maxsize=8)
def _load_sequences_cached(file_path, file_type, mtime, loader):
    return tuple(loader(file_path, file_type))


def get_sequence_records(sequence_file, loader):
    file_path = sequence_file.file.path
    file_mtime = os.path.getmtime(file_path)
    return _load_sequences_cached(file_path, sequence_file.file_type, file_mtime, loader)


def extract_record_features(record):
    features = []
    for feature in getattr(record, "features", []):
        try:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
        except (TypeError, ValueError):
            continue

        if end < start:
            continue

        qualifiers = getattr(feature, "qualifiers", {}) or {}
        label = (
            (qualifiers.get("label") or [None])[0]
            or (qualifiers.get("gene") or [None])[0]
            or (qualifiers.get("product") or [None])[0]
            or feature.type
        )
        description = (
            (qualifiers.get("note") or [None])[0]
            or (qualifiers.get("function") or [None])[0]
            or (qualifiers.get("product") or [None])[0]
            or (qualifiers.get("gene") or [None])[0]
            or (qualifiers.get("label") or [None])[0]
            or feature.type
        )

        features.append(
            {
                "start": start,
                "end": end,
                "type": feature.type,
                "strand": getattr(feature.location, "strand", None),
                "label": str(label),
                "description": str(description),
                "note": str(description),
            }
        )
    return features


def extract_user_features(sequence_file, record_id):
    return [
        {
            "start": int(feature.start),
            "end": int(feature.end),
            "type": str(feature.feature_type),
            "strand": int(feature.strand),
            "label": str(feature.label),
            "description": str(feature.get_feature_type_display()),
            "note": str(feature.get_feature_type_display()),
            "source": "user",
            "feature_id": feature.id,
            "primer_id": feature.primer_id,
        }
        for feature in SequenceFeature.objects.filter(
            sequence_file=sequence_file,
            record_id=str(record_id),
        ).select_related("primer")
    ]


def extract_restriction_sites_in_range(window_sequence, range_start, range_end, record_length):
    restriction_sites = []
    try:
        restriction_results = CommOnly.search(Seq(window_sequence), linear=True)
        for enzyme, cut_positions in restriction_results.items():
            site_length = len(getattr(enzyme, "site", "") or "")
            if site_length <= 0:
                continue
            cut_offset = int(getattr(enzyme, "fst5", 0))
            for cut_position in cut_positions:
                global_cut_position = range_start + int(cut_position) - 1
                start = global_cut_position - cut_offset
                end = start + site_length - 1
                if start < 1 or end > record_length:
                    continue
                if end < range_start or start > range_end:
                    continue
                restriction_sites.append(
                    {
                        "enzyme": str(enzyme),
                        "site": str(getattr(enzyme, "site", "")),
                        "cut_offset": cut_offset,
                        "start": start,
                        "end": end,
                    }
                )
    except Exception:
        return []
    return restriction_sites


def serialize_record_summary(record):
    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "length": len(record.seq),
    }


def serialize_record_region(record, range_start, range_end, sequence_file):
    record_length = len(record.seq)
    range_start = max(1, min(range_start, record_length))
    range_end = max(1, min(range_end, record_length))
    if range_end < range_start:
        range_start, range_end = range_end, range_start
    window_sequence = str(record.seq[range_start - 1:range_end]).upper()

    combined_features = extract_record_features(record) + extract_user_features(
        sequence_file=sequence_file,
        record_id=record.id,
    )
    combined_features.sort(
        key=lambda f: (
            int(f.get("start", 1)),
            int(f.get("end", 1)),
            str(f.get("label", "")),
        )
    )

    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "length": record_length,
        "region_start": range_start,
        "region_end": range_end,
        "sequence": window_sequence,
        "features": combined_features,
        "restriction_sites": extract_restriction_sites_in_range(
            window_sequence,
            range_start,
            range_end,
            record_length,
        ),
    }


def parse_json_or_form_payload(request):
    if request.body:
        try:
            return json.loads(request.body.decode("utf-8"))
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("Invalid JSON payload.") from exc
    return request.POST
