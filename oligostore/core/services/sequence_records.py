import json
import os
import re
from functools import lru_cache

from Bio.Restriction import CommOnly
from Bio.Seq import Seq

from ..models import SequenceFeature


IUPAC_DNA_BASES = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T",
    "R": "AG",
    "Y": "CT",
    "S": "GC",
    "W": "AT",
    "K": "GT",
    "M": "AC",
    "B": "CGT",
    "D": "AGT",
    "H": "ACT",
    "V": "ACG",
    "N": "ACGT",
}
IUPAC_COMPLEMENTS = str.maketrans("ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN")


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


def reverse_complement_site(site):
    return str(site).upper().translate(IUPAC_COMPLEMENTS)[::-1]


def site_regex(site):
    parts = []
    for base in str(site).upper():
        values = IUPAC_DNA_BASES.get(base)
        if values is None:
            parts.append(re.escape(base))
        elif len(values) == 1:
            parts.append(values)
        else:
            parts.append(f"[{values}]")
    return re.compile(f"(?=({''.join(parts)}))")


def find_recognition_matches(window_sequence, site):
    site = str(site or "").upper()
    reverse_site = reverse_complement_site(site)
    matches = []
    for strand, motif in ((1, site), (-1, reverse_site)):
        if not motif:
            continue
        if strand == -1 and motif == site:
            continue
        pattern = site_regex(motif)
        matches.extend(
            {
                "start": match.start() + 1,
                "end": match.start() + len(motif),
                "strand": strand,
            }
            for match in pattern.finditer(window_sequence)
        )
    return matches


def recognition_cut_boundaries(match, site_length, cut_offset, reverse_cut_offset):
    if match["strand"] == -1:
        return {
            "forward": match["start"] - reverse_cut_offset,
            "reverse": match["start"] - (cut_offset - site_length),
        }
    return {
        "forward": match["start"] + cut_offset,
        "reverse": match["end"] + reverse_cut_offset + 1,
    }


def extract_restriction_sites_in_range(window_sequence, range_start, range_end, record_length):
    restriction_sites = []
    try:
        restriction_results = CommOnly.search(Seq(window_sequence), linear=True)
        for enzyme, cut_positions in restriction_results.items():
            site = str(getattr(enzyme, "site", "") or "")
            site_length = len(site)
            if site_length <= 0:
                continue
            cut_offset = int(getattr(enzyme, "fst5", 0))
            reverse_cut_offset = int(getattr(enzyme, "fst3", cut_offset))
            recognition_matches = find_recognition_matches(window_sequence, site)
            matches_by_forward_cut = {}
            for match in recognition_matches:
                boundaries = recognition_cut_boundaries(
                    match,
                    site_length,
                    cut_offset,
                    reverse_cut_offset,
                )
                matches_by_forward_cut.setdefault(boundaries["forward"], []).append(
                    (match, boundaries)
                )
            for cut_position in cut_positions:
                global_cut_position = range_start + int(cut_position) - 1
                local_cut_position = int(cut_position)
                matched_sites = matches_by_forward_cut.get(local_cut_position)
                if matched_sites:
                    match, boundaries = matched_sites[0]
                    start = range_start + match["start"] - 1
                    end = range_start + match["end"] - 1
                    recognition_strand = match["strand"]
                    forward_cut_boundary = range_start + boundaries["forward"] - 1
                    reverse_cut_boundary = range_start + boundaries["reverse"] - 1
                else:
                    start = global_cut_position - cut_offset
                    end = start + site_length - 1
                    recognition_strand = 1
                    forward_cut_boundary = global_cut_position
                    reverse_cut_boundary = end + reverse_cut_offset + 1
                if start < 1 or end > record_length:
                    continue
                site_start = min(start, end, forward_cut_boundary, reverse_cut_boundary)
                site_end = max(start, end, forward_cut_boundary, reverse_cut_boundary)
                if site_end < range_start or site_start > range_end:
                    continue
                restriction_sites.append(
                    {
                        "enzyme": str(enzyme),
                        "site": site,
                        "cut_offset": cut_offset,
                        "cut_offset_3": reverse_cut_offset,
                        "start": start,
                        "end": end,
                        "recognition_start": start,
                        "recognition_end": end,
                        "recognition_strand": recognition_strand,
                        "cut_boundary_forward": forward_cut_boundary,
                        "cut_boundary_reverse": reverse_cut_boundary,
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
