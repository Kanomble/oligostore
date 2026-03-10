import os
import json
from functools import lru_cache

from celery.result import AsyncResult
from django.contrib.auth.decorators import login_required
from django.core.exceptions import ValidationError
from django.db.models import Case, IntegerField, Q, Value, When
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages
from Bio.Seq import Seq
from Bio.Restriction import CommOnly

from ..models import Primer, SequenceFeature, SequenceFile
from ..forms import clean_optional_sequence_value, clean_sequence_value
from ..services.primer_binding import analyze_primer_binding
from ..services.sequence_loader import load_sequences
from ..tasks import analyze_primer_binding_task
from .utils import paginate_queryset


@lru_cache(maxsize=8)
def _load_sequences_cached(file_path, file_type, mtime):
    # mtime participates in the cache key so file updates invalidate cache.
    return tuple(load_sequences(file_path, file_type))


def _get_sequence_records(sequence_file):
    file_path = sequence_file.file.path
    file_mtime = os.path.getmtime(file_path)
    return _load_sequences_cached(file_path, sequence_file.file_type, file_mtime)


def _extract_record_features(record):
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

        features.append(
            {
                "start": start,
                "end": end,
                "type": feature.type,
                "strand": getattr(feature.location, "strand", None),
                "label": str(label),
            }
        )
    return features


def _extract_user_features(sequence_file, record_id):
    features = []
    for feature in SequenceFeature.objects.filter(
        sequence_file=sequence_file,
        record_id=str(record_id),
    ).select_related("primer"):
        features.append(
            {
                "start": int(feature.start),
                "end": int(feature.end),
                "type": str(feature.feature_type),
                "strand": int(feature.strand),
                "label": str(feature.label),
                "source": "user",
                "primer_id": feature.primer_id,
            }
        )
    return features


def _extract_restriction_sites_in_range(window_sequence, range_start, range_end, record_length):
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
        restriction_sites = []
    return restriction_sites


def _serialize_record_summary(record):
    return {
        "id": record.id,
        "name": record.name,
        "description": record.description,
        "length": len(record.seq),
    }


def _parse_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _parse_json_or_form_payload(request):
    if request.body:
        try:
            return json.loads(request.body.decode("utf-8"))
        except (TypeError, ValueError, json.JSONDecodeError):
            raise ValueError("Invalid JSON payload.")
    return request.POST


def _normalize_feature_attachment_data(payload):
    attach_feature = _parse_bool(payload.get("attach_feature", False))
    save_to_primers = _parse_bool(payload.get("save_to_primers", True))
    record_id = str(payload.get("record_id", "")).strip()
    try:
        feature_start = int(payload.get("feature_start", 0))
        feature_end = int(payload.get("feature_end", 0))
    except (TypeError, ValueError):
        feature_start = 0
        feature_end = 0
    try:
        feature_strand = int(payload.get("feature_strand", 1))
    except (TypeError, ValueError):
        feature_strand = 1
    if feature_strand not in (-1, 1):
        feature_strand = 1
    return {
        "attach_feature": attach_feature,
        "save_to_primers": save_to_primers,
        "record_id": record_id,
        "feature_start": feature_start,
        "feature_end": feature_end,
        "feature_strand": feature_strand,
    }


def _validate_feature_attachment(sequence_file, attachment_data):
    if not attachment_data["attach_feature"]:
        return None

    record_id = attachment_data["record_id"]
    feature_start = attachment_data["feature_start"]
    feature_end = attachment_data["feature_end"]
    feature_strand = attachment_data["feature_strand"]

    if not record_id:
        raise ValueError("record_id is required to attach primer as feature.")
    if feature_start < 1 or feature_end < feature_start:
        raise ValueError("feature_start/feature_end are invalid.")

    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        raise ValueError("Could not parse the selected sequence file.")

    target_record = next((r for r in records if str(r.id) == record_id), None)
    if not target_record:
        raise ValueError("record_id does not exist in this sequence file.")
    if feature_end > len(target_record.seq):
        raise ValueError("Feature coordinates exceed record length.")

    return {
        "record_id": record_id,
        "start": feature_start,
        "end": feature_end,
        "strand": feature_strand,
    }


def _serialize_record_region(record, range_start, range_end, sequence_file):
    record_length = len(record.seq)
    range_start = max(1, min(range_start, record_length))
    range_end = max(1, min(range_end, record_length))
    if range_end < range_start:
        range_start, range_end = range_end, range_start
    window_sequence = str(record.seq[range_start - 1:range_end]).upper()

    combined_features = _extract_record_features(record) + _extract_user_features(
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
        "restriction_sites": _extract_restriction_sites_in_range(window_sequence, range_start, range_end, record_length),
    }


@login_required
def sequencefile_upload(request):
    """
    Upload a FASTA or GenBank sequence file and persist it.
    """

    if request.method == "POST":
        name = request.POST.get("name")
        file = request.FILES.get("file")
        file_type = request.POST.get("file_type")
        description = request.POST.get("description", "")

        if not name or not file or file_type not in ("fasta", "genbank"):
            return render(
                request,
                "core/sequencefile_upload.html",
                {
                    "error": "All required fields must be provided.",
                },
            )

        SequenceFile.objects.create(
            name=name,
            file=file,
            file_type=file_type,
            description=description,
            uploaded_by=request.user,
        )

        return redirect("sequencefile_list")

    return render(
        request,
        "core/sequencefile_upload.html",
    )


@login_required
def primer_binding_analysis_async(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    primer_id = request.POST.get("primer_id")
    sequence_file_id = request.POST.get("sequence_file_id")

    primer = get_object_or_404(
        Primer,
        id=primer_id,
        users=request.user,
    )

    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequence_file_id,
        uploaded_by=request.user,
    )

    task = analyze_primer_binding_task.delay(
        primer_id=primer.id,
        sequence_file_id=sequence_file.id,
        max_mismatches=2,
    )

    return JsonResponse({"task_id": task.id}, status=202)


@login_required
def primer_binding_status(request, task_id):
    result = AsyncResult(task_id)

    if result.state == "FAILURE":
        return JsonResponse(
            {"state": result.state, "error": str(result.info)},
            status=500,
        )

    if result.state == "SUCCESS":
        return JsonResponse(
            {"state": result.state, "result": result.result},
        )

    return JsonResponse({"state": result.state})


@login_required
def sequencefile_list(request):
    """
    List uploaded FASTA / GenBank sequence files for the current user.
    Supports optional filtering by search term and file type.
    """

    qs = SequenceFile.objects.filter(uploaded_by=request.user)

    q = request.GET.get("q")
    if q:
        qs = qs.filter(
            Q(name__icontains=q)
            | Q(description__icontains=q)
        )

    file_type = request.GET.get("type")
    if file_type in ("fasta", "genbank"):
        qs = qs.filter(file_type=file_type)

    order = request.GET.get("order", "uploaded_desc")
    allowed_orders = {
        "name": "name",
        "name_desc": "-name",
        "uploaded": "uploaded_at",
        "uploaded_desc": "-uploaded_at",
    }
    qs = qs.order_by(allowed_orders.get(order, "-uploaded_at"))

    page_obj, query_string = paginate_queryset(request, qs)
    return render(
        request,
        "core/sequencefile_list.html",
        {
            "sequence_files": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )


@login_required
def primer_binding_analysis(request):
    primers = Primer.objects.filter(users=request.user)
    sequence_files = SequenceFile.objects.filter(uploaded_by=request.user)

    preselected_primer = request.GET.get("primer")
    preselected_sequence_file = request.GET.get("sequence_file")

    if preselected_primer:
        preselected_primer = (
            primers.filter(id=preselected_primer)
            .values_list("id", flat=True)
            .first()
        )
        if not preselected_primer:
            messages.error(request, "Primer no longer available")

    if preselected_primer:
        primers = primers.order_by(
            Case(
                When(id=preselected_primer, then=Value(0)),
                default=Value(1),
                output_field=IntegerField(),
            ),
            "primer_name",
        )

    if preselected_sequence_file:
        preselected_sequence_file = (
            sequence_files.filter(id=preselected_sequence_file)
            .values_list("id", flat=True)
            .first()
        )

    if request.method == "POST":
        primer_id = request.POST.get("primer_id")
        sequence_file_id = request.POST.get("sequence_file_id")

        primer = get_object_or_404(
            Primer,
            id=primer_id,
            users=request.user,
        )

        sequence_file = get_object_or_404(
            SequenceFile,
            id=sequence_file_id,
            uploaded_by=request.user,
        )

        hits = analyze_primer_binding(
            primer_sequence=primer.sequence,
            sequence_file=sequence_file,
            max_mismatches=2,
        )

        return render(
            request,
            "core/primer_binding_results.html",
            {
                "primer": primer,
                "sequence_file": sequence_file,
                "hits": hits,
            },
        )

    return render(
        request,
        "core/primer_binding_upload.html",
        {
            "primers": primers,
            "sequence_files": sequence_files,
            "preselected_primer": preselected_primer,
            "preselected_sequence_file": preselected_sequence_file,
        },
    )


@login_required
def sequencefile_linear_view(request, sequencefile_id):
    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequencefile_id,
        uploaded_by=request.user,
    )

    records_payload = []
    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    try:
        for record in records:
            records_payload.append(_serialize_record_summary(record))
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    return render(
        request,
        "core/sequencefile_linear_view.html",
        {
            "sequence_file": sequence_file,
            "records_payload": records_payload,
        },
    )


@login_required
def sequencefile_linear_record_data(request, sequencefile_id):
    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequencefile_id,
        uploaded_by=request.user,
    )

    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        return JsonResponse({"error": "Could not parse the selected sequence file."}, status=400)

    if not records:
        return JsonResponse({"error": "No sequence records were found in this file."}, status=404)

    try:
        record_index = int(request.GET.get("record_index", 0))
    except (TypeError, ValueError):
        return JsonResponse({"error": "record_index must be an integer."}, status=400)

    if record_index < 0 or record_index >= len(records):
        return JsonResponse({"error": "record_index out of range."}, status=400)

    record = records[record_index]
    record_length = len(record.seq)
    try:
        range_start = int(request.GET.get("start", 1))
        range_end = int(request.GET.get("end", min(record_length, 5000)))
    except (TypeError, ValueError):
        return JsonResponse({"error": "start/end must be integers."}, status=400)

    try:
        payload = _serialize_record_region(record, range_start, range_end, sequence_file)
    except Exception:
        return JsonResponse({"error": "Could not parse selected sequence record."}, status=400)

    return JsonResponse(payload)


@login_required
def sequencefile_linear_create_primer(request, sequencefile_id):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequencefile_id,
        uploaded_by=request.user,
    )

    try:
        payload = _parse_json_or_form_payload(request)
    except ValueError as exc:
        return JsonResponse({"error": str(exc)}, status=400)

    primer_name = str(payload.get("primer_name", "")).strip()
    sequence = str(payload.get("sequence", "")).strip()
    overhang_sequence = str(payload.get("overhang_sequence", "")).strip()
    attachment_data = _normalize_feature_attachment_data(payload)

    if not primer_name:
        return JsonResponse({"error": "Primer name is required."}, status=400)
    if not attachment_data["attach_feature"] and not attachment_data["save_to_primers"]:
        return JsonResponse(
            {"error": "Select at least one destination: sequence file or oligostore primers."},
            status=400,
        )

    try:
        sequence = clean_sequence_value(sequence, allow_n=False, max_length=60)
        overhang_sequence = clean_optional_sequence_value(
            overhang_sequence,
            allow_n=False,
        )
    except ValidationError as exc:
        message = exc.messages[0] if exc.messages else "Invalid primer input."
        return JsonResponse({"error": message}, status=400)

    try:
        validated_attachment = _validate_feature_attachment(sequence_file, attachment_data)
    except ValueError as exc:
        return JsonResponse({"error": str(exc)}, status=400)

    primer = None
    if attachment_data["save_to_primers"]:
        primer = Primer.create_with_analysis(
            primer_name=primer_name,
            sequence=sequence,
            overhang_sequence=overhang_sequence,
            user=request.user,
        )

    attached_feature = None
    if validated_attachment:
        attached_feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id=validated_attachment["record_id"],
            start=validated_attachment["start"],
            end=validated_attachment["end"],
            strand=validated_attachment["strand"],
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label=primer.primer_name if primer else primer_name,
            created_by=request.user,
        )

    return JsonResponse(
        {
            "ok": True,
            "primer": (
                {
                    "id": primer.id,
                    "name": primer.primer_name,
                    "sequence": primer.sequence,
                    "overhang_sequence": primer.overhang_sequence or "",
                    "length": primer.length,
                }
                if primer
                else None
            ),
            "attached_feature": (
                {
                    "id": attached_feature.id,
                    "record_id": attached_feature.record_id,
                    "start": attached_feature.start,
                    "end": attached_feature.end,
                    "strand": attached_feature.strand,
                }
                if attached_feature
                else None
            ),
        },
        status=201,
    )
