from django.contrib.auth.decorators import login_required
from django.core.exceptions import ValidationError
from django.db.models import Case, IntegerField, Value, When
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

from ..access import (
    accessible_pcr_products,
    accessible_primers,
    accessible_sequence_files,
)
from ..models import AnalysisJob, SequenceFeature, SequenceFile
from ..forms import clean_optional_sequence_value, clean_sequence_value
from ..services.async_jobs import create_analysis_job, get_owned_job_or_404, mark_job_running
from ..services.creation import create_pcr_product, create_primer_and_optional_feature
from ..services.listing import apply_ordering, apply_search
from ..services.primer_binding import analyze_primer_binding
from ..services.sequence_loader import load_sequences
from ..services.sequence_records import (
    get_sequence_records,
    parse_json_or_form_payload,
    serialize_record_region,
    serialize_record_summary,
)
from ..tasks import analyze_primer_binding_task
from .utils import paginate_queryset


def _allowed_sequence_file_types():
    return {value for value, _label in SequenceFile.FILE_TYPE_CHOICES}


def _get_sequence_records(sequence_file):
    return get_sequence_records(sequence_file, load_sequences)


def _parse_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _parse_optional_int(value, default=None):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


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

@login_required
def sequencefile_upload(request):
    """
    Upload a FASTA, GenBank, or SnapGene sequence file and persist it.
    """
    allowed_file_types = _allowed_sequence_file_types()

    if request.method == "POST":
        name = request.POST.get("name")
        file = request.FILES.get("file")
        file_type = request.POST.get("file_type")
        description = request.POST.get("description", "")

        if not name or not file or file_type not in allowed_file_types:
            return render(
                request,
                "core/sequencefile_upload.html",
                {
                    "error": "All required fields must be provided.",
                },
            )

        created_file = SequenceFile.objects.create(
            name=name,
            file=file,
            file_type=file_type,
            description=description,
            uploaded_by=request.user,
        )
        created_file.users.add(request.user)

        return redirect("sequencefile_list")

    return render(
        request,
        "core/sequencefile_upload.html",
    )


@login_required
def sequencefile_update_type(request, sequencefile_id):
    if request.method != "POST":
        return redirect("sequencefile_list")

    sequence_file = get_object_or_404(
        accessible_sequence_files(request.user),
        id=sequencefile_id,
    )
    file_type = request.POST.get("file_type")
    allowed_file_types = _allowed_sequence_file_types()
    redirect_url = request.POST.get("next") or request.META.get("HTTP_REFERER") or redirect("sequencefile_list").url

    if file_type not in allowed_file_types:
        messages.error(request, "Invalid sequence file type.")
        return redirect(redirect_url)

    if sequence_file.file_type != file_type:
        sequence_file.file_type = file_type
        sequence_file.save(update_fields=["file_type"])
        messages.success(request, f"Updated sequence file type for {sequence_file.name}.")

    return redirect(redirect_url)


@login_required
def primer_binding_analysis_async(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    primer_id = request.POST.get("primer_id")
    sequence_file_id = request.POST.get("sequence_file_id")

    primer = get_object_or_404(accessible_primers(request.user), id=primer_id)

    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequence_file_id)
    job = create_analysis_job(
        owner=request.user,
        job_type=AnalysisJob.TYPE_PRIMER_BINDING,
        primer=primer,
        sequence_file=sequence_file,
    )

    task = analyze_primer_binding_task.delay(
        analysis_job_id=job.id,
        primer_id=primer.id,
        sequence_file_id=sequence_file.id,
        max_mismatches=2,
    )
    mark_job_running(job, task.id)
    return JsonResponse({"job_id": job.id, "task_id": task.id}, status=202)


@login_required
def primer_binding_status(request, task_id):
    job = get_owned_job_or_404(owner=request.user, job_id=task_id)
    status_map = {
        AnalysisJob.STATUS_PENDING: "PENDING",
        AnalysisJob.STATUS_RUNNING: "STARTED",
        AnalysisJob.STATUS_SUCCESS: "SUCCESS",
        AnalysisJob.STATUS_FAILURE: "FAILURE",
    }
    payload = {"state": status_map[job.status], "job_id": job.id}
    if job.status == AnalysisJob.STATUS_FAILURE:
        payload["error"] = job.error_message
        return JsonResponse(payload, status=500)
    if job.status == AnalysisJob.STATUS_SUCCESS:
        payload["result"] = job.result_payload
    return JsonResponse(payload)


@login_required
def sequencefile_list(request):
    """
    List uploaded FASTA / GenBank / SnapGene sequence files for the current user.
    Supports optional filtering by search term and file type.
    """

    qs = accessible_sequence_files(request.user)

    q = request.GET.get("q")
    qs = apply_search(qs, q, ["name", "description"])

    file_type = request.GET.get("type")
    if file_type in _allowed_sequence_file_types():
        qs = qs.filter(file_type=file_type)

    order = request.GET.get("order", "uploaded_desc")
    allowed_orders = {
        "name": "name",
        "name_desc": "-name",
        "uploaded": "uploaded_at",
        "uploaded_desc": "-uploaded_at",
    }
    qs = apply_ordering(qs, order, allowed_orders, "-uploaded_at")

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
def pcrproduct_list(request):
    qs = (
        accessible_pcr_products(request.user)
        .distinct()
        .select_related(
            "sequence_file",
            "forward_primer",
            "reverse_primer",
            "forward_feature",
            "reverse_feature",
        )
    )

    q = request.GET.get("q")
    qs = apply_search(
        qs,
        q,
        [
            "name",
            "record_id",
            "sequence_file__name",
            "forward_primer_label",
            "reverse_primer_label",
            "forward_primer__primer_name",
            "reverse_primer__primer_name",
        ],
    )

    order = request.GET.get("order", "created_desc")
    allowed_orders = {
        "created_desc": "-created_at",
        "created": "created_at",
        "name": "name",
        "name_desc": "-name",
        "length_desc": "-length",
        "length": "length",
    }
    qs = apply_ordering(qs, order, allowed_orders, "-created_at")

    page_obj, query_string = paginate_queryset(request, qs)
    return render(
        request,
        "core/pcrproduct_list.html",
        {
            "pcr_products": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )


@login_required
def primer_binding_analysis(request):
    primers = accessible_primers(request.user)
    sequence_files = accessible_sequence_files(request.user)

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

        primer = get_object_or_404(accessible_primers(request.user), id=primer_id)

        sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequence_file_id)

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
    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    records_payload = []
    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    try:
        for record in records:
            records_payload.append(serialize_record_summary(record))
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    initial_pcr_product = None
    requested_pcr_product_id = _parse_optional_int(request.GET.get("pcr_product"))
    if requested_pcr_product_id:
        pcr_product = (
            accessible_pcr_products(request.user)
            .select_related("forward_feature", "reverse_feature", "forward_primer", "reverse_primer")
            .filter(id=requested_pcr_product_id, sequence_file=sequence_file)
            .first()
        )
        if pcr_product:
            initial_pcr_product = {
                "id": pcr_product.id,
                "name": pcr_product.name,
                "record_id": pcr_product.record_id,
                "start": pcr_product.start,
                "end": pcr_product.end,
                "length": pcr_product.length,
                "sequence": pcr_product.sequence,
                "forward_feature_id": pcr_product.forward_feature_id,
                "reverse_feature_id": pcr_product.reverse_feature_id,
                "forward_primer_id": pcr_product.forward_primer_id,
                "reverse_primer_id": pcr_product.reverse_primer_id,
                "forward_primer_label": pcr_product.forward_primer_label,
                "reverse_primer_label": pcr_product.reverse_primer_label,
            }

    return render(
        request,
        "core/sequencefile_linear_view.html",
        {
            "sequence_file": sequence_file,
            "records_payload": records_payload,
            "initial_pcr_product": initial_pcr_product,
        },
    )


@login_required
def sequencefile_circular_view(request, sequencefile_id):
    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    records_payload = []
    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    try:
        for record in records:
            records_payload.append(serialize_record_summary(record))
    except Exception:
        messages.error(request, "Could not parse the selected sequence file.")
        return redirect("sequencefile_list")

    return render(
        request,
        "core/sequencefile_circular_view.html",
        {
            "sequence_file": sequence_file,
            "records_payload": records_payload,
        },
    )


@login_required
def sequencefile_linear_record_data(request, sequencefile_id):
    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

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
        payload = serialize_record_region(record, range_start, range_end, sequence_file)
    except Exception:
        return JsonResponse({"error": "Could not parse selected sequence record."}, status=400)

    return JsonResponse(payload)


@login_required
def sequencefile_linear_create_primer(request, sequencefile_id):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

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

    primer, attached_feature = create_primer_and_optional_feature(
        user=request.user,
        sequence_file=sequence_file,
        primer_name=primer_name,
        sequence=sequence,
        overhang_sequence=overhang_sequence,
        save_to_primers=attachment_data["save_to_primers"],
        feature_attachment=validated_attachment,
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


@login_required
def sequencefile_linear_save_pcr_product(request, sequencefile_id):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    try:
        payload = _parse_json_or_form_payload(request)
    except ValueError as exc:
        return JsonResponse({"error": str(exc)}, status=400)

    name = str(payload.get("name", "")).strip()
    record_id = str(payload.get("record_id", "")).strip()
    start = _parse_optional_int(payload.get("start"), 0)
    end = _parse_optional_int(payload.get("end"), 0)
    sequence = str(payload.get("sequence", "")).strip().upper()
    forward_primer_label = str(payload.get("forward_primer_label", "")).strip()
    reverse_primer_label = str(payload.get("reverse_primer_label", "")).strip()
    forward_primer_id = _parse_optional_int(payload.get("forward_primer_id"))
    reverse_primer_id = _parse_optional_int(payload.get("reverse_primer_id"))
    forward_feature_id = _parse_optional_int(payload.get("forward_feature_id"))
    reverse_feature_id = _parse_optional_int(payload.get("reverse_feature_id"))

    if not record_id:
        return JsonResponse({"error": "record_id is required."}, status=400)
    if start < 1 or end < start:
        return JsonResponse({"error": "start/end are invalid."}, status=400)

    try:
        sequence = clean_sequence_value(sequence, allow_n=True)
    except ValidationError as exc:
        message = exc.messages[0] if exc.messages else "Invalid PCR product sequence."
        return JsonResponse({"error": message}, status=400)

    try:
        records = _get_sequence_records(sequence_file)
    except Exception:
        return JsonResponse({"error": "Could not parse the selected sequence file."}, status=400)

    target_record = next((r for r in records if str(r.id) == record_id), None)
    if not target_record:
        return JsonResponse({"error": "record_id does not exist in this sequence file."}, status=400)
    if end > len(target_record.seq):
        return JsonResponse({"error": "PCR product coordinates exceed record length."}, status=400)

    expected_sequence = str(target_record.seq[start - 1:end]).upper()
    if expected_sequence != sequence:
        return JsonResponse({"error": "PCR product sequence does not match the selected sequence record."}, status=400)

    if not name:
        name = f"{sequence_file.name}:{record_id}:{start}-{end}"

    primer_queryset = accessible_primers(request.user)
    forward_primer = primer_queryset.filter(id=forward_primer_id).first() if forward_primer_id else None
    reverse_primer = primer_queryset.filter(id=reverse_primer_id).first() if reverse_primer_id else None
    forward_feature = (
        SequenceFeature.objects.filter(id=forward_feature_id, sequence_file=sequence_file).first()
        if forward_feature_id else None
    )
    reverse_feature = (
        SequenceFeature.objects.filter(id=reverse_feature_id, sequence_file=sequence_file).first()
        if reverse_feature_id else None
    )

    product = create_pcr_product(
        user=request.user,
        sequence_file=sequence_file,
        name=name,
        record_id=record_id,
        start=start,
        end=end,
        sequence=sequence,
        forward_primer=forward_primer,
        reverse_primer=reverse_primer,
        forward_feature=forward_feature,
        reverse_feature=reverse_feature,
        forward_primer_label=forward_primer_label,
        reverse_primer_label=reverse_primer_label,
    )

    return JsonResponse(
        {
            "ok": True,
            "pcr_product": {
                "id": product.id,
                "name": product.name,
                "record_id": product.record_id,
                "start": product.start,
                "end": product.end,
                "length": product.length,
            },
        },
        status=201,
    )


@login_required
def sequencefile_linear_delete_primer(request, sequencefile_id):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    try:
        payload = _parse_json_or_form_payload(request)
    except ValueError as exc:
        return JsonResponse({"error": str(exc)}, status=400)

    try:
        feature_id = int(payload.get("feature_id", 0))
    except (TypeError, ValueError):
        feature_id = 0
    delete_primer = _parse_bool(payload.get("delete_primer", False))

    if feature_id < 1:
        return JsonResponse({"error": "feature_id is required."}, status=400)

    feature = get_object_or_404(
        SequenceFeature,
        id=feature_id,
        sequence_file=sequence_file,
    )

    primer = feature.primer
    deleted_primer_id = None
    if delete_primer:
        if not primer:
            return JsonResponse({"error": "This sequence feature is not linked to an oligostore primer."}, status=400)
        if primer.creator_id != request.user.id:
            return JsonResponse({"error": "You do not have permission to delete this primer from oligostore."}, status=403)
        deleted_primer_id = primer.id

    deleted_feature_id = feature.id
    feature.delete()

    if delete_primer and primer:
        primer.delete()

    return JsonResponse(
        {
            "ok": True,
            "deleted_feature_id": deleted_feature_id,
            "deleted_primer_id": deleted_primer_id,
        }
    )
