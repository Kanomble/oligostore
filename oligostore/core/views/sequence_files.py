from celery.result import AsyncResult
from django.contrib.auth.decorators import login_required
from django.db.models import Case, IntegerField, Q, Value, When
from django.http import JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from django.contrib import messages

from ..models import Primer, SequenceFile
from ..services.primer_binding import analyze_primer_binding
from ..tasks import analyze_primer_binding_task
from .utils import paginate_queryset

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