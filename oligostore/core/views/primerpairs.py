from io import BytesIO
from django import forms
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook

from ..access import accessible_primer_pairs, accessible_primers, editable_primer_pairs
from ..forms import (
    PCRProductDiscoveryForm,
    PrimerPairForm,
    PrimerPairCreateCombinedForm,
)
from ..models import Primer, PrimerPair
from ..models import AnalysisJob
from ..services.async_jobs import create_analysis_job, mark_job_running
from ..services.creation import create_owned_primer_pair, create_primer_pair_with_new_primers
from ..services.export_helpers import build_primerpair_worksheet
from ..services.listing import apply_ordering, apply_search
from ..tasks import analyze_primerpair_products_task
from .utils import paginate_queryset


@login_required(login_url="login")
def primerpair_list(request):
    primer_pairs = accessible_primer_pairs(request.user)
    q = request.GET.get("q")
    primer_pairs = apply_search(
        primer_pairs,
        q,
        [
            "name",
            "forward_primer__primer_name",
            "forward_primer__sequence",
            "reverse_primer__primer_name",
            "reverse_primer__sequence",
        ],
    )

    order = request.GET.get("order", "name")
    allowed_orders = {
        "name": "name",
        "name_desc": "-name",
        "forward_name": "forward_primer__primer_name",
        "forward_name_desc": "-forward_primer__primer_name",
        "reverse_name": "reverse_primer__primer_name",
        "reverse_name_desc": "-reverse_primer__primer_name",
        "forward_tm": "forward_primer__tm",
        "forward_tm_desc": "-forward_primer__tm",
        "reverse_tm": "reverse_primer__tm",
        "reverse_tm_desc": "-reverse_primer__tm",
    }
    primer_pairs = apply_ordering(primer_pairs, order, allowed_orders, "name")

    page_obj, query_string = paginate_queryset(request, primer_pairs)
    return render(
        request,
        "core/primerpair_list.html",
        {
            "primer_pairs": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )


@login_required(login_url="login")
def primerpair_products(request):
    initial = {}
    preselected_pair = request.GET.get("primer_pair")
    if preselected_pair:
        initial["primer_pair"] = (
            accessible_primer_pairs(request.user)
            .filter(id=preselected_pair)
            .values_list("id", flat=True)
            .first()
        )

    form = PCRProductDiscoveryForm(user=request.user, initial=initial)

    return render(
        request,
        "core/primerpair_products.html",
        {
            "form": form,
        },
    )


@login_required(login_url="login")
def primerpair_products_async(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST required."}, status=405)

    form = PCRProductDiscoveryForm(request.POST, user=request.user)
    if not form.is_valid():
        return JsonResponse({"error": "Invalid PCR product request."}, status=400)

    primer_pair = form.cleaned_data["primer_pair"]
    sequence_file = form.cleaned_data["sequence_file"]
    job = create_analysis_job(
        owner=request.user,
        job_type=AnalysisJob.TYPE_PCR_PRODUCT_DISCOVERY,
        primer_pair=primer_pair,
        sequence_file=sequence_file,
    )
    task = analyze_primerpair_products_task.delay(
        analysis_job_id=job.id,
        primer_pair_id=primer_pair.id,
        sequence_file_id=sequence_file.id,
        max_mismatches=form.cleaned_data["max_mismatches"],
        block_3prime_mismatch=form.cleaned_data["block_3prime_mismatch"],
    )
    mark_job_running(job, task.id)
    return JsonResponse({"job_id": job.id, "task_id": task.id}, status=202)


@login_required(login_url="login")
def primerpair_create(request):

    all_primers = accessible_primers(request.user)
    primers = all_primers

    q = request.GET.get("q")
    primers = apply_search(primers, q, ["primer_name", "sequence"])

    order = request.GET.get("order", "created_desc")
    allowed_orders = {
        "created_desc": "-created_at",
        "created": "created_at",
        "name": "primer_name",
        "name_desc": "-primer_name",
        "length_desc": "-length",
        "length": "length",
        "gc_desc": "-gc_content",
        "gc": "gc_content",
        "tm_desc": "-tm",
        "tm": "tm",
    }
    primers = apply_ordering(primers, order, allowed_orders, "-created_at")

    page_obj, query_string = paginate_queryset(request, primers)

    selected_forward = None
    selected_reverse = None

    if request.method == "POST":
        form = PrimerPairForm(request.POST, user=request.user)
        if form.is_valid():
            pair = form.save(commit=False)
            create_owned_primer_pair(
                name=pair.name,
                forward_primer=pair.forward_primer,
                reverse_primer=pair.reverse_primer,
                user=request.user,
            )
            messages.success(request, "Primer pair created.")
            form = PrimerPairForm(user=request.user)
        else:
            forward_id = form.data.get("forward_primer")
            reverse_id = form.data.get("reverse_primer")
            if forward_id:
                selected_forward = all_primers.filter(id=forward_id).first()
            if reverse_id:
                selected_reverse = all_primers.filter(id=reverse_id).first()

    else:
        form = PrimerPairForm(user=request.user)

    form.fields["forward_primer"].widget = forms.HiddenInput()
    form.fields["reverse_primer"].widget = forms.HiddenInput()

    return render(
        request,
        "core/primerpair_form.html",
        {
            "form": form,
            "primers": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
            "selected_forward": selected_forward,
            "selected_reverse": selected_reverse,
            "search_query": q or "",
            "order": order,
        },
    )



@login_required
def primerpair_combined_create(request):
    if request.method == "POST":
        form = PrimerPairCreateCombinedForm(request.POST)
        if form.is_valid():
            _, _, pair = create_primer_pair_with_new_primers(
                pair_name=form.cleaned_data["pair_name"],
                forward_name=form.cleaned_data["forward_name"],
                forward_sequence=form.cleaned_data["forward_sequence"],
                reverse_name=form.cleaned_data["reverse_name"],
                reverse_sequence=form.cleaned_data["reverse_sequence"],
                forward_overhang=form.cleaned_data.get("forward_overhang", ""),
                reverse_overhang=form.cleaned_data.get("reverse_overhang", ""),
                user=request.user,
            )
            messages.success(request, "Primer pair created.")
            form = PrimerPairCreateCombinedForm()

    else:
        form = PrimerPairCreateCombinedForm()

    return render(request, "core/primerpair_combined_form.html", {"form": form})


@login_required
def primerpair_delete(request, primerpair_id):
    pair = get_object_or_404(editable_primer_pairs(request.user), id=primerpair_id)
    if pair.creator != request.user:
        raise PermissionDenied("You are not the creator of this primer.")
    pair.delete()
    return redirect("primerpair_list")


@login_required(login_url="login")
def download_selected_primerpairs(request):
    if request.method != "POST":
        return HttpResponse("POST only", status=405)

    primerpair_ids = request.POST.getlist("primerpair_ids")
    if not primerpair_ids:
        messages.error(request, "Select at least one primer pair to download.")
        return redirect("primerpair_list")

    primer_pairs = (
        accessible_primer_pairs(request.user).filter(id__in=primerpair_ids)
        .select_related("forward_primer", "reverse_primer")
        .order_by("name")
    )
    if not primer_pairs.exists():
        messages.error(request, "No primer pairs matched your selection.")
        return redirect("primerpair_list")

    workbook = Workbook()
    build_primerpair_worksheet(workbook, primer_pairs)

    output = BytesIO()
    workbook.save(output)
    output.seek(0)

    response = HttpResponse(
        output.getvalue(),
        content_type=(
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
        ),
    )
    response["Content-Disposition"] = 'attachment; filename="primer_pairs.xlsx"'
    return response
