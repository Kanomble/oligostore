from io import BytesIO

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook

from ..forms import PrimerPairForm, PrimerPairCreateCombinedForm
from ..models import Primer, PrimerPair
from ..services.user_assignment import assign_creator
from ..services.export_helpers import build_primerpair_worksheet
from .utils import paginate_queryset


@login_required(login_url="login")
def primerpair_list(request):
    primer_pairs = PrimerPair.objects.filter(users=request.user)

    q = request.GET.get("q")
    if q:
        primer_pairs = primer_pairs.filter(
            Q(name__icontains=q)
            | Q(forward_primer__primer_name__icontains=q)
            | Q(forward_primer__sequence__icontains=q)
            | Q(reverse_primer__primer_name__icontains=q)
            | Q(reverse_primer__sequence__icontains=q)
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
    primer_pairs = primer_pairs.order_by(allowed_orders.get(order, "name"))

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
def primerpair_create(request):
    if request.method == "POST":
        form = PrimerPairForm(request.POST, user=request.user)
        if form.is_valid():
            pair = form.save(commit=False)
            pair = assign_creator(pair, request.user)
            pair.save()

            return redirect("primerpair_list")
    else:
        form = PrimerPairForm(user=request.user)

    return render(request, "core/primerpair_form.html", {"form": form})


@login_required
def primerpair_combined_create(request):
    if request.method == "POST":
        form = PrimerPairCreateCombinedForm(request.POST)
        if form.is_valid():

            forward = Primer.create_with_analysis(
                primer_name=form.cleaned_data["forward_name"],
                sequence=form.cleaned_data["forward_sequence"],
                overhang_sequence=form.cleaned_data.get("forward_overhang", ""),
                user=request.user,
            )

            reverse = Primer.create_with_analysis(
                primer_name=form.cleaned_data["reverse_name"],
                sequence=form.cleaned_data["reverse_sequence"],
                overhang_sequence=form.cleaned_data.get("reverse_overhang", ""),
                user=request.user,
            )

            pair = PrimerPair(
                name=form.cleaned_data["pair_name"],
                forward_primer=forward,
                reverse_primer=reverse,
            )
            pair = assign_creator(pair, request.user)
            pair.save()
            return redirect("primerpair_list")

    else:
        form = PrimerPairCreateCombinedForm()

    return render(request, "core/primerpair_combined_form.html", {"form": form})


@login_required
def primerpair_delete(request, primerpair_id):
    pair = get_object_or_404(PrimerPair, id=primerpair_id)
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
        PrimerPair.objects.filter(users=request.user, id__in=primerpair_ids)
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