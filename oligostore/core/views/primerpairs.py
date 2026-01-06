from io import BytesIO

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook
from openpyxl.styles import Font

from ..forms import PrimerPairForm, PrimerPairCreateCombinedForm
from ..models import Primer, PrimerPair
from ..services.primer_analysis import analyze_primer
from ..services.user_assignment import assign_creator
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

            forward = Primer(
                primer_name=form.cleaned_data["forward_name"],
                sequence=form.cleaned_data["forward_sequence"],
            )
            forward = assign_creator(forward, request.user)
            analysis = analyze_primer(forward.sequence)
            forward.length = len(forward.sequence)
            forward.gc_content = analysis["gc_content"]
            forward.tm = analysis["tm"]
            forward.hairpin_dg = analysis["hairpin_dg"]
            forward.self_dimer_dg = analysis["self_dimer_dg"]
            forward.save()

            reverse = Primer(
                primer_name=form.cleaned_data["reverse_name"],
                sequence=form.cleaned_data["reverse_sequence"],
            )
            reverse = assign_creator(reverse, request.user)
            analysis = analyze_primer(reverse.sequence)
            reverse.length = len(reverse.sequence)
            reverse.gc_content = analysis["gc_content"]
            reverse.tm = analysis["tm"]
            reverse.hairpin_dg = analysis["hairpin_dg"]
            reverse.self_dimer_dg = analysis["self_dimer_dg"]

            reverse.save()

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
    sheet = workbook.active
    sheet.title = "Primer Pairs"

    headers = [
        "Pair Name",
        "Forward Name",
        "Forward Sequence",
        "Forward TM",
        "Reverse Name",
        "Reverse Sequence",
        "Reverse TM",
    ]
    sheet.append(headers)
    for cell in sheet[1]:
        cell.font = Font(bold=True)

    for pair in primer_pairs:
        sheet.append(
            [
                pair.name,
                pair.forward_primer.primer_name,
                pair.forward_primer.sequence,
                pair.forward_primer.tm,
                pair.reverse_primer.primer_name,
                pair.reverse_primer.sequence,
                pair.reverse_primer.tm,
            ]
        )

    for column_cells in sheet.columns:
        max_length = 0
        for cell in column_cells:
            cell_value = cell.value
            if cell_value is None:
                continue
            max_length = max(max_length, len(str(cell_value)))
        adjusted_width = min(max_length + 2, 50)
        column_letter = column_cells[0].column_letter
        sheet.column_dimensions[column_letter].width = max(adjusted_width, 12)

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