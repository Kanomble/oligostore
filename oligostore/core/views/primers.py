from io import BytesIO

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook
from openpyxl.styles import Font

from ..forms import PrimerForm
from ..models import Primer
from ..services.primer_analysis import analyze_primer
from ..services.user_assignment import assign_creator
from .utils import paginate_queryset


@login_required(login_url="login")
def primer_list(request):
    primers = Primer.objects.filter(users=request.user)

    q = request.GET.get("q")
    if q:
        primers = primers.filter(
            Q(primer_name__icontains=q)
            | Q(sequence__icontains=q)
            | Q(creator__username__icontains=q)
        )

    order = request.GET.get("order", "created_desc")
    allowed_orders = {
        "name": "primer_name",
        "name_desc": "-primer_name",
        "created": "created_at",
        "created_desc": "-created_at",
        "length": "length",
        "length_desc": "-length",
        "gc": "gc_content",
        "gc_desc": "-gc_content",
        "tm": "tm",
        "tm_desc": "-tm",
    }
    primers = primers.order_by(allowed_orders.get(order, "-created_at"))
    page_obj, query_string = paginate_queryset(request, primers)
    return render(
        request,
        "core/primer_list.html",
        {
            "primers": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )


@login_required(login_url="login")
def primer_create(request):
    if request.method == "POST":
        form = PrimerForm(request.POST)
        if form.is_valid():
            primer = form.save(commit=False)
            primer = assign_creator(primer, user=request.user)
            primer.length = len(primer.sequence)
            analysis = analyze_primer(primer.sequence)
            primer.gc_content = analysis["gc_content"]
            primer.tm = analysis["tm"]
            primer.hairpin_dg = analysis["hairpin_dg"]
            primer.self_dimer_dg = analysis["self_dimer_dg"]

            primer.save()
            return redirect("primer_list")
    else:
        form = PrimerForm()

    return render(request, "core/primer_create.html", {"form": form})


@login_required(login_url="login")
def primer_delete(request, primer_id):
    primer = get_object_or_404(Primer, id=primer_id)
    if primer.creator != request.user:
        raise PermissionDenied("You are not the creator of this primer.")

    primer.delete()
    return redirect("primer_list")


@login_required(login_url="login")
def download_selected_primers(request):
    if request.method != "POST":
        return HttpResponse("POST only", status=405)

    primer_ids = request.POST.getlist("primer_ids")
    if not primer_ids:
        messages.error(request, "Select at least one primer to download.")
        return redirect("primer_list")

    primers = (
        Primer.objects.filter(users=request.user, id__in=primer_ids)
        .order_by("primer_name")
    )
    if not primers.exists():
        messages.error(request, "No primers matched your selection.")
        return redirect("primer_list")

    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Primers"

    headers = [
        "Name",
        "Sequence",
        "Length",
        "GC Content",
        "Temperature",
        "Hairpin",
        "Self Dimer",
        "Creator",
        "Created",
    ]
    sheet.append(headers)
    for cell in sheet[1]:
        cell.font = Font(bold=True)

    for primer in primers:
        sheet.append(
            [
                primer.primer_name,
                primer.sequence,
                primer.length,
                primer.gc_content,
                primer.tm,
                primer.hairpin_dg,
                primer.self_dimer_dg,
                str(primer.creator),
                primer.created_at.strftime("%Y-%m-%d %H:%M"),
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
    response["Content-Disposition"] = 'attachment; filename="primers.xlsx"'
    return response