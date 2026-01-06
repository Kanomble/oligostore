from io import BytesIO

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook
from openpyxl.styles import Font
from ..services.export_helpers import build_primer_worksheet
from ..forms import PrimerForm
from ..models import Primer
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
            Primer.create_with_analysis(
                primer_name=form.cleaned_data["primer_name"],
                sequence=form.cleaned_data["sequence"],
                overhang_sequence=form.cleaned_data.get("overhang_sequence", ""),
                user=request.user,
            )
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
    build_primer_worksheet(workbook, primers)

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