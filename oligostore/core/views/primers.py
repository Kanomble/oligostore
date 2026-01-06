from io import BytesIO
from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from django.http import HttpResponse
from django.shortcuts import render, redirect, get_object_or_404
from openpyxl import Workbook
from ..services.export_helpers import build_primer_worksheet
from ..forms import (
    PrimerExcelColumnMapForm,
    PrimerExcelUploadForm,
    PrimerForm,
    clean_optional_sequence_value,
    clean_sequence_value,
)
from ..models import Primer
from .utils import paginate_queryset
from django.core.exceptions import ValidationError
import pandas as pd
import json

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
def delete_selected_primers(request):
    if request.method != "POST":
        return HttpResponse("POST only", status=405)

    primer_ids = request.POST.getlist("primer_ids")
    if not primer_ids:
        messages.error(request, "Select at least one primer to delete.")
        return redirect("primer_list")

    primers = Primer.objects.filter(id__in=primer_ids, creator=request.user)
    deletable_count = primers.count()
    if deletable_count == 0:
        messages.error(
            request, "No primers matched your selection or you lack delete access."
        )
        return redirect("primer_list")

    primers.delete()
    messages.success(request, f"Deleted {deletable_count} primer(s).")
    if deletable_count < len(primer_ids):
        messages.warning(
            request, "Some primers were not deleted because you do not own them."
        )
    return redirect("primer_list")

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

@login_required(login_url="login")
def primer_import_excel(request):
    column_session_key = "primer_excel_columns"
    row_session_key = "primer_excel_rows"
    columns = request.session.get(column_session_key, [])
    rows = request.session.get(row_session_key, [])
    upload_form = PrimerExcelUploadForm()
    map_form = PrimerExcelColumnMapForm(columns=columns) if columns else None

    if request.method == "POST":
        if "upload_excel" in request.POST:
            upload_form = PrimerExcelUploadForm(request.POST, request.FILES)
            if upload_form.is_valid():
                excel_file = upload_form.cleaned_data["excel_file"]
                try:
                    df = pd.read_excel(excel_file)
                except Exception:
                    messages.error(request, "Unable to read that Excel file.")
                else:
                    if df.empty or not list(df.columns):
                        messages.error(request, "The Excel file does not contain any rows.")
                    else:
                        columns = [str(col) for col in df.columns]
                        rows = df.fillna("").to_dict(orient="records")
                        request.session[column_session_key] = columns
                        request.session[row_session_key] = rows
                        map_form = PrimerExcelColumnMapForm(columns=columns)
        elif "map_columns" in request.POST:
            if not columns or not rows:
                messages.error(request, "Upload an Excel file before mapping columns.")
            else:
                map_form = PrimerExcelColumnMapForm(request.POST, columns=columns)
                if map_form.is_valid():
                    name_col = map_form.cleaned_data["name_column"]
                    sequence_col = map_form.cleaned_data["sequence_column"]
                    overhang_col = map_form.cleaned_data.get("overhang_column") or ""
                    created = 0
                    skipped = 0
                    errors = []

                    edited_rows = request.POST.get("edited_rows")
                    if edited_rows:
                        try:
                            parsed_rows = json.loads(edited_rows)
                        except json.JSONDecodeError:
                            parsed_rows = None
                            messages.error(
                                request, "Preview data could not be read. Please try again."
                            )
                    else:
                        parsed_rows = None

                    source_rows = parsed_rows if isinstance(parsed_rows, list) else rows

                    for idx, row in enumerate(source_rows, start=1):
                        if isinstance(row, dict) and parsed_rows is not None:
                            name = str(row.get("name", "")).strip()
                            sequence = str(row.get("sequence", "")).strip()
                            overhang = str(row.get("overhang", "")).strip()
                        else:
                            name = str(row.get(name_col, "")).strip()
                            sequence = str(row.get(sequence_col, "")).strip()
                            overhang = (
                                str(row.get(overhang_col, "")).strip()
                                if overhang_col
                                else ""
                            )

                        if not name or not sequence:
                            skipped += 1
                            continue

                        try:
                            sequence = clean_sequence_value(sequence, allow_n=False)
                            overhang = clean_optional_sequence_value(
                                overhang, allow_n=False
                            )
                        except ValidationError as exc:
                            errors.append(
                                f"Row {idx}: {exc.messages[0] if exc.messages else 'Invalid sequence.'}"
                            )
                            continue

                        Primer.create_with_analysis(
                            primer_name=name,
                            sequence=sequence,
                            overhang_sequence=overhang,
                            user=request.user,
                        )
                        created += 1

                    request.session.pop(column_session_key, None)
                    request.session.pop(row_session_key, None)

                    if errors:
                        for error in errors[:5]:
                            messages.error(request, error)
                        if len(errors) > 5:
                            messages.error(
                                request,
                                f"{len(errors) - 5} additional rows were skipped due to errors.",
                            )

                    if skipped:
                        messages.warning(
                            request, f"Skipped {skipped} row(s) without required data."
                        )

                    if created:
                        messages.success(
                            request, f"Imported {created} primer(s) successfully."
                        )
                        return redirect("primer_list")
                    messages.error(request, "No primers were imported.")

    return render(
        request,
        "core/primer_import.html",
        {
            "upload_form": upload_form,
            "map_form": map_form,
            "columns": columns,
            "rows": rows,
        },
    )