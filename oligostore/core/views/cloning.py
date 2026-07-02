from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render

from ..access import (
    accessible_cloning_constructs,
    accessible_pcr_products,
    accessible_sequence_files,
    editable_cloning_constructs,
)
from ..forms import (
    CloningConstructAssetForm,
    CloningConstructAssemblyForm,
    CloningConstructSequenceFileForm,
)
from ..services.cloning import (
    build_cloning_assembly_visual_preview,
    build_pcr_product_asset_choice,
    build_cloning_construct_detail_display,
    build_sequence_file_asset_choice,
    preview_cloning_construct,
    resolve_cloning_assets,
    save_cloning_construct,
)
from ..services.cloning_exports import (
    build_cloning_construct_record,
    export_cloning_construct_genbank,
    save_cloning_construct_sequence_file,
)
from ..services.listing import apply_ordering, apply_search
from ..services.sequence_records import extract_record_features
from .utils import paginate_queryset


def _build_construct_asset_form_values(construct):
    if construct.vector_source_type == construct.SOURCE_SEQUENCE_FILE:
        vector_asset = build_sequence_file_asset_choice(
            sequence_file_id=construct.vector_sequence_file_id,
            record_id=construct.vector_record_id,
        ).encoded_value
    elif construct.vector_source_type == construct.SOURCE_TEMPLATE:
        vector_asset = f"{construct.SOURCE_TEMPLATE}:{construct.vector_template_name}:{construct.vector_record_id}"
    else:
        vector_asset = build_pcr_product_asset_choice(
            pcr_product_id=construct.vector_pcr_product_id,
        ).encoded_value

    if construct.insert_source_type == construct.SOURCE_SEQUENCE_FILE:
        insert_asset = build_sequence_file_asset_choice(
            sequence_file_id=construct.insert_sequence_file_id,
            record_id=construct.insert_record_id,
        ).encoded_value
    elif construct.insert_source_type == construct.SOURCE_TEMPLATE:
        insert_asset = f"{construct.SOURCE_TEMPLATE}:{construct.insert_template_name}:{construct.insert_record_id}"
    else:
        insert_asset = build_pcr_product_asset_choice(
            pcr_product_id=construct.insert_pcr_product_id,
        ).encoded_value

    return {
        "name": construct.name,
        "description": construct.description,
        "assembly_strategy": construct.assembly_strategy,
        "vector_asset": vector_asset,
        "insert_asset": insert_asset,
        "left_enzyme": construct.left_enzyme,
        "right_enzyme": construct.right_enzyme,
        "is_circular": "1" if getattr(construct, "is_circular", False) else "0",
        "selected_enzymes": [
            enzyme_name
            for enzyme_name in (construct.left_enzyme, construct.right_enzyme)
            if enzyme_name
        ],
        "vector_fragment_index": "" if construct.vector_fragment_index is None else str(construct.vector_fragment_index),
        "insert_fragment_index": "" if construct.insert_fragment_index is None else str(construct.insert_fragment_index),
    }


def _build_construct_sequence_file_form(construct, data=None):
    if data is None:
        return CloningConstructSequenceFileForm(construct=construct)
    return CloningConstructSequenceFileForm(data, construct=construct)


def _get_accessible_construct(request, construct_id):
    return get_object_or_404(
        accessible_cloning_constructs(request.user).select_related(
            "vector_sequence_file",
            "vector_pcr_product",
            "vector_pcr_product__sequence_file",
            "insert_sequence_file",
            "insert_pcr_product",
            "insert_pcr_product__sequence_file",
        ),
        id=construct_id,
    )


def _build_sequence_lines(sequence, *, line_length=60, group_length=10):
    lines = []
    normalized_sequence = str(sequence or "").upper()
    for offset in range(0, len(normalized_sequence), line_length):
        chunk = normalized_sequence[offset : offset + line_length]
        groups = [chunk[index : index + group_length] for index in range(0, len(chunk), group_length)]
        lines.append(
            {
                "start": offset + 1,
                "end": offset + len(chunk),
                "text": " ".join(groups),
            }
        )
    return lines


def _build_construct_linear_context(construct):
    record = build_cloning_construct_record(construct)
    feature_rows = []
    for feature in extract_record_features(record):
        start = int(feature.get("start", 1))
        end = int(feature.get("end", start))
        feature_rows.append(
            {
                "label": feature.get("label", ""),
                "type": feature.get("type", ""),
                "start": start,
                "end": end,
                "length": max(0, end - start + 1),
                "strand": feature.get("strand"),
                "description": feature.get("description", ""),
            }
        )
    feature_rows.sort(key=lambda feature: (feature["start"], feature["end"], feature["label"]))
    sequence = str(record.seq).upper()
    return {
        "record": {
            "id": record.id,
            "name": record.name,
            "description": record.description,
            "sequence": sequence,
            "length": len(sequence),
            "features": feature_rows,
            "annotations": dict(getattr(record, "annotations", {}) or {}),
        },
        "sequence_lines": _build_sequence_lines(sequence),
    }


def _build_construct_detail_context(request, construct, *, save_sequence_file_form=None):
    return {
        "construct": construct,
        "detail_display": build_cloning_construct_detail_display(construct),
        "can_delete": construct.creator_id == request.user.id,
        "save_sequence_file_form": save_sequence_file_form or _build_construct_sequence_file_form(construct),
    }


def _build_construct_linear_page_context(request, construct, *, save_sequence_file_form=None):
    detail_display = build_cloning_construct_detail_display(construct)
    linear_context = _build_construct_linear_context(construct)
    return {
        "construct": construct,
        "detail_display": detail_display,
        "construct_record": linear_context["record"],
        "sequence_lines": linear_context["sequence_lines"],
        "can_delete": construct.creator_id == request.user.id,
        "save_sequence_file_form": save_sequence_file_form or _build_construct_sequence_file_form(construct),
    }


def _form_field_value(form, field_name):
    try:
        value = form[field_name].value()
    except KeyError:
        value = None
    return value if value is not None else ""


def _form_field_values(form, field_name):
    value = _form_field_value(form, field_name)
    if isinstance(value, (list, tuple)):
        return list(value)
    if value in (None, ""):
        return []
    return [value]


def _build_assembly_visual_preview(assembly_form, assets):
    if assembly_form is None or assets is None:
        return None
    selected_left_enzyme = _form_field_value(assembly_form, "left_enzyme")
    selected_right_enzyme = _form_field_value(assembly_form, "right_enzyme")
    return build_cloning_assembly_visual_preview(
        vector_asset=assets.vector_asset,
        insert_asset=assets.insert_asset,
        selected_left_enzyme=selected_left_enzyme,
        selected_right_enzyme=selected_right_enzyme,
        map_enzyme_names=_form_field_values(assembly_form, "selected_enzymes"),
        vector_fragment_index=_form_field_value(assembly_form, "vector_fragment_index"),
        insert_fragment_index=_form_field_value(assembly_form, "insert_fragment_index"),
    )


@login_required
def cloning_asset_list(request):
    sequence_files = accessible_sequence_files(request.user)
    pcr_products = accessible_pcr_products(request.user).select_related("sequence_file")
    q = request.GET.get("q")
    sequence_files = apply_search(sequence_files, q, ["name", "description"])
    pcr_products = apply_search(
        pcr_products,
        q,
        ["name", "record_id", "sequence_file__name", "forward_primer_label", "reverse_primer_label"],
    )
    return render(
        request,
        "core/cloning_asset_list.html",
        {
            "sequence_files": sequence_files.order_by("name"),
            "pcr_products": pcr_products.order_by("name"),
        },
    )


@login_required
def cloning_construct_list(request):
    constructs = accessible_cloning_constructs(request.user).select_related(
        "vector_sequence_file",
        "vector_pcr_product",
        "insert_sequence_file",
        "insert_pcr_product",
    )
    q = request.GET.get("q")
    constructs = apply_search(
        constructs,
        q,
        [
            "name",
            "description",
            "vector_sequence_file__name",
            "vector_template_name",
            "vector_pcr_product__name",
            "insert_sequence_file__name",
            "insert_template_name",
            "insert_pcr_product__name",
            "left_enzyme",
            "right_enzyme",
        ],
    )
    order = request.GET.get("order", "created_desc")
    allowed_orders = {
        "created_desc": "-created_at",
        "created": "created_at",
        "name": "name",
        "name_desc": "-name",
    }
    constructs = apply_ordering(constructs, order, allowed_orders, "-created_at")
    page_obj, query_string = paginate_queryset(request, constructs)
    return render(
        request,
        "core/cloning_construct_list.html",
        {"constructs": page_obj, "page_obj": page_obj, "query_string": query_string},
    )


@login_required
def cloning_construct_create(request):
    asset_form = CloningConstructAssetForm(user=request.user)
    assembly_form = None
    preview_data = None
    visual_preview = None
    assets = None
    review_construct = None

    review_construct_id = request.GET.get("review_construct")
    if request.method == "GET" and review_construct_id:
        review_construct = get_object_or_404(
            accessible_cloning_constructs(request.user).select_related(
                "vector_sequence_file",
                "vector_pcr_product",
                "insert_sequence_file",
                "insert_pcr_product",
            ),
            id=review_construct_id,
            is_valid=False,
        )
        initial = _build_construct_asset_form_values(review_construct)
        asset_form = CloningConstructAssetForm(initial=initial, user=request.user)
        assembly_form = CloningConstructAssemblyForm(initial=initial, user=request.user)
        try:
            assets = resolve_cloning_assets(
                user=request.user,
                vector_asset_choice=initial["vector_asset"],
                insert_asset_choice=initial["insert_asset"],
            )
        except ValueError as exc:
            messages.error(request, str(exc))
            return redirect("cloning_construct_detail", construct_id=review_construct.id)

    if request.method == "POST":
        step = request.POST.get("step")
        if not step and request.POST.get("left_enzyme") and request.POST.get("right_enzyme"):
            step = "save"
        if step == "assets":
            asset_form = CloningConstructAssetForm(request.POST, user=request.user)
            if asset_form.is_valid():
                try:
                    assets = resolve_cloning_assets(
                        user=request.user,
                        vector_asset_choice=asset_form.cleaned_data["vector_asset"],
                        insert_asset_choice=asset_form.cleaned_data["insert_asset"],
                    )
                except ValueError as exc:
                    asset_form.add_error(None, str(exc))
                else:
                    initial = dict(asset_form.cleaned_data)
                    initial["is_circular"] = "1" if assets.vector_asset.is_circular else "0"
                    initial["selected_enzymes"] = []
                    assembly_form = CloningConstructAssemblyForm(
                        initial=initial,
                        user=request.user,
                    )
        elif step in {"preview", "save"}:
            assembly_form = CloningConstructAssemblyForm(request.POST, user=request.user)
            if assembly_form.is_valid():
                try:
                    assets = resolve_cloning_assets(
                        user=request.user,
                        vector_asset_choice=assembly_form.cleaned_data["vector_asset"],
                        insert_asset_choice=assembly_form.cleaned_data["insert_asset"],
                    )
                except ValueError as exc:
                    assembly_form.add_error(None, str(exc))
                else:
                    left_enzyme = assembly_form.cleaned_data["left_enzyme"]
                    right_enzyme = assembly_form.cleaned_data["right_enzyme"]
                    if not left_enzyme or not right_enzyme:
                        if step == "save":
                            assembly_form.add_error(
                                None,
                                "Select both left and right enzymes before saving this construct.",
                            )
                    else:
                        try:
                            preview_data = preview_cloning_construct(
                                vector_asset=assets.vector_asset,
                                insert_asset=assets.insert_asset,
                                assembly_strategy=assembly_form.cleaned_data["assembly_strategy"],
                                left_enzyme=left_enzyme,
                                right_enzyme=right_enzyme,
                                is_circular=assembly_form.cleaned_data["is_circular"],
                                vector_fragment_index=assembly_form.cleaned_data.get("vector_fragment_index"),
                                insert_fragment_index=assembly_form.cleaned_data.get("insert_fragment_index"),
                            )
                        except ValueError as exc:
                            assembly_form.add_error(None, str(exc))
                        else:
                            if step == "save":
                                if not preview_data.is_valid:
                                    assembly_form.add_error(
                                        None,
                                        "Resolve the cloning validation messages before saving this construct.",
                                    )
                                else:
                                    construct = save_cloning_construct(
                                        name=assembly_form.cleaned_data["name"],
                                        description=assembly_form.cleaned_data["description"],
                                        preview_data=preview_data,
                                        user=request.user,
                                    )
                                    messages.success(request, "Construct assembled successfully.")
                                    return redirect("cloning_construct_detail", construct_id=construct.id)
            if assembly_form is not None:
                asset_form = CloningConstructAssetForm(
                    initial={
                        "name": assembly_form.data.get("name", ""),
                        "description": assembly_form.data.get("description", ""),
                        "assembly_strategy": assembly_form.data.get("assembly_strategy", ""),
                        "vector_asset": assembly_form.data.get("vector_asset", ""),
                        "insert_asset": assembly_form.data.get("insert_asset", ""),
                        "left_enzyme": assembly_form.data.get("left_enzyme", ""),
                        "right_enzyme": assembly_form.data.get("right_enzyme", ""),
                        "is_circular": assembly_form.data.get("is_circular", ""),
                        "selected_enzymes": assembly_form.data.getlist("selected_enzymes"),
                        "vector_fragment_index": assembly_form.data.get("vector_fragment_index", ""),
                        "insert_fragment_index": assembly_form.data.get("insert_fragment_index", ""),
                    },
                    user=request.user,
                )

    visual_preview = _build_assembly_visual_preview(assembly_form, assets)

    return render(
        request,
        "core/cloning_construct_form.html",
        {
            "asset_form": asset_form,
            "assembly_form": assembly_form,
            "preview_data": preview_data,
            "visual_preview": visual_preview,
            "assets": assets,
            "review_construct": review_construct,
        },
    )


@login_required
def cloning_construct_detail(request, construct_id):
    construct = _get_accessible_construct(request, construct_id)
    return render(
        request,
        "core/cloning_construct_detail.html",
        _build_construct_detail_context(request, construct),
    )


@login_required
def cloning_construct_linear_view(request, construct_id):
    construct = _get_accessible_construct(request, construct_id)
    try:
        context = _build_construct_linear_page_context(request, construct)
    except ValueError as exc:
        messages.error(request, str(exc))
        return redirect("cloning_construct_detail", construct_id=construct.id)

    return render(
        request,
        "core/cloning_construct_linear_view.html",
        context,
    )


@login_required
def cloning_construct_save_sequence_file(request, construct_id):
    if request.method != "POST":
        return HttpResponse("POST only", status=405)

    construct = _get_accessible_construct(request, construct_id)
    save_sequence_file_form = _build_construct_sequence_file_form(construct, request.POST)
    return_to = str(request.POST.get("return_to", "detail")).strip().lower()

    if save_sequence_file_form.is_valid():
        try:
            sequence_file = save_cloning_construct_sequence_file(
                construct=construct,
                user=request.user,
                name=save_sequence_file_form.cleaned_data["name"],
                description=save_sequence_file_form.cleaned_data["description"],
                file_type=save_sequence_file_form.cleaned_data["file_type"],
            )
        except ValueError as exc:
            save_sequence_file_form.add_error(None, str(exc))
        else:
            messages.success(request, f"Saved construct as sequence file {sequence_file.name}.")
            return redirect("sequencefile_linear_view", sequencefile_id=sequence_file.id)

    if return_to == "linear":
        try:
            context = _build_construct_linear_page_context(
                request,
                construct,
                save_sequence_file_form=save_sequence_file_form,
            )
        except ValueError as exc:
            messages.error(request, str(exc))
            context = _build_construct_detail_context(
                request,
                construct,
                save_sequence_file_form=save_sequence_file_form,
            )
            return render(
                request,
                "core/cloning_construct_detail.html",
                context,
                status=400,
            )
        return render(
            request,
            "core/cloning_construct_linear_view.html",
            context,
            status=400,
        )

    context = _build_construct_detail_context(
        request,
        construct,
        save_sequence_file_form=save_sequence_file_form,
    )
    return render(
        request,
        "core/cloning_construct_detail.html",
        context,
        status=400,
    )


@login_required
def cloning_construct_delete(request, construct_id):
    if request.method != "POST":
        return HttpResponse("POST only", status=405)

    construct = get_object_or_404(editable_cloning_constructs(request.user), id=construct_id)
    construct_name = construct.name
    construct.delete()
    messages.success(request, f"Deleted cloning construct {construct_name}.")
    return redirect("cloning_construct_list")


@login_required
def cloning_construct_download_genbank(request, construct_id):
    construct = _get_accessible_construct(request, construct_id)
    try:
        exported_genbank = export_cloning_construct_genbank(construct)
    except ValueError as exc:
        return HttpResponse(str(exc), status=400)

    response = HttpResponse(exported_genbank, content_type="application/genbank")
    filename = f"{str(construct.name).strip().replace(' ', '_') or 'construct'}.gb"
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response
