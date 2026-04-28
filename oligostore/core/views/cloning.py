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
from ..forms import CloningConstructAssetForm, CloningConstructAssemblyForm
from ..services.cloning import (
    build_pcr_product_asset_choice,
    build_cloning_construct_detail_display,
    build_sequence_file_asset_choice,
    preview_cloning_construct,
    resolve_cloning_assets,
    save_cloning_construct,
)
from ..services.cloning_exports import export_cloning_construct_genbank
from ..services.listing import apply_ordering, apply_search
from .utils import paginate_queryset


def _build_construct_asset_form_values(construct):
    if construct.vector_source_type == construct.SOURCE_SEQUENCE_FILE:
        vector_asset = build_sequence_file_asset_choice(
            sequence_file_id=construct.vector_sequence_file_id,
            record_id=construct.vector_record_id,
        ).encoded_value
    else:
        vector_asset = build_pcr_product_asset_choice(
            pcr_product_id=construct.vector_pcr_product_id,
        ).encoded_value

    if construct.insert_source_type == construct.SOURCE_SEQUENCE_FILE:
        insert_asset = build_sequence_file_asset_choice(
            sequence_file_id=construct.insert_sequence_file_id,
            record_id=construct.insert_record_id,
        ).encoded_value
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
    }


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
            "vector_pcr_product__name",
            "insert_sequence_file__name",
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
                    assembly_form = CloningConstructAssemblyForm(
                        initial=asset_form.cleaned_data,
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
                    try:
                        preview_data = preview_cloning_construct(
                            vector_asset=assets.vector_asset,
                            insert_asset=assets.insert_asset,
                            assembly_strategy=assembly_form.cleaned_data["assembly_strategy"],
                            left_enzyme=assembly_form.cleaned_data["left_enzyme"],
                            right_enzyme=assembly_form.cleaned_data["right_enzyme"],
                        )
                    except ValueError as exc:
                        assembly_form.add_error(None, str(exc))
                    else:
                        if step == "save":
                            construct = save_cloning_construct(
                                name=assembly_form.cleaned_data["name"],
                                description=assembly_form.cleaned_data["description"],
                                preview_data=preview_data,
                                user=request.user,
                            )
                            if construct.is_valid:
                                messages.success(request, "Construct assembled successfully.")
                            else:
                                messages.warning(request, "Construct saved with validation warnings.")
                            return redirect("cloning_construct_detail", construct_id=construct.id)
            if assembly_form is not None:
                asset_form = CloningConstructAssetForm(
                    initial={
                        "name": assembly_form.data.get("name", ""),
                        "description": assembly_form.data.get("description", ""),
                        "assembly_strategy": assembly_form.data.get("assembly_strategy", ""),
                        "vector_asset": assembly_form.data.get("vector_asset", ""),
                        "insert_asset": assembly_form.data.get("insert_asset", ""),
                    },
                    user=request.user,
                )

    return render(
        request,
        "core/cloning_construct_form.html",
        {
            "asset_form": asset_form,
            "assembly_form": assembly_form,
            "preview_data": preview_data,
            "assets": assets,
            "review_construct": review_construct,
        },
    )


@login_required
def cloning_construct_detail(request, construct_id):
    construct = get_object_or_404(
        accessible_cloning_constructs(request.user).select_related(
            "vector_sequence_file",
            "vector_pcr_product",
            "insert_sequence_file",
            "insert_pcr_product",
        ),
        id=construct_id,
    )
    detail_display = build_cloning_construct_detail_display(construct)
    return render(
        request,
        "core/cloning_construct_detail.html",
        {
            "construct": construct,
            "detail_display": detail_display,
            "can_delete": construct.creator_id == request.user.id,
        },
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
    construct = get_object_or_404(
        accessible_cloning_constructs(request.user).select_related(
            "vector_sequence_file",
            "vector_pcr_product",
            "insert_sequence_file",
            "insert_pcr_product",
            "insert_pcr_product__sequence_file",
            "vector_pcr_product__sequence_file",
        ),
        id=construct_id,
    )
    try:
        exported_genbank = export_cloning_construct_genbank(construct)
    except ValueError as exc:
        return HttpResponse(str(exc), status=400)

    response = HttpResponse(exported_genbank, content_type="application/genbank")
    filename = f"{str(construct.name).strip().replace(' ', '_') or 'construct'}.gb"
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response
