from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect, render

from ..access import (
    accessible_cloning_constructs,
    accessible_pcr_products,
    accessible_sequence_files,
)
from ..forms import CloningConstructForm
from ..services.cloning import create_cloning_construct, resolve_asset_choice
from ..services.listing import apply_ordering, apply_search
from .utils import paginate_queryset


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
    if request.method == "POST":
        form = CloningConstructForm(request.POST, user=request.user)
        if form.is_valid():
            try:
                vector_asset = resolve_asset_choice(
                    user=request.user,
                    choice=form.cleaned_data["vector_asset"],
                )
                insert_asset = resolve_asset_choice(
                    user=request.user,
                    choice=form.cleaned_data["insert_asset"],
                )
            except ValueError as exc:
                form.add_error(None, str(exc))
            else:
                try:
                    construct = create_cloning_construct(
                        name=form.cleaned_data["name"],
                        description=form.cleaned_data["description"],
                        vector_asset=vector_asset,
                        insert_asset=insert_asset,
                        assembly_strategy=form.cleaned_data["assembly_strategy"],
                        left_enzyme=form.cleaned_data["left_enzyme"],
                        right_enzyme=form.cleaned_data["right_enzyme"],
                        user=request.user,
                    )
                except ValueError as exc:
                    form.add_error(None, str(exc))
                else:
                    if construct.is_valid:
                        messages.success(request, "Construct assembled successfully.")
                    else:
                        messages.warning(request, "Construct saved with validation warnings.")
                    return redirect("cloning_construct_detail", construct_id=construct.id)
    else:
        form = CloningConstructForm(user=request.user)
    return render(request, "core/cloning_construct_form.html", {"form": form})


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
    return render(request, "core/cloning_construct_detail.html", {"construct": construct})
