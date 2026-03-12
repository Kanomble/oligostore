import os
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.http import HttpResponse, HttpResponseForbidden
from django.shortcuts import render, redirect, get_object_or_404
from django.utils.text import slugify

from ..access import (
    accessible_pcr_products,
    accessible_primer_pairs,
    accessible_primers,
    accessible_projects,
    accessible_sequence_files,
)
from ..forms import ProjectForm
from ..models import Primer, PrimerPair, Project, SequenceFile
from ..services.user_assignment import assign_creator
from .utils import paginate_queryset


def _get_member_project(request, project_id):
    project = get_object_or_404(
        Project.objects.prefetch_related(
            "users",
            "primerpairs__forward_primer",
            "primerpairs__reverse_primer",
            "sequence_files",
            "pcr_products__sequence_file",
            "pcr_products__forward_primer",
            "pcr_products__reverse_primer",
        ),
        id=project_id,
    )
    if request.user not in project.users.all():
        return None
    return project


def _project_member_forbidden(message):
    return HttpResponseForbidden(message)


@login_required
def project_add_primerpair(request, project_id, pair_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to add these primer to your project.")
    pair = get_object_or_404(accessible_primer_pairs(request.user), id=pair_id)

    project.primerpairs.add(pair)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_remove_primerpair(request, project_id, pair_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to remove these primer pairs from your project.")
    pair = get_object_or_404(accessible_primer_pairs(request.user), id=pair_id)

    project.primerpairs.remove(pair)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_add_sequencefile(request, project_id, sequencefile_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to add sequence files to this project.")
    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    project.sequence_files.add(sequence_file)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_remove_sequencefile(request, project_id, sequencefile_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to remove sequence files from this project.")
    sequence_file = get_object_or_404(accessible_sequence_files(request.user), id=sequencefile_id)

    project.sequence_files.remove(sequence_file)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_add_pcr_product(request, project_id, pcr_product_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to add PCR products to this project.")
    pcr_product = get_object_or_404(accessible_pcr_products(request.user), id=pcr_product_id)

    project.pcr_products.add(pcr_product)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_remove_pcr_product(request, project_id, pcr_product_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You are not allowed to remove PCR products from this project.")
    pcr_product = get_object_or_404(accessible_pcr_products(request.user), id=pcr_product_id)

    project.pcr_products.remove(pcr_product)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_dashboard(request, project_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You do not have access to this project.")

    attached_pair_ids = set(project.primerpairs.values_list("id", flat=True))
    attached_sequence_file_ids = set(project.sequence_files.values_list("id", flat=True))
    attached_pcr_product_ids = set(project.pcr_products.values_list("id", flat=True))

    all_primers = accessible_primers(request.user)
    all_pairs = accessible_primer_pairs(request.user).select_related("forward_primer", "reverse_primer").order_by("name")
    all_sequence_files = accessible_sequence_files(request.user).order_by("name")
    all_pcr_products = accessible_pcr_products(request.user).select_related(
        "sequence_file",
        "forward_primer",
        "reverse_primer",
    ).order_by("-created_at", "name")

    return render(
        request,
        "core/project_dashboard.html",
        {
            "project": project,
            "all_primers": all_primers,
            "attached_pairs": project.primerpairs.all(),
            "available_pairs": [pair for pair in all_pairs if pair.id not in attached_pair_ids],
            "attached_sequence_files": project.sequence_files.all(),
            "available_sequence_files": [item for item in all_sequence_files if item.id not in attached_sequence_file_ids],
            "attached_pcr_products": project.pcr_products.all(),
            "available_pcr_products": [item for item in all_pcr_products if item.id not in attached_pcr_product_ids],
        },
    )


@login_required
def project_primer_list(request, project_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You do not have access to this project.")

    primers = accessible_primers(request.user).filter(
        Q(as_forward_primer__projects=project)
        | Q(as_reverse_primer__projects=project)
    ).distinct()

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
        "core/project_primer_list.html",
        {
            "project": project,
            "primers": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )


@login_required
def project_download_sequence_files(request, project_id):
    project = _get_member_project(request, project_id)
    if project is None:
        return _project_member_forbidden("You do not have access to this project.")

    sequence_files = project.sequence_files.all()
    if not sequence_files.exists():
        messages.error(request, "No sequence files are associated with this project.")
        return redirect("project_dashboard", project_id=project_id)

    archive_buffer = BytesIO()
    with ZipFile(archive_buffer, "w", ZIP_DEFLATED) as archive:
        for sequence_file in sequence_files:
            if not sequence_file.file:
                continue
            filename = os.path.basename(sequence_file.file.name)
            with sequence_file.file.open("rb") as file_handle:
                archive.writestr(filename, file_handle.read())

    archive_buffer.seek(0)
    filename = f"{slugify(project.name) or 'project'}-sequence-files.zip"
    response = HttpResponse(archive_buffer.getvalue(), content_type="application/zip")
    response["Content-Disposition"] = f'attachment; filename="{filename}"'
    return response


@login_required
def project_create(request):
    if request.method == "POST":
        form = ProjectForm(request.POST)
        if form.is_valid():
            project = form.save(commit=False)
            project = assign_creator(project, request.user)
            project.save()

            return redirect("project_list")
    else:
        form = ProjectForm()

    return render(request, "core/project_form.html", {"form": form})


@login_required(login_url="login")
def project_list(request):
    projects = accessible_projects(request.user)

    q = request.GET.get("q")
    if q:
        projects = projects.filter(
            Q(name__icontains=q)
            | Q(description__icontains=q)
        )

    order = request.GET.get("order", "created_desc")
    allowed_orders = {
        "name": "name",
        "name_desc": "-name",
        "created": "created_at",
        "created_desc": "-created_at",
    }
    projects = projects.order_by(allowed_orders.get(order, "-created_at"))
    page_obj, query_string = paginate_queryset(request, projects)
    return render(
        request,
        "core/project_list.html",
        {
            "projects": page_obj,
            "page_obj": page_obj,
            "query_string": query_string,
        },
    )
