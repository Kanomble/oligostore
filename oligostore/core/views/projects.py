import os
from io import BytesIO
from zipfile import ZipFile, ZIP_DEFLATED

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.http import HttpResponse, HttpResponseForbidden
from django.shortcuts import render, redirect, get_object_or_404
from django.utils.text import slugify

from ..forms import ProjectForm
from ..models import Primer, PrimerPair, Project, SequenceFile
from ..services.user_assignment import assign_creator
from .utils import paginate_queryset


@login_required
def project_add_primerpair(request, project_id, pair_id):
    project = get_object_or_404(Project, id=project_id, users=request.user)
    pair = get_object_or_404(PrimerPair, id=pair_id, users=request.user)

    if request.user not in project.users.all():
        return HttpResponseForbidden("You are not allowed to add these primer to your project.")

    project.primerpairs.add(pair)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_remove_primerpair(request, project_id, pair_id):
    project = get_object_or_404(Project, id=project_id, users=request.user)
    pair = get_object_or_404(PrimerPair, id=pair_id, users=request.user)

    if request.user not in project.users.all():
        return HttpResponseForbidden("You are not allowed to remove these primer pairs from your project.")

    project.primerpairs.remove(pair)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_add_sequencefile(request, project_id, sequencefile_id):
    project = get_object_or_404(Project, id=project_id, users=request.user)
    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequencefile_id,
        uploaded_by=request.user,
    )

    if request.user not in project.users.all():
        return HttpResponseForbidden("You are not allowed to add sequence files to this project.")

    project.sequence_files.add(sequence_file)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_remove_sequencefile(request, project_id, sequencefile_id):
    project = get_object_or_404(Project, id=project_id, users=request.user)
    sequence_file = get_object_or_404(
        SequenceFile,
        id=sequencefile_id,
        uploaded_by=request.user,
    )

    if request.user not in project.users.all():
        return HttpResponseForbidden("You are not allowed to remove sequence files from this project.")

    project.sequence_files.remove(sequence_file)
    return redirect("project_dashboard", project_id=project_id)


@login_required
def project_dashboard(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.user not in project.users.all():
        return HttpResponseForbidden("You do not have access to this project.")

    all_primers = Primer.objects.filter(users=request.user)
    all_pairs = PrimerPair.objects.filter(users=request.user)
    all_sequence_files = SequenceFile.objects.filter(uploaded_by=request.user)

    return render(
        request,
        "core/project_dashboard.html",
        {
            "project": project,
            "all_primers": all_primers,
            "all_pairs": all_pairs,
            "all_sequence_files": all_sequence_files,
        },
    )


@login_required
def project_primer_list(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    if request.user not in project.users.all():
        return HttpResponseForbidden("You do not have access to this project.")

    primers = Primer.objects.filter(users=request.user).filter(
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
    project = get_object_or_404(Project, id=project_id)

    if request.user not in project.users.all():
        return HttpResponseForbidden("You do not have access to this project.")

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
    projects = request.user.project_access.all()

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