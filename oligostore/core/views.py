from django.shortcuts import render, redirect, get_object_or_404
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse, HttpResponse, HttpResponseForbidden
from django.contrib import messages
from django.core.exceptions import PermissionDenied
from django.db.models import Q
from .models import Primer, PrimerPair, Project, SequenceFile
from .forms import PrimerForm, CustomUserCreationForm, \
    PrimerPairForm, PrimerPairCreateCombinedForm, \
    ProjectForm, Primer3GlobalArgsForm
from .services.primer_analysis import analyze_primer, analyze_cross_dimer, \
    analyze_sequence, sanitize_sequence, reverse_complement, find_binding_site,\
    window_sequence, render_windowed_line, highlight_binding
from .services.primer_binding import analyze_primer_binding
from .services.user_assignment import assign_creator
import re

@login_required
def sequencefile_upload(request):
    """
    Upload a FASTA or GenBank sequence file and persist it.
    """

    if request.method == "POST":
        name = request.POST.get("name")
        file = request.FILES.get("file")
        file_type = request.POST.get("file_type")
        description = request.POST.get("description", "")

        if not name or not file or file_type not in ("fasta", "genbank"):
            return render(
                request,
                "core/sequencefile_upload.html",
                {
                    "error": "All required fields must be provided.",
                },
            )

        SequenceFile.objects.create(
            name=name,
            file=file,
            file_type=file_type,
            description=description,
            uploaded_by=request.user,
        )

        return redirect("sequencefile_list")

    return render(
        request,
        "core/sequencefile_upload.html",
    )

@login_required
def sequencefile_list(request):
    """
    List uploaded FASTA / GenBank sequence files for the current user.
    Supports optional filtering by search term and file type.
    """

    qs = SequenceFile.objects.filter(uploaded_by=request.user)

    # Search filter
    q = request.GET.get("q")
    if q:
        qs = qs.filter(
            Q(name__icontains=q) |
            Q(description__icontains=q)
        )

    # File type filter
    file_type = request.GET.get("type")
    if file_type in ("fasta", "genbank"):
        qs = qs.filter(file_type=file_type)

    qs = qs.order_by("-uploaded_at")

    return render(
        request,
        "core/sequencefile_list.html",
        {
            "sequence_files": qs,
        },
    )

@login_required
def primer_binding_analysis(request):
    primers = Primer.objects.filter(users=request.user)
    sequence_files = SequenceFile.objects.filter(uploaded_by=request.user)

    preselected_primer = request.GET.get("primer")
    preselected_sequence_file = request.GET.get("sequence_file")

    if preselected_primer:
        preselected_primer = (
            primers.filter(id=preselected_primer)
            .values_list("id", flat=True)
            .first()
        )

    if preselected_sequence_file:
        preselected_sequence_file = (
            sequence_files.filter(id=preselected_sequence_file)
            .values_list("id", flat=True)
            .first()
        )

    if request.method == "POST":
        primer_id = request.POST.get("primer_id")
        sequence_file_id = request.POST.get("sequence_file_id")

        primer = get_object_or_404(
            Primer,
            id=primer_id,
            users=request.user,
        )

        sequence_file = get_object_or_404(
            SequenceFile,
            id=sequence_file_id,
            uploaded_by=request.user,
        )

        hits = analyze_primer_binding(
            primer_sequence=primer.sequence,
            sequence_file=sequence_file,
            max_mismatches=2,
        )

        return render(
            request,
            "core/primer_binding_results.html",
            {
                "primer": primer,
                "sequence_file": sequence_file,
                "hits": hits,
            },
        )

    return render(
        request,
        "core/primer_binding_upload.html",
        {
            "primers": primers,
            "sequence_files": sequence_files,
            "preselected_primer": preselected_primer,
            "preselected_sequence_file": preselected_sequence_file,
        },
    )

@login_required
def download_product_sequence(request):
    product_seq = request.POST.get("product_sequence")
    pair_index = request.POST.get("pair_index", "unknown")

    if not product_seq:
        return HttpResponse("No product sequence provided.", status=400)

    fasta = (
        f">PCR_product_pair_{pair_index}\n"
        f"{product_seq}\n"
    )

    response = HttpResponse(fasta, content_type="text/plain")
    response["Content-Disposition"] = (
        f"attachment; filename=pcr_product_pair_{pair_index}.fasta"
    )
    return response

@login_required
def save_generated_primerpair(request):
    if request.method != "POST":
        return redirect("analyze_sequence")

    try:
        # Input fields
        left_seq = request.POST.get("left_seq")
        right_seq = request.POST.get("right_seq")
        left_name = request.POST.get("forward_name")
        right_name = request.POST.get("reverse_name")
        pair_name = request.POST.get("pair_name")

        # Cross-dimer score
        # currently not used
        # hetero_dimer = analyze_cross_dimer(left_seq, right_seq)

        # ---------------------------
        # Create FORWARD primer
        # ---------------------------
        forward = Primer(
            primer_name=left_name,
            sequence=left_seq,
        )
        forward = assign_creator(forward, request.user)
        analysis = analyze_primer(forward.sequence)
        forward.length = len(forward.sequence)
        forward.gc_content = analysis["gc_content"]
        forward.tm = analysis["tm"]
        forward.hairpin_dg = analysis["hairpin_dg"]
        forward.self_dimer_dg = analysis["self_dimer_dg"]

        forward.save()
        # ---------------------------
        # Create REVERSE primer
        # ---------------------------
        reverse = Primer(
            primer_name=right_name,
            sequence=right_seq,
        )
        reverse = assign_creator(reverse, request.user)
        analysis = analyze_primer(reverse.sequence)
        reverse.length = len(reverse.sequence)
        reverse.gc_content = analysis["gc_content"]
        reverse.tm = analysis["tm"]
        reverse.hairpin_dg = analysis["hairpin_dg"]
        reverse.self_dimer_dg = analysis["self_dimer_dg"]

        reverse.save()

        # ---------------------------
        # Create PRIMER PAIR
        # ---------------------------
        pair = PrimerPair.objects.create(
            name=pair_name,
            forward_primer=forward,
            reverse_primer=reverse,
        )
        pair = assign_creator(pair, request.user)
        pair.save()
        messages.success(request, f"Primer Pair '{pair.name}' saved successfully!")
        return redirect("primerpair_list")

    except Exception as e:
        messages.error(request, f"ERROR saving primer pair: {e}")
        return redirect("analyze_sequence")

@login_required
def analyze_sequence_view(request):
    if request.method == "POST":
        form = Primer3GlobalArgsForm(request.POST)
        raw_seq = request.POST.get("sequence")

        if form.is_valid():
            try:
                sequence = sanitize_sequence(raw_seq)
            except ValueError as e:
                return render(
                    request,
                    "core/analyze_sequence.html",
                    {
                        "form": form,
                        "error_message": str(e),
                        "product_size_fields": ["PRIMER_PRODUCT_SIZE_RANGE"],
                        "primer_length_fields": ["PRIMER_OPT_SIZE", "PRIMER_MIN_SIZE", "PRIMER_MAX_SIZE"],
                        "tm_fields": ["PRIMER_OPT_TM", "PRIMER_MIN_TM", "PRIMER_MAX_TM"],
                        "gc_fields": ["PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_OPT_GC_PERCENT"],
                        "self_comp_fields": ["PRIMER_MAX_SELF_ANY", "PRIMER_MAX_SELF_END", "PRIMER_MAX_HAIRPIN_TH"],
                        "product_tm_fields": [
                            "PRIMER_PRODUCT_OPT_TM",
                            "PRIMER_PRODUCT_MIN_TM",
                            "PRIMER_PRODUCT_MAX_TM",
                        ],
                        "chemistry_fields": [
                            "PRIMER_SALT_MONOVALENT",
                            "PRIMER_SALT_DIVALENT",
                            "PRIMER_DNTP_CONC",
                            "PRIMER_DNA_CONC",
                        ],
                        "misc_fields": [
                            "PRIMER_MAX_POLY_X",
                            "PRIMER_EXPLAIN_FLAG",
                            "PRIMER_NUM_RETURN",
                            "PRIMER_GC_CLAMP",
                        ],
                    },
                )

            # ---------------- GLOBAL ARGS ----------------

            global_args = form.cleaned_data.copy()

            primer_sides = global_args.pop("PRIMER_SIDES", ["LEFT", "RIGHT"])
            if not primer_sides:
                return render(
                    request,
                    "core/analyze_sequence.html",
                    {
                        "form": form,
                        "error_message": "Select at least one primer direction.",
                    },
                )

            global_args["PRIMER_PICK_LEFT_PRIMER"] = 1 if "LEFT" in primer_sides else 0
            global_args["PRIMER_PICK_RIGHT_PRIMER"] = 1 if "RIGHT" in primer_sides else 0

            # Normalize boolean flags to Primer3 ints
            for key in (
                "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT",
                "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT",
            ):
                global_args[key] = 1 if global_args.get(key) else 0

            # Fill defaults for optional fields
            for field_name, field in form.fields.items():
                if field_name not in global_args or global_args[field_name] is None:
                    global_args[field_name] = field.initial

            # Remove unset product TM constraints
            for key in (
                "PRIMER_PRODUCT_MIN_TM",
                "PRIMER_PRODUCT_MAX_TM",
                "PRIMER_PRODUCT_OPT_TM",
            ):
                if not global_args.get(key):
                    global_args.pop(key, None)

            # ---------------- RUN PRIMER3 ----------------

            primer_list, raw_results, mode = analyze_sequence(sequence, global_args)

            # ---------------- VISUALIZATION ENRICHMENT ----------------
            FLANK_SIZE = 25

            for p in primer_list:
                if p["mode"] == "PAIR":
                    fwd = p["left_seq"].upper()
                    rev = p["right_seq"].upper()
                    rev_rc = reverse_complement(rev)

                    fwd_pos = find_binding_site(sequence, fwd)
                    rev_pos = find_binding_site(sequence, rev_rc)

                    fwd_window, fwd_start, fwd_len = window_sequence(
                        sequence, fwd_pos, len(fwd), flank=FLANK_SIZE
                    )
                    rev_window, rev_start, rev_len = window_sequence(
                        sequence, rev_pos, len(rev_rc), flank=FLANK_SIZE
                    )

                    p["forward_window"] = highlight_binding(
                        fwd_window, fwd_start, fwd_len
                    )
                    p["reverse_window"] = highlight_binding(
                        rev_window, rev_start, rev_len
                    )

                    p["forward_window_line"] = render_windowed_line(
                        fwd_window, fwd_start, fwd_len
                    )
                    p["reverse_window_line"] = render_windowed_line(
                        rev_window, rev_start, rev_len
                    )
                    p["forward_pos"] = fwd_pos
                    p["reverse_pos"] = rev_pos

                    # for downloading product
                    product_start = fwd_pos
                    product_end = rev_pos + len(rev)

                    p["product_sequence"] = sequence[product_start:product_end]
                    p["product_length"] = len(p["product_sequence"])

                else:
                    primer = p["seq"].upper()
                    if p["mode"] == "RIGHT":
                        primer = reverse_complement(primer)

                    pos = find_binding_site(sequence, primer)

                    window, start, length = window_sequence(
                        sequence, pos, len(primer), flank=FLANK_SIZE
                    )

                    p["window"] = highlight_binding(
                        window, start, length
                    )

                    p["window_line"] = render_windowed_line(
                        window, start, length
                    )
                    p["pos"] = pos

            # ---------------- RENDER RESULTS ----------------

            return render(
                request,
                "core/analyze_sequence_results.html",
                {
                    "sequence": sequence,
                    "primer_list": primer_list,
                    "raw_results": raw_results,
                    "mode": mode,
                },
            )

    # ---------------- GET REQUEST ----------------

    form = Primer3GlobalArgsForm()
    return render(
        request,
        "core/analyze_sequence.html",
        {
            "form": form,
            "product_size_fields": ["PRIMER_PRODUCT_SIZE_RANGE"],
            "primer_length_fields": ["PRIMER_OPT_SIZE", "PRIMER_MIN_SIZE", "PRIMER_MAX_SIZE"],
            "tm_fields": ["PRIMER_OPT_TM", "PRIMER_MIN_TM", "PRIMER_MAX_TM"],
            "gc_fields": ["PRIMER_MIN_GC", "PRIMER_MAX_GC", "PRIMER_OPT_GC_PERCENT"],
            "self_comp_fields": ["PRIMER_MAX_SELF_ANY", "PRIMER_MAX_SELF_END", "PRIMER_MAX_HAIRPIN_TH"],
            "product_tm_fields": [
                "PRIMER_PRODUCT_OPT_TM",
                "PRIMER_PRODUCT_MIN_TM",
                "PRIMER_PRODUCT_MAX_TM",
            ],
            "chemistry_fields": [
                "PRIMER_SALT_MONOVALENT",
                "PRIMER_SALT_DIVALENT",
                "PRIMER_DNTP_CONC",
                "PRIMER_DNA_CONC",
            ],
            "misc_fields": [
                "PRIMER_MAX_POLY_X",
                "PRIMER_EXPLAIN_FLAG",
                "PRIMER_NUM_RETURN",
                "PRIMER_GC_CLAMP",
            ],
        },
    )

@login_required
def analyze_primer_view(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST only"}, status=400)

    seq = request.POST.get("sequence", "")
    if not seq:
        return JsonResponse({"error":"Sequence is empty"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", seq):
        return JsonResponse({"error": "Sequence contains invalid characters (allowed: A C G T N)"}, status=400)

    primer_analysis_dict = analyze_primer(seq)

    return JsonResponse(primer_analysis_dict)

def register(request):
    if request.method == "POST":
        form = CustomUserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            return redirect("login")
    else:
        form = CustomUserCreationForm()

    return render(request, "registration/register.html", {"form": form})

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
def project_dashboard(request, project_id):
    project = get_object_or_404(Project, id=project_id)

    # Ensure the user has access
    if request.user not in project.users.all():
        return HttpResponseForbidden("You do not have access to this project.")

    # All objects user can link
    all_primers = Primer.objects.filter(users=request.user)
    all_pairs = PrimerPair.objects.filter(users=request.user)

    return render(request, "core/project_dashboard.html", {
        "project": project,
        "all_primers": all_primers,
        "all_pairs": all_pairs,
    })

@login_required
def project_create(request):
    if request.method == "POST":
        form = ProjectForm(request.POST)
        if form.is_valid():
            project = form.save(commit=False)
            project = assign_creator(project, request.user)
            project.save()

            # Creator automatically has access
            return redirect("project_list")
    else:
        form = ProjectForm()

    return render(request, "core/project_form.html", {"form": form})


@login_required(login_url="login")
def project_list(request):
    projects = request.user.project_access.all()
    return render(request, "core/project_list.html",{"projects":projects})

@login_required
def primerpair_delete(request, primerpair_id):
    pair = get_object_or_404(PrimerPair, id=primerpair_id)
    if pair.creator != request.user:
        raise PermissionDenied("You are not the creator of this primer.")
    pair.delete()
    return redirect("primerpair_list")

@login_required
def analyze_primerpair_view(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST required"}, status=400)

    fwd = request.POST.get("forward_sequence", "").strip()
    rev = request.POST.get("reverse_sequence", "").strip()

    # Validate input
    if not fwd or not rev:
        return JsonResponse({"error": "Both sequences are required"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", fwd):
        return JsonResponse({"error": "Forward sequence contains invalid characters"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", rev):
        return JsonResponse({"error": "Reverse sequence contains invalid characters"}, status=400)

    try:
        # Analyze forward primer
        forward_analysis = analyze_primer(fwd)

        # Analyze reverse primer
        reverse_analysis = analyze_primer(rev)

        hetero_dimer = analyze_cross_dimer(fwd, rev)

        return JsonResponse({
            "forward": {
                "tm": forward_analysis["tm"],
                "gc": forward_analysis["gc_content"],
                "hairpin_dg": forward_analysis["hairpin_dg"],
                "self_dimer_dg": forward_analysis["self_dimer_dg"],
            },
            "reverse": {
                "tm": reverse_analysis["tm"],
                "gc": reverse_analysis["gc_content"],
                "hairpin_dg": reverse_analysis["hairpin_dg"],
                "self_dimer_dg": reverse_analysis["self_dimer_dg"],
            },
            "hetero_dimer_dg": hetero_dimer.dg,
        })

    except Exception as e:
        return JsonResponse({"error": f"Analysis failed: {str(e)}"}, status=500)


@login_required
def primerpair_combined_create(request):
    if request.method == "POST":
        form = PrimerPairCreateCombinedForm(request.POST)
        if form.is_valid():

            # 1) Create forward primer
            forward = Primer(
                primer_name=form.cleaned_data["forward_name"],
                sequence=form.cleaned_data["forward_sequence"],
            )
            forward = assign_creator(forward, request.user)
            analysis = analyze_primer(forward.sequence)
            forward.length = len(forward.sequence)
            # run analysis
            forward.gc_content = analysis["gc_content"]
            forward.tm = analysis["tm"]
            forward.hairpin_dg = analysis["hairpin_dg"]
            forward.self_dimer_dg = analysis["self_dimer_dg"]
            forward.save()

            # 2) Create reverse primer
            reverse = Primer(
                primer_name=form.cleaned_data["reverse_name"],
                sequence=form.cleaned_data["reverse_sequence"],
            )
            reverse = assign_creator(reverse, request.user)
            analysis = analyze_primer(reverse.sequence)
            reverse.length = len(reverse.sequence)
            # run analysis
            reverse.gc_content = analysis["gc_content"]
            reverse.tm = analysis["tm"]
            reverse.hairpin_dg = analysis["hairpin_dg"]
            reverse.self_dimer_dg = analysis["self_dimer_dg"]

            reverse.save()

            # 3) Create primer pair
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

@login_required(login_url="login")
def primerpair_list(request):
    primer_pairs = PrimerPair.objects.filter(users=request.user)
    return render(request, "core/primerpair_list.html",{"primer_pairs":primer_pairs})

@login_required(login_url="login")
def primer_list(request):
    primers = Primer.objects.filter(users=request.user).order_by('-created_at')
    return render(request, "core/primer_list.html", {"primers": primers})

@login_required(login_url="login")
def primer_create(request):
    if request.method == "POST":
        form = PrimerForm(request.POST)
        if form.is_valid():
            primer = form.save(commit=False)
            primer = assign_creator(primer, user=request.user)
            primer.length = len(primer.sequence)
            # run analysis
            analysis = analyze_primer(primer.sequence)
            # set analysis variables to primer model
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

def home(request):
    return render(request, "core/home.html")
