import re

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.http import JsonResponse, HttpResponse
from django.shortcuts import render, redirect

from ..forms import Primer3GlobalArgsForm
from ..models import Primer, PrimerPair
from ..services.primer_analysis import (
    analyze_primer,
    analyze_cross_dimer,
    analyze_sequence,
    sanitize_sequence,
    reverse_complement,
    find_binding_site,
    window_sequence,
    render_windowed_line,
    highlight_binding,
)
from ..services.user_assignment import assign_creator


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

            for key in (
                "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT",
                "PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT",
            ):
                global_args[key] = 1 if global_args.get(key) else 0

            for field_name, field in form.fields.items():
                if field_name not in global_args or global_args[field_name] is None:
                    global_args[field_name] = field.initial

            for key in (
                "PRIMER_PRODUCT_MIN_TM",
                "PRIMER_PRODUCT_MAX_TM",
                "PRIMER_PRODUCT_OPT_TM",
            ):
                if not global_args.get(key):
                    global_args.pop(key, None)

            primer_list, raw_results, mode = analyze_sequence(sequence, global_args)

            flank_size = 25

            for p in primer_list:
                if p["mode"] == "PAIR":
                    fwd = p["left_seq"].upper()
                    rev = p["right_seq"].upper()
                    rev_rc = reverse_complement(rev)

                    fwd_pos = find_binding_site(sequence, fwd)
                    rev_pos = find_binding_site(sequence, rev_rc)

                    fwd_window, fwd_start, fwd_len = window_sequence(
                        sequence, fwd_pos, len(fwd), flank=flank_size
                    )
                    rev_window, rev_start, rev_len = window_sequence(
                        sequence, rev_pos, len(rev_rc), flank=flank_size
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
                        sequence, pos, len(primer), flank=flank_size
                    )

                    p["window"] = highlight_binding(
                        window, start, length
                    )

                    p["window_line"] = render_windowed_line(
                        window, start, length
                    )
                    p["pos"] = pos

            request.session["last_primer_results"] = {
                "sequence": sequence,
                "primer_list": primer_list,
                "mode": mode,
            }
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
        return JsonResponse({"error": "Sequence is empty"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", seq):
        return JsonResponse({"error": "Sequence contains invalid characters (allowed: A C G T N)"}, status=400)

    primer_analysis_dict = analyze_primer(seq)

    return JsonResponse(primer_analysis_dict)


@login_required
def analyze_primerpair_view(request):
    if request.method != "POST":
        return JsonResponse({"error": "POST required"}, status=400)

    fwd = request.POST.get("forward_sequence", "").strip()
    rev = request.POST.get("reverse_sequence", "").strip()

    if not fwd or not rev:
        return JsonResponse({"error": "Both sequences are required"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", fwd):
        return JsonResponse({"error": "Forward sequence contains invalid characters"}, status=400)

    if not re.fullmatch(r"[ACGTNacgtn]+", rev):
        return JsonResponse({"error": "Reverse sequence contains invalid characters"}, status=400)

    try:
        forward_analysis = analyze_primer(fwd)

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
        left_seq = request.POST.get("left_seq")
        right_seq = request.POST.get("right_seq")
        left_name = request.POST.get("forward_name")
        right_name = request.POST.get("reverse_name")
        pair_name = request.POST.get("pair_name")

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

        pair = PrimerPair(
            name=pair_name,
            forward_primer=forward,
            reverse_primer=reverse,
        )
        pair = assign_creator(pair, request.user)
        pair.save()
        messages.success(request, f"Primer Pair '{pair.name}' saved successfully!")
        last_results = request.session.get("last_primer_results")
        if last_results:
            return render(
                request,
                "core/analyze_sequence_results.html",
                last_results,
            )
        return redirect("analyze_sequence")

    except Exception as e:
        messages.error(request, f"ERROR saving primer pair: {e}")
        return redirect("analyze_sequence")