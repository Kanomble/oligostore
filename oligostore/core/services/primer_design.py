from .sequence_utils import reverse_complement
from .primer_analysis import (
    find_binding_site,
    highlight_binding,
    render_windowed_line,
    window_sequence,
)


def enrich_primer_design_results(sequence, primer_list, flank_size=25):
    for primer_result in primer_list:
        if primer_result["mode"] == "PAIR":
            _enrich_pair_result(sequence, primer_result, flank_size)
        else:
            _enrich_single_result(sequence, primer_result, flank_size)
    return primer_list


def _enrich_pair_result(sequence, primer_result, flank_size):
    forward_sequence = primer_result["left_seq"].upper()
    reverse_sequence = primer_result["right_seq"].upper()
    reverse_complement_sequence = reverse_complement(reverse_sequence)

    forward_pos = find_binding_site(sequence, forward_sequence)
    reverse_pos = find_binding_site(sequence, reverse_complement_sequence)

    forward_window, forward_start, forward_len = window_sequence(
        sequence,
        forward_pos,
        len(forward_sequence),
        flank=flank_size,
    )
    reverse_window, reverse_start, reverse_len = window_sequence(
        sequence,
        reverse_pos,
        len(reverse_complement_sequence),
        flank=flank_size,
    )

    primer_result["forward_window"] = highlight_binding(
        forward_window,
        forward_start,
        forward_len,
    )
    primer_result["reverse_window"] = highlight_binding(
        reverse_window,
        reverse_start,
        reverse_len,
    )
    primer_result["forward_window_line"] = render_windowed_line(
        forward_window,
        forward_start,
        forward_len,
    )
    primer_result["reverse_window_line"] = render_windowed_line(
        reverse_window,
        reverse_start,
        reverse_len,
    )
    primer_result["forward_pos"] = forward_pos
    primer_result["reverse_pos"] = reverse_pos

    product_start = forward_pos
    product_end = reverse_pos + len(reverse_sequence)
    primer_result["product_sequence"] = sequence[product_start:product_end]
    primer_result["product_length"] = len(primer_result["product_sequence"])


def _enrich_single_result(sequence, primer_result, flank_size):
    primer_sequence = primer_result["seq"].upper()
    if primer_result["mode"] == "RIGHT":
        primer_sequence = reverse_complement(primer_sequence)

    pos = find_binding_site(sequence, primer_sequence)
    window, start, length = window_sequence(
        sequence,
        pos,
        len(primer_sequence),
        flank=flank_size,
    )
    primer_result["window"] = highlight_binding(window, start, length)
    primer_result["window_line"] = render_windowed_line(window, start, length)
    primer_result["pos"] = pos
