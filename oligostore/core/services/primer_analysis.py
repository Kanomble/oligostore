import primer3
import re
from core.services.sequence_utils import reverse_complement

def analyze_primer(seq: str) -> dict:
    seq = seq.upper()

    # primer3 objects
    hairpin = primer3.calc_hairpin(seq)
    dimer = primer3.calc_homodimer(seq)

    gc_content = (seq.count("G") + seq.count("C"))/len(seq)

    return {
        "gc_content": round(gc_content, 2),
        "tm": round(primer3.calcTm(seq),2),
        "hairpin_dg": round(hairpin.dg / 1000, 2),
        "self_dimer_dg": round(dimer.dg / 1000, 2),
        "hairpin":hairpin.structure_found,
        "self_dimer":dimer.structure_found
    }

def analyze_cross_dimer(fwd, rev):
    dimer = primer3.calcHeterodimer(fwd, rev)
    return dimer

def analyze_sequence(sequence: str, global_args: dict):
    primer3_input = {
        "SEQUENCE_ID": "user_sequence",
        "SEQUENCE_TEMPLATE": sequence,
    }

    results = primer3.bindings.design_primers(primer3_input, global_args)

    primer_list = []

    left_n = results.get("PRIMER_LEFT_NUM_RETURNED", 0)
    right_n = results.get("PRIMER_RIGHT_NUM_RETURNED", 0)
    pair_n = results.get("PRIMER_PAIR_NUM_RETURNED", 0)

    if pair_n > 0:
        mode = "PAIR"
        for i in range(pair_n):
            primer_list.append({
                "mode": "PAIR",
                "left_seq": results.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
                "right_seq": results.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
                "left_tm": results.get(f"PRIMER_LEFT_{i}_TM"),
                "right_tm": results.get(f"PRIMER_RIGHT_{i}_TM"),
                "product_size": results.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE"),
                "penalty": results.get(f"PRIMER_PAIR_{i}_PENALTY"),
            })

    elif left_n > 0:
        mode = "LEFT"
        for i in range(left_n):
            primer_list.append({
                "mode": "LEFT",
                "seq": results.get(f"PRIMER_LEFT_{i}_SEQUENCE"),
                "tm": results.get(f"PRIMER_LEFT_{i}_TM"),
                "penalty": results.get(f"PRIMER_LEFT_{i}_PENALTY"),
            })

    elif right_n > 0:
        mode = "RIGHT"
        for i in range(right_n):
            primer_list.append({
                "mode": "RIGHT",
                "seq": results.get(f"PRIMER_RIGHT_{i}_SEQUENCE"),
                "tm": results.get(f"PRIMER_RIGHT_{i}_TM"),
                "penalty": results.get(f"PRIMER_RIGHT_{i}_PENALTY"),
            })

    else:
        mode = "NONE"

    return primer_list, results, mode

def sanitize_sequence(raw: str) -> str:
    """
    Cleans user input so Primer3 receives a valid DNA sequence.
    - Removes whitespace, tabs, newlines
    - Converts to uppercase
    - Validates allowed characters
    """

    if not raw:
        raise ValueError("No sequence provided.")

    # Remove whitespace/newlines/tabs/spaces
    cleaned = re.sub(r"\s+", "", raw).upper()

    # Validate: A, C, G, T, (optionally N)
    if not re.fullmatch(r"[ACGTN]+", cleaned):
        invalid_chars = sorted(set(re.sub(r"[ACGTN]", "", cleaned)))
        raise ValueError(
            f"Sequence contains invalid characters: {', '.join(invalid_chars)}"
        )

    return cleaned

def find_binding_site(seq:str, primer:str) -> int:
    """Return start index where primer binds. Return None if not found."""
    idx = seq.find(primer)
    return idx if idx != -1 else None

def render_binding_line(seq:str, primer:str, pos:int)->str:
    if pos is None:
        return "No binding site found."
    left = "." * pos
    right = "." * (len(seq) - pos - len(primer))
    return left + primer + right

def window_sequence(seq, pos, primer_len, flank=25):
    """Return shortened window with ellipsis and primer start."""
    if pos is None:
        return None, None, None

    start = max(0, pos - flank)
    end = min(len(seq), pos + primer_len + flank)

    window = seq[start:end]

    left_offset = 0
    if start > 0:
        window = "..." + window
        left_offset = 3

    if end < len(seq):
        window = window + "..."

    primer_start_in_window = (pos - start) + left_offset
    return window, primer_start_in_window, primer_len


def render_windowed_line(seq_window, primer_start, primer_len):
    if seq_window is None:
        return "No binding site found."
    left = "." * primer_start
    right = "." * (len(seq_window) - primer_start - primer_len)
    primer_chars = seq_window[primer_start: primer_start + primer_len]
    return left + primer_chars + right

def highlight_binding(window: str, start: int, length: int) -> str:
    before = window[:start]
    bound = window[start:start + length]
    after = window[start + length:]

    return (
        f"{before}"
        f"<span class='bg-primary text-primary-content font-bold px-1 rounded'>"
        f"{bound}"
        f"</span>"
        f"{after}"
    )
