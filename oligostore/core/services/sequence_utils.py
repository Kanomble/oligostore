def reverse_complement(seq: str) -> str:
    complement = str.maketrans("ACGT", "TGCA")
    return seq.translate(complement)[::-1]