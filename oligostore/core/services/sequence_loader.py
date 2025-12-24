from Bio import SeqIO
from typing import Iterable


def load_sequences(file_path: str, file_type: str) -> Iterable:
    """
    Returns an iterator of SeqRecord objects.
    """
    if file_type == "fasta":
        return SeqIO.parse(file_path, "fasta")

    if file_type == "genbank":
        return SeqIO.parse(file_path, "genbank")

    raise ValueError(f"Unsupported file type: {file_type}")
