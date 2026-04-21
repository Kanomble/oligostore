from typing import Iterable

from Bio import SeqIO


FORMAT_BY_FILE_TYPE = {
    "fasta": "fasta",
    "genbank": "genbank",
    "snapgene": "snapgene",
}


def resolve_sequence_format(file_type: str) -> str:
    normalized_file_type = str(file_type or "").strip().lower()
    sequence_format = FORMAT_BY_FILE_TYPE.get(normalized_file_type)
    if not sequence_format:
        raise ValueError(f"Unsupported file type: {file_type}")
    return sequence_format


def load_sequences(file_path: str, file_type: str) -> Iterable:
    """
    Returns an iterator of SeqRecord objects.
    """
    return SeqIO.parse(file_path, resolve_sequence_format(file_type))


def load_sequences_from_sequence_file(sequence_file) -> Iterable:
    return load_sequences(sequence_file.file.path, sequence_file.file_type)
