from Bio.Restriction import CommOnly
from django.db import transaction

from ..access import accessible_pcr_products, accessible_sequence_files
from ..models import CloningConstruct
from .ownership import assign_creator, grant_user_access
from .sequence_loader import load_sequences


def get_common_enzyme_choices():
    return sorted((str(enzyme), str(enzyme)) for enzyme in CommOnly)


def _get_enzyme_by_name(name: str):
    normalized = (name or "").strip()
    for enzyme in CommOnly:
        if str(enzyme) == normalized:
            return enzyme
    return None


def _find_site_positions(sequence: str, site: str):
    positions = []
    start = 0
    while True:
        index = sequence.find(site, start)
        if index == -1:
            break
        positions.append(index)
        start = index + 1
    return positions


def _resolve_sequence_file_asset(sequence_file):
    records = list(load_sequences(sequence_file.file.path, sequence_file.file_type))
    if not records:
        raise ValueError(f"Sequence file '{sequence_file.name}' does not contain any records.")
    record = records[0]
    if len(records) > 1:
        message = f"Sequence file '{sequence_file.name}' contains multiple records; using the first record '{record.id}'."
    else:
        message = None
    return {
        "source_type": CloningConstruct.SOURCE_SEQUENCE_FILE,
        "sequence_file": sequence_file,
        "pcr_product": None,
        "record_id": str(record.id),
        "name": sequence_file.name,
        "sequence": str(record.seq).upper(),
        "message": message,
    }


def _resolve_pcr_product_asset(pcr_product):
    return {
        "source_type": CloningConstruct.SOURCE_PCR_PRODUCT,
        "sequence_file": None,
        "pcr_product": pcr_product,
        "record_id": pcr_product.record_id,
        "name": pcr_product.name,
        "sequence": pcr_product.sequence,
        "message": None,
    }


def parse_asset_choice(choice: str):
    parts = str(choice or "").split(":", 1)
    if len(parts) != 2:
        raise ValueError("Invalid cloning asset selection.")
    source_type, object_id = parts
    try:
        return source_type, int(object_id)
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid cloning asset id.") from exc


def resolve_asset_choice(*, user, choice: str):
    source_type, object_id = parse_asset_choice(choice)
    if source_type == CloningConstruct.SOURCE_SEQUENCE_FILE:
        sequence_file = accessible_sequence_files(user).filter(id=object_id).first()
        if sequence_file is None:
            raise ValueError("Selected sequence file is not available.")
        return _resolve_sequence_file_asset(sequence_file)
    if source_type == CloningConstruct.SOURCE_PCR_PRODUCT:
        pcr_product = accessible_pcr_products(user).filter(id=object_id).first()
        if pcr_product is None:
            raise ValueError("Selected PCR product is not available.")
        return _resolve_pcr_product_asset(pcr_product)
    raise ValueError("Unknown cloning asset type.")


def _validate_restriction_ligation(vector_sequence: str, insert_sequence: str, left_enzyme_name: str, right_enzyme_name: str):
    messages = []
    left_enzyme = _get_enzyme_by_name(left_enzyme_name)
    right_enzyme = _get_enzyme_by_name(right_enzyme_name)

    if left_enzyme is None:
        messages.append(f"Unknown left enzyme: {left_enzyme_name}.")
    if right_enzyme is None:
        messages.append(f"Unknown right enzyme: {right_enzyme_name}.")
    if messages:
        return False, messages, ""

    if left_enzyme_name == right_enzyme_name:
        messages.append("Choose two different enzymes for this workflow.")
        return False, messages, ""

    left_site = str(getattr(left_enzyme, "site", "") or "")
    right_site = str(getattr(right_enzyme, "site", "") or "")
    left_hits = _find_site_positions(vector_sequence, left_site)
    right_hits = _find_site_positions(vector_sequence, right_site)

    if len(left_hits) != 1:
        messages.append(f"Vector must contain exactly one {left_enzyme_name} site; found {len(left_hits)}.")
    if len(right_hits) != 1:
        messages.append(f"Vector must contain exactly one {right_enzyme_name} site; found {len(right_hits)}.")
    if _find_site_positions(insert_sequence, left_site):
        messages.append(f"Insert contains an internal {left_enzyme_name} site.")
    if _find_site_positions(insert_sequence, right_site):
        messages.append(f"Insert contains an internal {right_enzyme_name} site.")
    if messages:
        return False, messages, ""

    left_start = left_hits[0]
    right_start = right_hits[0]
    left_end = left_start + len(left_site)

    if right_start <= left_end:
        messages.append("Restriction sites are overlapping or out of order in the vector.")
        return False, messages, ""

    assembled_sequence = (
        vector_sequence[:left_start]
        + insert_sequence
        + vector_sequence[right_start + len(right_site):]
    )
    return True, ["Construct assembled successfully."], assembled_sequence


@transaction.atomic
def create_cloning_construct(
    *,
    name: str,
    description: str,
    vector_asset,
    insert_asset,
    assembly_strategy: str,
    left_enzyme: str,
    right_enzyme: str,
    user,
):
    validation_messages = []
    if vector_asset.get("message"):
        validation_messages.append(vector_asset["message"])
    if insert_asset.get("message"):
        validation_messages.append(insert_asset["message"])

    is_valid, strategy_messages, assembled_sequence = _validate_restriction_ligation(
        vector_asset["sequence"],
        insert_asset["sequence"],
        left_enzyme,
        right_enzyme,
    )
    validation_messages.extend(strategy_messages)

    construct = CloningConstruct(
        name=name,
        description=description,
        vector_source_type=vector_asset["source_type"],
        vector_sequence_file=vector_asset["sequence_file"],
        vector_pcr_product=vector_asset["pcr_product"],
        vector_record_id=vector_asset["record_id"],
        insert_source_type=insert_asset["source_type"],
        insert_sequence_file=insert_asset["sequence_file"],
        insert_pcr_product=insert_asset["pcr_product"],
        insert_record_id=insert_asset["record_id"],
        assembly_strategy=assembly_strategy,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        assembled_sequence=assembled_sequence,
        is_valid=is_valid,
        validation_messages=validation_messages,
    )
    assign_creator(construct, user)
    construct.save()
    grant_user_access(construct, user)
    return construct
