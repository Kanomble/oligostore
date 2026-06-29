from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio.Restriction import CommOnly
from Bio.Seq import Seq
from django.conf import settings
from django.db import transaction

from ..access import accessible_pcr_products, accessible_sequence_files
from ..models import CloningConstruct
from .ownership import assign_creator, grant_user_access
from .sequence_records import get_sequence_records
from .sequence_loader import load_sequences


@dataclass(frozen=True)
class CloningAssetChoice:
    source_type: str
    object_id: int | str
    record_id: Optional[str] = None

    @property
    def encoded_value(self) -> str:
        if self.record_id is None:
            return f"{self.source_type}:{self.object_id}"
        return f"{self.source_type}:{self.object_id}:{self.record_id}"


@dataclass(frozen=True)
class ResolvedCloningAsset:
    source_type: str
    sequence_file: Optional[object]
    pcr_product: Optional[object]
    template_name: Optional[str]
    record_id: str
    name: str
    sequence: str
    message: Optional[str] = None


@dataclass(frozen=True)
class ResolvedCloningAssets:
    vector_asset: ResolvedCloningAsset
    insert_asset: ResolvedCloningAsset


@dataclass(frozen=True)
class CloningValidationResult:
    is_valid: bool
    validation_messages: tuple[str, ...] = field(default_factory=tuple)
    assembled_sequence: str = ""
    insertion_start: Optional[int] = None
    inserted_length: int = 0
    vector_fragment_start: Optional[int] = None
    vector_fragment_end: Optional[int] = None
    insert_fragment_start: Optional[int] = None
    insert_fragment_end: Optional[int] = None


@dataclass(frozen=True)
class CloningCutSitePreview:
    enzyme_name: str
    site_sequence: str
    vector_cut_positions: tuple[int, ...] = field(default_factory=tuple)
    insert_recognition_site_positions: tuple[int, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class CloningJunctionContext:
    label: str
    left_context: str
    right_context: str

    @property
    def display(self) -> str:
        return f"{self.left_context}|{self.right_context}"


@dataclass(frozen=True)
class CloningConstructPreview:
    vector_asset: ResolvedCloningAsset
    insert_asset: ResolvedCloningAsset
    assembly_strategy: str
    left_enzyme: str
    right_enzyme: str
    assembled_sequence: str
    assembled_length: int
    is_valid: bool
    validation_messages: tuple[str, ...] = field(default_factory=tuple)
    cut_site_previews: tuple[CloningCutSitePreview, ...] = field(default_factory=tuple)
    vector_fragment_index: Optional[int] = None
    insert_fragment_index: Optional[int] = None


@dataclass(frozen=True)
class CloningConstructDetailDisplay:
    junction_context_window: int
    cut_site_previews: tuple[CloningCutSitePreview, ...] = field(default_factory=tuple)
    junction_contexts: tuple[CloningJunctionContext, ...] = field(default_factory=tuple)
    source_errors: tuple[str, ...] = field(default_factory=tuple)


TEMPLATE_SEQUENCE_FILENAMES = (
    "LvL25_without_J23100.gb",
    "tProm_leader_sequence_CRISPR_1.gb",
)


def build_sequence_file_asset_choice(*, sequence_file_id: int, record_id: str) -> CloningAssetChoice:
    return CloningAssetChoice(
        source_type=CloningConstruct.SOURCE_SEQUENCE_FILE,
        object_id=sequence_file_id,
        record_id=str(record_id),
    )


def build_pcr_product_asset_choice(*, pcr_product_id: int) -> CloningAssetChoice:
    return CloningAssetChoice(
        source_type=CloningConstruct.SOURCE_PCR_PRODUCT,
        object_id=pcr_product_id,
    )


def build_template_asset_choice(*, template_name: str, record_id: str) -> CloningAssetChoice:
    return CloningAssetChoice(
        source_type=CloningConstruct.SOURCE_TEMPLATE,
        object_id=str(template_name),
        record_id=str(record_id),
    )


def _template_directory() -> Path:
    return Path(settings.MEDIA_ROOT) / "sequence_files"


def _resolve_template_path(template_name: str) -> Path:
    normalized_template_name = str(template_name or "").strip()
    if not normalized_template_name:
        raise ValueError("A template file name is required.")

    template_directory = _template_directory()
    if not template_directory.exists():
        raise ValueError("The template media directory is not available.")

    for candidate in template_directory.iterdir():
        if not candidate.is_file():
            continue
        if candidate.name.lower() == normalized_template_name.lower() or candidate.stem.lower() == normalized_template_name.lower():
            return candidate

    raise ValueError(f"Template '{normalized_template_name}' is not available.")


def _load_template_records(template_name: str):
    template_path = _resolve_template_path(template_name)
    try:
        template_format = {
            ".gb": "genbank",
            ".gbk": "genbank",
            ".gbff": "genbank",
            ".fasta": "fasta",
            ".fa": "fasta",
            ".dna": "snapgene",
        }.get(template_path.suffix.lower())
        if template_format is None:
            raise ValueError(f"Template '{template_path.name}' uses an unsupported file type.")
        records = list(load_sequences(str(template_path), template_format))
    except Exception as exc:
        raise ValueError(f"Template '{template_name}' could not be parsed for cloning.") from exc
    if not records:
        raise ValueError(f"Template '{template_name}' does not contain any records.")
    return template_path, records


def get_template_asset_choices():
    choices = []
    for template_name in TEMPLATE_SEQUENCE_FILENAMES:
        try:
            template_path, records = _load_template_records(template_name)
        except ValueError:
            continue
        for record in records:
            asset_choice = build_template_asset_choice(
                template_name=template_path.name,
                record_id=str(record.id),
            )
            choices.append(
                (
                    asset_choice.encoded_value,
                    (
                        f"Template | {template_path.stem} | "
                        f"record {record.id} | {len(record.seq)} bp | template file"
                    ),
                )
            )
    return choices


def get_detected_enzyme_choices(*, user, selected_asset_choice=None):
    if not selected_asset_choice:
        return []

    try:
        asset = resolve_asset_choice(user=user, choice=selected_asset_choice)
    except ValueError:
        return []

    detected = set()
    for sequence in [asset.sequence]:
        try:
            results = CommOnly.search(Seq(sequence), linear=True)
        except Exception:
            continue
        for enzyme, cut_positions in results.items():
            if cut_positions:
                detected.add(str(enzyme))

    return [(name, name) for name in sorted(detected)]


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


def _build_cut_site_preview(*, vector_sequence: str, insert_sequence: str, enzyme_name: str) -> CloningCutSitePreview:
    enzyme = _get_enzyme_by_name(enzyme_name)
    if enzyme is None:
        return CloningCutSitePreview(
            enzyme_name=enzyme_name,
            site_sequence="",
        )
    site_sequence = str(getattr(enzyme, "site", "") or "")
    return CloningCutSitePreview(
        enzyme_name=enzyme_name,
        site_sequence=site_sequence,
        vector_cut_positions=tuple(position + 1 for position in _find_cut_positions(vector_sequence, enzyme)),
        insert_recognition_site_positions=tuple(
            position + 1 for position in _find_site_positions(insert_sequence, site_sequence)
        ),
    )


def _find_cut_positions(sequence: str, enzyme):
    try:
        results = CommOnly.search(Seq(sequence), linear=True)
    except Exception:
        return []
    return [int(position) - 1 for position in results.get(enzyme, [])]


@dataclass(frozen=True)
class DigestedFragment:
    index: int
    start: int
    end: int
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)


def _digest_sequence_fragments(sequence: str, enzyme) -> tuple[DigestedFragment, ...]:
    cut_positions = sorted(
        set(position for position in _find_cut_positions(sequence, enzyme) if 0 < position < len(sequence))
    )
    boundaries = [0, *cut_positions, len(sequence)]
    fragments = []
    for index, (start, end) in enumerate(zip(boundaries, boundaries[1:]), start=1):
        if end <= start:
            continue
        fragments.append(
            DigestedFragment(
                index=index,
                start=start,
                end=end,
                sequence=sequence[start:end],
            )
        )
    return tuple(fragments)


def _select_digest_fragment(
    *,
    sequence: str,
    enzyme,
    fragment_index: Optional[int],
    minimum_length: int = 1,
    label: str,
) -> DigestedFragment:
    fragments = _digest_sequence_fragments(sequence, enzyme)
    if not fragments:
        raise ValueError(f"{label} does not yield any fragments after digestion.")

    eligible_fragments = [fragment for fragment in fragments if fragment.length >= minimum_length]
    if fragment_index is not None:
        selected_fragment = next((fragment for fragment in eligible_fragments if fragment.index == fragment_index), None)
        if selected_fragment is None:
            raise ValueError(f"Selected {label.lower()} fragment is not available after digestion.")
        return selected_fragment

    if len(eligible_fragments) == 1:
        return eligible_fragments[0]
    if not eligible_fragments:
        raise ValueError(f"{label} does not produce any fragments with the selected enzyme.")
    raise ValueError(
        f"Select which {label.lower()} fragment to use after digestion; {len(eligible_fragments)} fragments are available."
    )


def build_digest_fragment_choices(*, sequence: str, enzyme_name: str, minimum_length: int = 1):
    enzyme = _get_enzyme_by_name(enzyme_name)
    if enzyme is None:
        return []
    fragments = _digest_sequence_fragments(sequence, enzyme)
    eligible_fragments = sorted(
        (fragment for fragment in fragments if fragment.length >= minimum_length),
        key=lambda fragment: (-fragment.length, fragment.start, fragment.index),
    )
    return [
        (
            str(fragment.index),
            f"Fragment {fragment.index} | {fragment.length} bp | bases {fragment.start + 1}-{fragment.end}",
        )
        for fragment in eligible_fragments
    ]


def _normalize_fragment_index(value) -> Optional[int]:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid fragment selection.") from exc


def _load_sequence_file_records(sequence_file):
    try:
        records = list(get_sequence_records(sequence_file, load_sequences))
    except Exception as exc:
        raise ValueError(
            f"Sequence file '{sequence_file.name}' could not be parsed for cloning."
        ) from exc
    if not records:
        raise ValueError(f"Sequence file '{sequence_file.name}' does not contain any records.")
    return records


def _resolve_sequence_file_asset(sequence_file):
    records = _load_sequence_file_records(sequence_file)
    if len(records) > 1:
        raise ValueError(
            f"Sequence file '{sequence_file.name}' contains multiple records. Select a specific record for cloning."
        )
    record = records[0]
    return ResolvedCloningAsset(
        source_type=CloningConstruct.SOURCE_SEQUENCE_FILE,
        sequence_file=sequence_file,
        pcr_product=None,
        template_name=None,
        record_id=str(record.id),
        name=sequence_file.name,
        sequence=str(record.seq).upper(),
    )


def _resolve_sequence_file_record_asset(sequence_file, record_id: str):
    records = _load_sequence_file_records(sequence_file)
    normalized_record_id = str(record_id or "").strip()
    if not normalized_record_id:
        raise ValueError("A sequence record must be selected for cloning.")
    record = next((record for record in records if str(record.id) == normalized_record_id), None)
    if record is None:
        raise ValueError(
            f"Record '{normalized_record_id}' is not available in sequence file '{sequence_file.name}'."
        )
    return ResolvedCloningAsset(
        source_type=CloningConstruct.SOURCE_SEQUENCE_FILE,
        sequence_file=sequence_file,
        pcr_product=None,
        template_name=None,
        record_id=str(record.id),
        name=sequence_file.name,
        sequence=str(record.seq).upper(),
    )


def _resolve_pcr_product_asset(pcr_product):
    return ResolvedCloningAsset(
        source_type=CloningConstruct.SOURCE_PCR_PRODUCT,
        sequence_file=None,
        pcr_product=pcr_product,
        template_name=None,
        record_id=pcr_product.record_id,
        name=pcr_product.name,
        sequence=pcr_product.sequence,
    )


def _resolve_template_asset(template_name, record_id: Optional[str] = None):
    template_path, records = _load_template_records(template_name)
    normalized_record_id = str(record_id or "").strip()
    if normalized_record_id:
        record = next((record for record in records if str(record.id) == normalized_record_id), None)
        if record is None:
            raise ValueError(
                f"Record '{normalized_record_id}' is not available in template '{template_path.name}'."
            )
    elif len(records) == 1:
        record = records[0]
    else:
        raise ValueError(
            f"Template '{template_path.name}' contains multiple records. Select a specific record for cloning."
        )

    return ResolvedCloningAsset(
        source_type=CloningConstruct.SOURCE_TEMPLATE,
        sequence_file=None,
        pcr_product=None,
        template_name=template_path.name,
        record_id=str(record.id),
        name=template_path.name,
        sequence=str(record.seq).upper(),
    )


def _resolve_construct_asset(*, source_type, sequence_file, pcr_product, template_name, record_id, label) -> ResolvedCloningAsset:
    if source_type == CloningConstruct.SOURCE_PCR_PRODUCT:
        if pcr_product is None:
            raise ValueError(f"{label} PCR product is no longer available.")
        return _resolve_pcr_product_asset(pcr_product)
    if source_type == CloningConstruct.SOURCE_SEQUENCE_FILE:
        if sequence_file is None:
            raise ValueError(f"{label} sequence file is no longer available.")
        if record_id:
            return _resolve_sequence_file_record_asset(sequence_file, record_id)
        return _resolve_sequence_file_asset(sequence_file)
    if source_type == CloningConstruct.SOURCE_TEMPLATE:
        if not template_name:
            raise ValueError(f"{label} template is no longer available.")
        return _resolve_template_asset(template_name, record_id)
    raise ValueError(f"Unknown {label.lower()} source type.")


def parse_asset_choice(choice: str):
    parts = str(choice or "").split(":", 2)
    if len(parts) not in {2, 3}:
        raise ValueError("Invalid cloning asset selection.")
    source_type, object_id = parts[:2]
    record_id = parts[2] if len(parts) == 3 else None
    try:
        if source_type == CloningConstruct.SOURCE_TEMPLATE:
            return CloningAssetChoice(
                source_type=source_type,
                object_id=object_id,
                record_id=record_id,
            )
        return CloningAssetChoice(
            source_type=source_type,
            object_id=int(object_id),
            record_id=record_id,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("Invalid cloning asset id.") from exc


def resolve_asset_choice(*, user, choice: str) -> ResolvedCloningAsset:
    asset_choice = parse_asset_choice(choice)
    if asset_choice.source_type == CloningConstruct.SOURCE_SEQUENCE_FILE:
        sequence_file = accessible_sequence_files(user).filter(id=asset_choice.object_id).first()
        if sequence_file is None:
            raise ValueError("Selected sequence file is not available.")
        if asset_choice.record_id is not None:
            return _resolve_sequence_file_record_asset(sequence_file, asset_choice.record_id)
        return _resolve_sequence_file_asset(sequence_file)
    if asset_choice.source_type == CloningConstruct.SOURCE_PCR_PRODUCT:
        pcr_product = accessible_pcr_products(user).filter(id=asset_choice.object_id).first()
        if pcr_product is None:
            raise ValueError("Selected PCR product is not available.")
        return _resolve_pcr_product_asset(pcr_product)
    if asset_choice.source_type == CloningConstruct.SOURCE_TEMPLATE:
        return _resolve_template_asset(asset_choice.object_id, asset_choice.record_id)
    raise ValueError("Unknown cloning asset type.")


def resolve_cloning_assets(*, user, vector_asset_choice: str, insert_asset_choice: str) -> ResolvedCloningAssets:
    return ResolvedCloningAssets(
        vector_asset=resolve_asset_choice(user=user, choice=vector_asset_choice),
        insert_asset=resolve_asset_choice(user=user, choice=insert_asset_choice),
    )


def _invalid_cloning_result(*messages: str) -> CloningValidationResult:
    return CloningValidationResult(is_valid=False, validation_messages=tuple(messages))


def _validate_same_enzyme_fragment_ligation(
    *,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_name: str,
    enzyme,
    vector_fragment_index: Optional[int],
    insert_fragment_index: Optional[int],
) -> CloningValidationResult:
    try:
        vector_fragment = _select_digest_fragment(
            sequence=vector_sequence,
            enzyme=enzyme,
            fragment_index=vector_fragment_index,
            label="Vector",
        )
        insert_fragment = _select_digest_fragment(
            sequence=insert_sequence,
            enzyme=enzyme,
            fragment_index=insert_fragment_index,
            label="Insert",
        )
    except ValueError as exc:
        return _invalid_cloning_result(str(exc))

    assembled_sequence = vector_fragment.sequence + insert_fragment.sequence
    return CloningValidationResult(
        is_valid=True,
        validation_messages=(
            f"Construct assembled successfully using {enzyme_name} fragment selection (vector fragment {vector_fragment.index} and insert fragment {insert_fragment.index}).",
        ),
        assembled_sequence=assembled_sequence,
        insertion_start=len(vector_fragment.sequence),
        inserted_length=len(insert_fragment.sequence),
        vector_fragment_start=vector_fragment.start,
        vector_fragment_end=vector_fragment.end,
        insert_fragment_start=insert_fragment.start,
        insert_fragment_end=insert_fragment.end,
    )


def _validate_same_enzyme_ligation(
    *,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_name: str,
    enzyme,
) -> CloningValidationResult:
    vector_cut_hits = _find_cut_positions(vector_sequence, enzyme)
    insert_cut_hits = _find_cut_positions(insert_sequence, enzyme)

    if len(vector_cut_hits) == 1 and len(insert_cut_hits) == 0:
        insertion_point = vector_cut_hits[0]
        assembled_sequence = (
            vector_sequence[:insertion_point]
            + insert_sequence
            + vector_sequence[insertion_point:]
        )
        return CloningValidationResult(
            is_valid=True,
            validation_messages=(
                f"Construct assembled successfully using {enzyme_name} single-cut insertion mode.",
            ),
            assembled_sequence=assembled_sequence,
            insertion_start=insertion_point,
            inserted_length=len(insert_sequence),
        )

    if len(vector_cut_hits) == 2 and len(insert_cut_hits) == 2:
        vector_left_cut, vector_right_cut = sorted(vector_cut_hits)
        insert_left_cut, insert_right_cut = sorted(insert_cut_hits)
        messages = []
        if vector_right_cut <= vector_left_cut:
            messages.append(f"{enzyme_name} cut positions in the vector are overlapping or out of order.")
        if insert_right_cut <= insert_left_cut:
            messages.append(f"{enzyme_name} cut positions in the insert are overlapping or out of order.")
        if messages:
            return _invalid_cloning_result(*messages)

        insert_fragment = insert_sequence[insert_left_cut:insert_right_cut]
        assembled_sequence = (
            vector_sequence[:vector_left_cut]
            + insert_fragment
            + vector_sequence[vector_right_cut:]
        )
        return CloningValidationResult(
            is_valid=True,
            validation_messages=(
                f"Construct assembled successfully using {enzyme_name} double-cut fragment replacement mode.",
            ),
            assembled_sequence=assembled_sequence,
            insertion_start=vector_left_cut,
            inserted_length=len(insert_fragment),
        )

    messages = []
    if len(vector_cut_hits) not in {1, 2}:
        messages.append(
            f"Vector must contain either 1 cut site (direct insertion) or 2 cut sites (fragment replacement) for same-enzyme cloning with {enzyme_name}; found {len(vector_cut_hits)}."
        )
    if len(insert_cut_hits) not in {0, 2}:
        messages.append(
            f"Insert must contain either 0 cut sites (insert the full asset) or 2 cut sites (excise the fragment between them) for same-enzyme cloning with {enzyme_name}; found {len(insert_cut_hits)}."
        )
    if len(vector_cut_hits) in {1, 2} and len(insert_cut_hits) in {0, 2}:
        messages.append(
            f"The selected same-enzyme topology is not supported. Use vector 1 cut / insert 0 cuts for direct insertion, or vector 2 cuts / insert 2 cuts for fragment replacement. Found vector {len(vector_cut_hits)} / insert {len(insert_cut_hits)}."
        )
    return _invalid_cloning_result(*messages)


def _validate_two_enzyme_ligation(
    *,
    vector_sequence: str,
    insert_sequence: str,
    left_enzyme_name: str,
    right_enzyme_name: str,
    left_enzyme,
    right_enzyme,
) -> CloningValidationResult:
    messages = []
    left_site = str(getattr(left_enzyme, "site", "") or "")
    right_site = str(getattr(right_enzyme, "site", "") or "")
    left_hits = _find_cut_positions(vector_sequence, left_enzyme)
    right_hits = _find_cut_positions(vector_sequence, right_enzyme)

    if len(left_hits) != 1:
        messages.append(f"Vector must contain exactly one {left_enzyme_name} cut site; found {len(left_hits)}.")
    if len(right_hits) != 1:
        messages.append(f"Vector must contain exactly one {right_enzyme_name} cut site; found {len(right_hits)}.")
    if _find_site_positions(insert_sequence, left_site):
        messages.append(f"Insert contains an internal {left_enzyme_name} site.")
    if _find_site_positions(insert_sequence, right_site):
        messages.append(f"Insert contains an internal {right_enzyme_name} site.")
    if messages:
        return _invalid_cloning_result(*messages)

    left_cut = left_hits[0]
    right_cut = right_hits[0]
    if right_cut <= left_cut:
        return _invalid_cloning_result("Restriction sites are overlapping or out of order in the vector.")

    assembled_sequence = (
        vector_sequence[:left_cut]
        + insert_sequence
        + vector_sequence[right_cut:]
    )
    return CloningValidationResult(
        is_valid=True,
        validation_messages=("Construct assembled successfully using two-enzyme cloning.",),
        assembled_sequence=assembled_sequence,
        insertion_start=left_cut,
        inserted_length=len(insert_sequence),
    )


def _validate_restriction_ligation(
    vector_sequence: str,
    insert_sequence: str,
    left_enzyme_name: str,
    right_enzyme_name: str,
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
) -> CloningValidationResult:
    left_enzyme = _get_enzyme_by_name(left_enzyme_name)
    right_enzyme = _get_enzyme_by_name(right_enzyme_name)

    messages = []
    if left_enzyme is None:
        messages.append(f"Unknown left enzyme: {left_enzyme_name}.")
    if right_enzyme is None:
        messages.append(f"Unknown right enzyme: {right_enzyme_name}.")
    if messages:
        return _invalid_cloning_result(*messages)

    if left_enzyme_name == right_enzyme_name:
        if vector_fragment_index is not None or insert_fragment_index is not None:
            return _validate_same_enzyme_fragment_ligation(
                vector_sequence=vector_sequence,
                insert_sequence=insert_sequence,
                enzyme_name=left_enzyme_name,
                enzyme=left_enzyme,
                vector_fragment_index=vector_fragment_index,
                insert_fragment_index=insert_fragment_index,
            )
        # Preserve the legacy same-enzyme topology for already-saved constructs that do not store fragment indexes.
        return _validate_same_enzyme_ligation(
            vector_sequence=vector_sequence,
            insert_sequence=insert_sequence,
            enzyme_name=left_enzyme_name,
            enzyme=left_enzyme,
        )

    return _validate_two_enzyme_ligation(
        vector_sequence=vector_sequence,
        insert_sequence=insert_sequence,
        left_enzyme_name=left_enzyme_name,
        right_enzyme_name=right_enzyme_name,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
    )


def preview_cloning_construct(
    *,
    vector_asset: ResolvedCloningAsset,
    insert_asset: ResolvedCloningAsset,
    assembly_strategy: str,
    left_enzyme: str,
    right_enzyme: str,
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
) -> CloningConstructPreview:
    if assembly_strategy != CloningConstruct.STRATEGY_RESTRICTION_LIGATION:
        raise ValueError("Only restriction ligation is currently supported.")

    vector_fragment_index = _normalize_fragment_index(vector_fragment_index)
    insert_fragment_index = _normalize_fragment_index(insert_fragment_index)

    validation_messages = []
    if vector_asset.message:
        validation_messages.append(vector_asset.message)
    if insert_asset.message:
        validation_messages.append(insert_asset.message)

    validation_result = _validate_restriction_ligation(
        vector_asset.sequence,
        insert_asset.sequence,
        left_enzyme,
        right_enzyme,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
    )
    validation_messages.extend(validation_result.validation_messages)
    unique_enzyme_names = []
    for enzyme_name in (left_enzyme, right_enzyme):
        if enzyme_name and enzyme_name not in unique_enzyme_names:
            unique_enzyme_names.append(enzyme_name)

    return CloningConstructPreview(
        vector_asset=vector_asset,
        insert_asset=insert_asset,
        assembly_strategy=assembly_strategy,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
        assembled_sequence=validation_result.assembled_sequence,
        assembled_length=len(validation_result.assembled_sequence or ""),
        is_valid=validation_result.is_valid,
        validation_messages=tuple(validation_messages),
        cut_site_previews=tuple(
            _build_cut_site_preview(
                vector_sequence=vector_asset.sequence,
                insert_sequence=insert_asset.sequence,
                enzyme_name=enzyme_name,
            )
            for enzyme_name in unique_enzyme_names
        ),
    )


def _build_junction_context(*, sequence: str, boundary_index: int, window: int, label: str) -> CloningJunctionContext:
    left_context = sequence[max(0, boundary_index - window):boundary_index]
    right_context = sequence[boundary_index:boundary_index + window]
    return CloningJunctionContext(
        label=label,
        left_context=left_context,
        right_context=right_context,
    )


def _serialize_detail_display(
    detail_display: CloningConstructDetailDisplay,
) -> dict:
    return {
        "junction_context_window": detail_display.junction_context_window,
        "cut_site_previews": [
            {
                "enzyme_name": preview.enzyme_name,
                "site_sequence": preview.site_sequence,
                "vector_cut_positions": list(preview.vector_cut_positions),
                "insert_recognition_site_positions": list(preview.insert_recognition_site_positions),
            }
            for preview in detail_display.cut_site_previews
        ],
        "junction_contexts": [
            {
                "label": junction.label,
                "left_context": junction.left_context,
                "right_context": junction.right_context,
            }
            for junction in detail_display.junction_contexts
        ],
        "source_errors": list(detail_display.source_errors),
    }


def _deserialize_detail_display(snapshot: dict) -> CloningConstructDetailDisplay:
    return CloningConstructDetailDisplay(
        junction_context_window=int(
            snapshot.get("junction_context_window", CloningConstruct.JUNCTION_CONTEXT_WINDOW)
        ),
        cut_site_previews=tuple(
            CloningCutSitePreview(
                enzyme_name=str(preview.get("enzyme_name", "")),
                site_sequence=str(preview.get("site_sequence", "")),
                vector_cut_positions=tuple(
                    int(position) for position in preview.get("vector_cut_positions", [])
                ),
                insert_recognition_site_positions=tuple(
                    int(position)
                    for position in preview.get(
                        "insert_recognition_site_positions",
                        preview.get("insert_site_positions", []),
                    )
                ),
            )
            for preview in snapshot.get("cut_site_previews", [])
        ),
        junction_contexts=tuple(
            CloningJunctionContext(
                label=str(junction.get("label", "")),
                left_context=str(junction.get("left_context", "")),
                right_context=str(junction.get("right_context", "")),
            )
            for junction in snapshot.get("junction_contexts", [])
        ),
        source_errors=tuple(str(message) for message in snapshot.get("source_errors", [])),
    )


def _build_detail_display_from_preview(
    preview_data: CloningConstructPreview,
) -> CloningConstructDetailDisplay:
    junction_contexts = tuple()
    if preview_data.assembled_sequence:
        validation_result = _validate_restriction_ligation(
            preview_data.vector_asset.sequence,
            preview_data.insert_asset.sequence,
            preview_data.left_enzyme,
            preview_data.right_enzyme,
            vector_fragment_index=preview_data.vector_fragment_index,
            insert_fragment_index=preview_data.insert_fragment_index,
        )
        if validation_result.is_valid and validation_result.insertion_start is not None:
            right_boundary = (
                validation_result.insertion_start + validation_result.inserted_length
            )
            if right_boundary <= len(preview_data.assembled_sequence):
                junction_contexts = (
                    _build_junction_context(
                        sequence=preview_data.assembled_sequence,
                        boundary_index=validation_result.insertion_start,
                        window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
                        label="Vector -> insert junction",
                    ),
                    _build_junction_context(
                        sequence=preview_data.assembled_sequence,
                        boundary_index=right_boundary,
                        window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
                        label="Insert -> vector junction",
                    ),
                )

    return CloningConstructDetailDisplay(
        junction_context_window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
        cut_site_previews=preview_data.cut_site_previews,
        junction_contexts=junction_contexts,
        source_errors=tuple(),
    )


def build_cloning_construct_detail_display(
    construct: CloningConstruct,
) -> CloningConstructDetailDisplay:
    if construct.detail_display_snapshot:
        return _deserialize_detail_display(construct.detail_display_snapshot)

    source_errors = []
    try:
        vector_asset = _resolve_construct_asset(
            source_type=construct.vector_source_type,
            sequence_file=construct.vector_sequence_file,
            pcr_product=construct.vector_pcr_product,
            template_name=construct.vector_template_name,
            record_id=construct.vector_record_id,
            label="Vector",
        )
    except ValueError as exc:
        source_errors.append(str(exc))
        vector_asset = None

    try:
        insert_asset = _resolve_construct_asset(
            source_type=construct.insert_source_type,
            sequence_file=construct.insert_sequence_file,
            pcr_product=construct.insert_pcr_product,
            template_name=construct.insert_template_name,
            record_id=construct.insert_record_id,
            label="Insert",
        )
    except ValueError as exc:
        source_errors.append(str(exc))
        insert_asset = None

    if vector_asset is None or insert_asset is None:
        return CloningConstructDetailDisplay(
            junction_context_window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
            source_errors=tuple(source_errors),
        )

    enzyme_names = []
    for enzyme_name in (construct.left_enzyme, construct.right_enzyme):
        if enzyme_name and enzyme_name not in enzyme_names:
            enzyme_names.append(enzyme_name)

    cut_site_previews = tuple(
        _build_cut_site_preview(
            vector_sequence=vector_asset.sequence,
            insert_sequence=insert_asset.sequence,
            enzyme_name=enzyme_name,
        )
        for enzyme_name in enzyme_names
    )

    junction_contexts = tuple()
    if (
        construct.assembly_strategy == CloningConstruct.STRATEGY_RESTRICTION_LIGATION
        and construct.assembled_sequence
    ):
        validation_result = _validate_restriction_ligation(
            vector_asset.sequence,
            insert_asset.sequence,
            construct.left_enzyme,
            construct.right_enzyme,
            vector_fragment_index=construct.vector_fragment_index,
            insert_fragment_index=construct.insert_fragment_index,
        )
        if validation_result.is_valid and validation_result.insertion_start is not None:
            right_boundary = (
                validation_result.insertion_start + validation_result.inserted_length
            )
            if right_boundary <= len(construct.assembled_sequence):
                junction_contexts = (
                    _build_junction_context(
                        sequence=construct.assembled_sequence,
                        boundary_index=validation_result.insertion_start,
                        window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
                        label="Vector -> insert junction",
                    ),
                    _build_junction_context(
                        sequence=construct.assembled_sequence,
                        boundary_index=right_boundary,
                        window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
                        label="Insert -> vector junction",
                    ),
                )

    return CloningConstructDetailDisplay(
        junction_context_window=CloningConstruct.JUNCTION_CONTEXT_WINDOW,
        cut_site_previews=cut_site_previews,
        junction_contexts=junction_contexts,
        source_errors=tuple(source_errors),
    )


@transaction.atomic
def save_cloning_construct(*, name: str, description: str, preview_data: CloningConstructPreview, user):
    detail_display = _build_detail_display_from_preview(preview_data)
    construct = CloningConstruct(
        name=name,
        description=description,
        vector_source_type=preview_data.vector_asset.source_type,
        vector_sequence_file=preview_data.vector_asset.sequence_file,
        vector_template_name=preview_data.vector_asset.template_name or "",
        vector_pcr_product=preview_data.vector_asset.pcr_product,
        vector_record_id=preview_data.vector_asset.record_id,
        vector_fragment_index=preview_data.vector_fragment_index,
        insert_source_type=preview_data.insert_asset.source_type,
        insert_sequence_file=preview_data.insert_asset.sequence_file,
        insert_template_name=preview_data.insert_asset.template_name or "",
        insert_pcr_product=preview_data.insert_asset.pcr_product,
        insert_record_id=preview_data.insert_asset.record_id,
        insert_fragment_index=preview_data.insert_fragment_index,
        assembly_strategy=preview_data.assembly_strategy,
        left_enzyme=preview_data.left_enzyme,
        right_enzyme=preview_data.right_enzyme,
        assembled_sequence=preview_data.assembled_sequence,
        is_valid=preview_data.is_valid,
        validation_messages=list(preview_data.validation_messages),
        detail_display_snapshot=_serialize_detail_display(detail_display),
    )
    assign_creator(construct, user)
    construct.save()
    grant_user_access(construct, user)
    return construct


def create_cloning_construct(
    *,
    name: str,
    description: str,
    vector_asset: ResolvedCloningAsset,
    insert_asset: ResolvedCloningAsset,
    assembly_strategy: str,
    left_enzyme: str,
    right_enzyme: str,
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
    user,
) -> CloningConstruct:
    preview_data = preview_cloning_construct(
        vector_asset=vector_asset,
        insert_asset=insert_asset,
        assembly_strategy=assembly_strategy,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
    )
    return save_cloning_construct(
        name=name,
        description=description,
        preview_data=preview_data,
        user=user,
    )
