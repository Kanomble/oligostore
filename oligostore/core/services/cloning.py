from dataclasses import dataclass, field
import math
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
    is_circular: bool = False
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
    vector_fragment_segments: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    insert_fragment_segments: tuple[tuple[int, int], ...] = field(default_factory=tuple)
    vector_assembly_segments: tuple[tuple[int, int, int], ...] = field(default_factory=tuple)
    insert_assembly_segments: tuple[tuple[int, int, int], ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class CloningCutSitePreview:
    enzyme_name: str
    site_sequence: str
    vector_cut_positions: tuple[int, ...] = field(default_factory=tuple)
    insert_recognition_site_positions: tuple[int, ...] = field(default_factory=tuple)


@dataclass(frozen=True)
class CloningRestrictionCutMarker:
    enzyme_name: str
    position: int

    @property
    def display_position(self) -> int:
        return self.position + 1

    @property
    def label(self) -> str:
        return f"{self.enzyme_name} cut {self.display_position}"


@dataclass(frozen=True)
class CloningDoubleStrandCutView:
    enzyme_name: str
    top_cut_position: int
    bottom_cut_position: int
    top_line: str
    bottom_line: str
    overhang_sequence: str
    overhang_type: str

    @property
    def top_cut_display(self) -> int:
        return self.top_cut_position + 1

    @property
    def bottom_cut_display(self) -> int:
        return self.bottom_cut_position + 1


@dataclass(frozen=True)
class CloningDigestSegment:
    start: int
    end: int
    label: str

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

    @property
    def display_start(self) -> int:
        return self.start + 1

    @property
    def display_end(self) -> int:
        return self.end


@dataclass(frozen=True)
class CloningDigestFragmentOption:
    field_name: str
    fragment_index: int
    start: int
    end: int
    sequence_length: int
    selected: bool = False
    wraps_origin: bool = False

    @property
    def length(self) -> int:
        if self.wraps_origin:
            return max(0, self.sequence_length - self.start + self.end)
        return max(0, self.end - self.start)

    @property
    def display_start(self) -> int:
        return self.start + 1

    @property
    def display_end(self) -> int:
        return self.end

    @property
    def label(self) -> str:
        if self.wraps_origin:
            return (
                f"Fragment {self.fragment_index}: "
                f"{self.display_start}-{self.sequence_length} + 1-{self.display_end} ({self.length} bp)"
            )
        return f"Fragment {self.fragment_index}: {self.display_start}-{self.display_end} ({self.length} bp)"


@dataclass(frozen=True)
class CloningVisualRegion:
    label: str
    field_name: str
    fragment_index: Optional[int]
    start: int
    end: int
    sequence_length: int
    selected: bool = False
    selectable: bool = False
    wraps_origin: bool = False

    @property
    def length(self) -> int:
        if self.wraps_origin:
            return max(0, self.sequence_length - self.start + self.end)
        return max(0, self.end - self.start)

    @property
    def display_start(self) -> int:
        return self.start + 1

    @property
    def display_end(self) -> int:
        return self.end

    @property
    def x_percent(self) -> str:
        if self.sequence_length <= 0:
            return "0.000"
        return f"{(self.start / self.sequence_length) * 100:.3f}"

    @property
    def width_percent(self) -> str:
        if self.sequence_length <= 0:
            return "0.000"
        return f"{max(1.25, (self.length / self.sequence_length) * 100):.3f}"

    @property
    def dash_length(self) -> str:
        if self.sequence_length <= 0:
            return "0.000"
        return f"{max(1.0, (self.length / self.sequence_length) * 100):.3f}"

    @property
    def dash_offset(self) -> str:
        if self.sequence_length <= 0:
            return "0.000"
        return f"{-(self.start / self.sequence_length) * 100:.3f}"

    @property
    def title(self) -> str:
        if self.wraps_origin:
            coordinates = f"{self.display_start}-{self.sequence_length} + 1-{self.display_end}"
        else:
            coordinates = f"{self.display_start}-{self.display_end}"
        state = "selected" if self.selected else "available"
        return f"{self.label}: {coordinates}, {self.length} bp, {state}"


@dataclass(frozen=True)
class CloningVisualRestrictionSite:
    enzyme_name: str
    site_sequence: str
    role: str
    position: int
    sequence_length: int
    vector_cut_count: int
    insert_site_count: int
    selected_left: bool = False
    selected_right: bool = False
    compatible: bool = True

    @property
    def display_position(self) -> int:
        return self.position + 1

    @property
    def is_selected(self) -> bool:
        return self.selected_left or self.selected_right

    @property
    def selected_label(self) -> str:
        if self.selected_left and self.selected_right:
            return "left and right enzyme"
        if self.selected_left:
            return "left enzyme"
        if self.selected_right:
            return "right enzyme"
        return ""

    @property
    def is_unique(self) -> bool:
        return self.vector_cut_count == 1

    @property
    def x_percent(self) -> str:
        if self.sequence_length <= 0:
            return "0.000"
        return f"{(self.position / self.sequence_length) * 100:.3f}"

    @property
    def marker_x(self) -> str:
        angle = self._angle_radians
        return f"{50 + 38 * math.cos(angle):.3f}"

    @property
    def marker_y(self) -> str:
        angle = self._angle_radians
        return f"{50 + 38 * math.sin(angle):.3f}"

    @property
    def label_x(self) -> str:
        angle = self._angle_radians
        return f"{50 + 47 * math.cos(angle):.3f}"

    @property
    def label_y(self) -> str:
        angle = self._angle_radians
        return f"{50 + 47 * math.sin(angle):.3f}"

    @property
    def text_anchor(self) -> str:
        angle = self._angle_radians
        x_component = math.cos(angle)
        if x_component > 0.25:
            return "start"
        if x_component < -0.25:
            return "end"
        return "middle"

    @property
    def stable_id(self) -> str:
        enzyme_part = "".join(
            character if character.isalnum() else "-"
            for character in self.enzyme_name.lower()
        ).strip("-")
        return f"{self.role}-{enzyme_part}-{self.position}"

    @property
    def compatibility_label(self) -> str:
        if self.vector_cut_count == 1 and self.insert_site_count == 0:
            return "unique vector cut, no insert site"
        if self.vector_cut_count == 1:
            return f"unique vector cut, {self.insert_site_count} insert site(s)"
        return f"{self.vector_cut_count} vector cuts; fragment selection may be needed"

    @property
    def title(self) -> str:
        selection = f", selected as {self.selected_label}" if self.is_selected else ""
        return (
            f"{self.enzyme_name} at {self.role} base {self.display_position}; "
            f"site {self.site_sequence or 'unavailable'}; {self.compatibility_label}{selection}"
        )

    @property
    def _angle_radians(self) -> float:
        if self.sequence_length <= 0:
            return -math.pi / 2
        return ((self.position / self.sequence_length) * math.tau) - (math.pi / 2)


@dataclass(frozen=True)
class CloningVisualDigestFragment:
    index: int
    start: int
    end: int
    sequence_length: int
    wraps_origin: bool = False
    boundary_labels: tuple[str, ...] = field(default_factory=tuple)

    @property
    def length(self) -> int:
        if self.wraps_origin:
            return max(0, self.sequence_length - self.start + self.end)
        return max(0, self.end - self.start)

    @property
    def display_start(self) -> int:
        return self.start + 1

    @property
    def display_end(self) -> int:
        return self.end

    @property
    def coordinate_label(self) -> str:
        if self.wraps_origin:
            return f"{self.display_start}-{self.sequence_length} + 1-{self.display_end}"
        return f"{self.display_start}-{self.display_end}"

    @property
    def label(self) -> str:
        return f"Fragment {self.index}: {self.coordinate_label} ({self.length} bp)"


@dataclass(frozen=True)
class CloningVisualEnzymeSummary:
    enzyme_name: str
    site_sequence: str
    vector_cut_count: int = 0
    insert_cut_count: int = 0

    @property
    def has_any_cut(self) -> bool:
        return bool(self.vector_cut_count or self.insert_cut_count)


@dataclass(frozen=True)
class CloningVisualSequenceMap:
    role: str
    asset_name: str
    sequence_length: int
    is_circular: bool = False
    map_shape: str = "linear"
    regions: tuple[CloningVisualRegion, ...] = field(default_factory=tuple)
    restriction_sites: tuple[CloningVisualRestrictionSite, ...] = field(default_factory=tuple)
    digest_fragments: tuple[CloningVisualDigestFragment, ...] = field(default_factory=tuple)

    @property
    def is_circular_map(self) -> bool:
        return self.map_shape == "circular"

    @property
    def map_shape_label(self) -> str:
        return "circular map" if self.is_circular_map else "linear map"

    @property
    def source_topology_label(self) -> str:
        return "circular source" if self.is_circular else "linear source"


@dataclass(frozen=True)
class CloningAssemblyVisualPreview:
    vector_map: CloningVisualSequenceMap
    insert_map: CloningVisualSequenceMap
    selected_enzyme_names: tuple[str, ...] = field(default_factory=tuple)
    enzyme_summaries: tuple[CloningVisualEnzymeSummary, ...] = field(default_factory=tuple)
    selected_left_enzyme: str = ""
    selected_right_enzyme: str = ""
    helper_text: str = "Select enzymes to overlay cut sites and digest fragments. Left and right enzymes still define the cloning strategy."


@dataclass(frozen=True)
class CloningSequencePreviewPart:
    kind: str
    text: str
    label: str = ""


@dataclass(frozen=True)
class CloningDigestSequenceView:
    role: str
    asset_name: str
    sequence_length: int
    used_segments: tuple[CloningDigestSegment, ...] = field(default_factory=tuple)
    cut_markers: tuple[CloningRestrictionCutMarker, ...] = field(default_factory=tuple)
    double_strand_cut_views: tuple[CloningDoubleStrandCutView, ...] = field(default_factory=tuple)
    fragment_options: tuple[CloningDigestFragmentOption, ...] = field(default_factory=tuple)
    preview_parts: tuple[CloningSequencePreviewPart, ...] = field(default_factory=tuple)


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
    is_circular: bool
    assembled_sequence: str
    assembled_length: int
    is_valid: bool
    validation_messages: tuple[str, ...] = field(default_factory=tuple)
    cut_site_previews: tuple[CloningCutSitePreview, ...] = field(default_factory=tuple)
    digest_sequence_views: tuple[CloningDigestSequenceView, ...] = field(default_factory=tuple)
    vector_fragment_index: Optional[int] = None
    insert_fragment_index: Optional[int] = None
    vector_fragment_start: Optional[int] = None
    vector_fragment_end: Optional[int] = None
    insert_fragment_start: Optional[int] = None
    insert_fragment_end: Optional[int] = None


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
DNA_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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
    normalized_template_stem = Path(normalized_template_name).stem.lower()

    template_directory = _template_directory()
    if not template_directory.exists():
        raise ValueError("The template media directory is not available.")

    for candidate in template_directory.iterdir():
        if not candidate.is_file():
            continue
        if (
            candidate.name.lower() == normalized_template_name.lower()
            or candidate.stem.lower() == normalized_template_name.lower()
            or candidate.stem.lower() == normalized_template_stem
        ):
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
                template_name=template_name,
                record_id=str(record.id),
            )
            choices.append(
                (
                    asset_choice.encoded_value,
                    (
                        f"Template | {Path(template_name).stem} | "
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
    if not site:
        return []
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


def _build_cut_markers(*, sequence: str, enzyme_names: tuple[str, ...]) -> tuple[CloningRestrictionCutMarker, ...]:
    markers = []
    for enzyme_name in enzyme_names:
        enzyme = _get_enzyme_by_name(enzyme_name)
        if enzyme is None:
            continue
        markers.extend(
            CloningRestrictionCutMarker(
                enzyme_name=enzyme_name,
                position=position,
            )
            for position in _find_cut_positions(sequence, enzyme)
        )
    return tuple(sorted(markers, key=lambda marker: (marker.position, marker.enzyme_name)))


def _insert_cut_boundary(line: str, *, boundary: int, window_start: int, window_end: int) -> str:
    boundary = max(window_start, min(window_end, boundary))
    offset = boundary - window_start
    return f"{line[:offset]}|{line[offset:]}"


def _build_double_strand_cut_view(
    *,
    sequence: str,
    enzyme_name: str,
    enzyme,
    top_cut_position: int,
    context: int = 14,
) -> CloningDoubleStrandCutView:
    sequence = str(sequence or "").upper()
    site_sequence = str(getattr(enzyme, "site", "") or "")
    site_length = len(site_sequence)
    cut_offset = int(getattr(enzyme, "fst5", 0))
    reverse_cut_offset = int(getattr(enzyme, "fst3", cut_offset))
    recognition_start = max(0, top_cut_position - cut_offset)
    recognition_end = min(len(sequence), recognition_start + site_length)
    bottom_cut_position = recognition_start + site_length + reverse_cut_offset
    bottom_cut_position = max(0, min(len(sequence), bottom_cut_position))
    top_cut_position = max(0, min(len(sequence), top_cut_position))

    overhang_start = min(top_cut_position, bottom_cut_position)
    overhang_end = max(top_cut_position, bottom_cut_position)
    overhang_sequence = sequence[overhang_start:overhang_end]
    if top_cut_position == bottom_cut_position:
        overhang_type = "Blunt end"
    elif top_cut_position < bottom_cut_position:
        overhang_type = "5' overhang"
    else:
        overhang_type = "3' overhang"

    window_start = max(
        0,
        min(top_cut_position, bottom_cut_position, recognition_start) - context,
    )
    window_end = min(
        len(sequence),
        max(top_cut_position, bottom_cut_position, recognition_end) + context,
    )
    top_window = sequence[window_start:window_end]
    bottom_window = top_window.translate(DNA_COMPLEMENT)

    return CloningDoubleStrandCutView(
        enzyme_name=enzyme_name,
        top_cut_position=top_cut_position,
        bottom_cut_position=bottom_cut_position,
        top_line=_insert_cut_boundary(
            top_window,
            boundary=top_cut_position,
            window_start=window_start,
            window_end=window_end,
        ),
        bottom_line=_insert_cut_boundary(
            bottom_window,
            boundary=bottom_cut_position,
            window_start=window_start,
            window_end=window_end,
        ),
        overhang_sequence=overhang_sequence,
        overhang_type=overhang_type,
    )


def _build_double_strand_cut_views(
    *,
    sequence: str,
    enzyme_names: tuple[str, ...],
) -> tuple[CloningDoubleStrandCutView, ...]:
    cut_views = []
    for enzyme_name in enzyme_names:
        enzyme = _get_enzyme_by_name(enzyme_name)
        if enzyme is None:
            continue
        for position in _find_cut_positions(sequence, enzyme):
            cut_views.append(
                _build_double_strand_cut_view(
                    sequence=sequence,
                    enzyme_name=enzyme_name,
                    enzyme=enzyme,
                    top_cut_position=position,
                )
            )
    return tuple(sorted(cut_views, key=lambda view: (view.top_cut_position, view.enzyme_name)))


def _build_fragment_options(
    *,
    sequence: str,
    enzyme_name: str,
    field_name: str,
    selected_fragment_index: Optional[int],
    is_circular: bool = False,
) -> tuple[CloningDigestFragmentOption, ...]:
    enzyme = _get_enzyme_by_name(enzyme_name)
    if enzyme is None:
        return tuple()
    return tuple(
        CloningDigestFragmentOption(
            field_name=field_name,
            fragment_index=fragment.index,
            start=fragment.start,
            end=fragment.end,
            sequence_length=len(sequence),
            selected=selected_fragment_index == fragment.index,
            wraps_origin=fragment.wraps_origin,
        )
        for fragment in _digest_sequence_fragments(sequence, enzyme, is_circular=is_circular)
    )


def _dedupe_enzyme_names(enzyme_names) -> tuple[str, ...]:
    deduped = []
    for enzyme_name in enzyme_names or ():
        normalized = str(enzyme_name or "").strip()
        if normalized and normalized not in deduped:
            deduped.append(normalized)
    return tuple(deduped)


def _visual_regions_from_fragment_options(
    fragment_options: tuple[CloningDigestFragmentOption, ...],
) -> tuple[CloningVisualRegion, ...]:
    return tuple(
        CloningVisualRegion(
            label=f"Fragment {fragment.fragment_index}",
            field_name=fragment.field_name,
            fragment_index=fragment.fragment_index,
            start=fragment.start,
            end=fragment.end,
            sequence_length=fragment.sequence_length,
            selected=fragment.selected,
            selectable=True,
            wraps_origin=fragment.wraps_origin,
        )
        for fragment in fragment_options
    )


def _visual_regions_from_used_segments(
    *,
    used_segments: tuple[CloningDigestSegment, ...],
    sequence_length: int,
    fallback_label: str,
) -> tuple[CloningVisualRegion, ...]:
    if used_segments:
        return tuple(
            CloningVisualRegion(
                label=segment.label,
                field_name="",
                fragment_index=None,
                start=segment.start,
                end=segment.end,
                sequence_length=sequence_length,
                selected=True,
                selectable=False,
            )
            for segment in used_segments
        )
    if sequence_length <= 0:
        return tuple()
    return (
        CloningVisualRegion(
            label=fallback_label,
            field_name="",
            fragment_index=None,
            start=0,
            end=sequence_length,
            sequence_length=sequence_length,
            selected=False,
            selectable=False,
        ),
    )


def _build_visual_restriction_sites(
    *,
    role: str,
    sequence: str,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_names: tuple[str, ...],
    selected_left_enzyme: str,
    selected_right_enzyme: str,
) -> tuple[CloningVisualRestrictionSite, ...]:
    sequence = str(sequence or "").upper()
    vector_sequence = str(vector_sequence or "").upper()
    insert_sequence = str(insert_sequence or "").upper()
    sites = []
    for enzyme_name in enzyme_names:
        enzyme = _get_enzyme_by_name(enzyme_name)
        if enzyme is None:
            continue
        site_sequence = str(getattr(enzyme, "site", "") or "")
        vector_cut_positions = tuple(sorted(_find_cut_positions(vector_sequence, enzyme)))
        insert_site_count = len(_find_site_positions(insert_sequence, site_sequence))
        for position in sorted(_find_cut_positions(sequence, enzyme)):
            sites.append(
                CloningVisualRestrictionSite(
                    enzyme_name=enzyme_name,
                    site_sequence=site_sequence,
                    role=role,
                    position=position,
                    sequence_length=len(sequence),
                    vector_cut_count=len(vector_cut_positions),
                    insert_site_count=insert_site_count,
                    selected_left=enzyme_name == selected_left_enzyme,
                    selected_right=enzyme_name == selected_right_enzyme,
                    compatible=bool(vector_cut_positions),
                )
            )
    return tuple(sorted(sites, key=lambda site: (site.position, site.enzyme_name)))


def _cut_positions_for_enzyme_names(sequence: str, enzyme_names: tuple[str, ...]) -> dict[str, tuple[int, ...]]:
    cut_positions = {}
    for enzyme_name in enzyme_names:
        enzyme = _get_enzyme_by_name(enzyme_name)
        if enzyme is None:
            cut_positions[enzyme_name] = tuple()
            continue
        cut_positions[enzyme_name] = tuple(sorted(_find_cut_positions(sequence, enzyme)))
    return cut_positions


def _boundary_labels_by_position(cut_positions_by_enzyme: dict[str, tuple[int, ...]]) -> dict[int, tuple[str, ...]]:
    labels_by_position = {}
    for enzyme_name, positions in cut_positions_by_enzyme.items():
        for position in positions:
            labels_by_position.setdefault(position, []).append(enzyme_name)
    return {
        position: tuple(sorted(enzyme_names))
        for position, enzyme_names in labels_by_position.items()
    }


def _build_visual_digest_fragments(
    *,
    sequence: str,
    enzyme_names: tuple[str, ...],
    is_circular: bool,
) -> tuple[CloningVisualDigestFragment, ...]:
    sequence = str(sequence or "").upper()
    sequence_length = len(sequence)
    if sequence_length <= 0 or not enzyme_names:
        return tuple()

    cut_positions_by_enzyme = _cut_positions_for_enzyme_names(sequence, enzyme_names)
    labels_by_position = _boundary_labels_by_position(cut_positions_by_enzyme)
    cut_positions = sorted(labels_by_position)
    if not cut_positions:
        return tuple()

    if is_circular:
        fragments = []
        for index, start in enumerate(cut_positions, start=1):
            end = cut_positions[index % len(cut_positions)]
            wraps_origin = end <= start
            fragments.append(
                CloningVisualDigestFragment(
                    index=index,
                    start=start,
                    end=end,
                    sequence_length=sequence_length,
                    wraps_origin=wraps_origin,
                    boundary_labels=tuple(
                        sorted(set(labels_by_position.get(start, ())) | set(labels_by_position.get(end, ())))
                    ),
                )
            )
        return tuple(fragments)

    boundaries = [0, *cut_positions, sequence_length]
    fragments = []
    for index, (start, end) in enumerate(zip(boundaries, boundaries[1:]), start=1):
        if end <= start:
            continue
        fragments.append(
            CloningVisualDigestFragment(
                index=len(fragments) + 1,
                start=start,
                end=end,
                sequence_length=sequence_length,
                boundary_labels=tuple(
                    sorted(set(labels_by_position.get(start, ())) | set(labels_by_position.get(end, ())))
                ),
            )
        )
    return tuple(fragments)


def _build_visual_enzyme_summaries(
    *,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_names: tuple[str, ...],
) -> tuple[CloningVisualEnzymeSummary, ...]:
    summaries = []
    for enzyme_name in enzyme_names:
        enzyme = _get_enzyme_by_name(enzyme_name)
        if enzyme is None:
            summaries.append(CloningVisualEnzymeSummary(enzyme_name=enzyme_name, site_sequence=""))
            continue
        summaries.append(
            CloningVisualEnzymeSummary(
                enzyme_name=enzyme_name,
                site_sequence=str(getattr(enzyme, "site", "") or ""),
                vector_cut_count=len(_find_cut_positions(vector_sequence, enzyme)),
                insert_cut_count=len(_find_cut_positions(insert_sequence, enzyme)),
            )
        )
    return tuple(summaries)


def _visual_map_shape_for_topology(*, is_circular: bool) -> str:
    return "circular" if is_circular else "linear"


def _visual_overlay_enzyme_names(
    *,
    map_enzyme_names=(),
    selected_left_enzyme: str = "",
    selected_right_enzyme: str = "",
) -> tuple[str, ...]:
    overlay_enzyme_names = _dedupe_enzyme_names(map_enzyme_names)
    if overlay_enzyme_names:
        return overlay_enzyme_names
    return _dedupe_enzyme_names((selected_left_enzyme, selected_right_enzyme))


def build_cloning_assembly_visual_preview(
    *,
    vector_asset: ResolvedCloningAsset,
    insert_asset: ResolvedCloningAsset,
    selected_left_enzyme: str = "",
    selected_right_enzyme: str = "",
    map_enzyme_names=(),
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
) -> CloningAssemblyVisualPreview:
    try:
        vector_fragment_index = _normalize_fragment_index(vector_fragment_index)
    except ValueError:
        vector_fragment_index = None
    try:
        insert_fragment_index = _normalize_fragment_index(insert_fragment_index)
    except ValueError:
        insert_fragment_index = None
    selected_left_enzyme = str(selected_left_enzyme or "").strip()
    selected_right_enzyme = str(selected_right_enzyme or "").strip()
    overlay_enzyme_names = _visual_overlay_enzyme_names(
        map_enzyme_names=map_enzyme_names,
        selected_left_enzyme=selected_left_enzyme,
        selected_right_enzyme=selected_right_enzyme,
    )

    validation_result = CloningValidationResult(is_valid=False)
    if selected_left_enzyme and selected_right_enzyme:
        validation_result = _validate_restriction_ligation(
            vector_asset.sequence,
            insert_asset.sequence,
            selected_left_enzyme,
            selected_right_enzyme,
            vector_fragment_index=vector_fragment_index,
            insert_fragment_index=insert_fragment_index,
            vector_is_circular=vector_asset.is_circular,
            insert_is_circular=insert_asset.is_circular,
        )
    vector_used_segments, insert_used_segments = _infer_used_digest_segments(
        vector_sequence=vector_asset.sequence,
        insert_sequence=insert_asset.sequence,
        left_enzyme_name=selected_left_enzyme,
        right_enzyme_name=selected_right_enzyme,
        validation_result=validation_result,
        explicit_fragment_selection=vector_fragment_index is not None or insert_fragment_index is not None,
    )

    same_enzyme_name = (
        selected_left_enzyme
        if selected_left_enzyme and selected_left_enzyme == selected_right_enzyme
        else None
    )
    vector_fragment_options = (
        _build_fragment_options(
            sequence=vector_asset.sequence,
            enzyme_name=same_enzyme_name,
            field_name="vector_fragment_index",
            selected_fragment_index=vector_fragment_index,
            is_circular=vector_asset.is_circular,
        )
        if same_enzyme_name
        else tuple()
    )
    insert_fragment_options = (
        _build_fragment_options(
            sequence=insert_asset.sequence,
            enzyme_name=same_enzyme_name,
            field_name="insert_fragment_index",
            selected_fragment_index=insert_fragment_index,
            is_circular=insert_asset.is_circular,
        )
        if same_enzyme_name
        else tuple()
    )

    vector_regions = _visual_regions_from_fragment_options(vector_fragment_options)
    if not vector_regions:
        vector_regions = _visual_regions_from_used_segments(
            used_segments=vector_used_segments,
            sequence_length=len(vector_asset.sequence),
            fallback_label="Full vector",
        )
    insert_regions = _visual_regions_from_fragment_options(insert_fragment_options)
    if not insert_regions:
        insert_regions = _visual_regions_from_used_segments(
            used_segments=insert_used_segments,
            sequence_length=len(insert_asset.sequence),
            fallback_label="Full insert",
        )

    return CloningAssemblyVisualPreview(
        selected_enzyme_names=overlay_enzyme_names,
        enzyme_summaries=_build_visual_enzyme_summaries(
            vector_sequence=vector_asset.sequence,
            insert_sequence=insert_asset.sequence,
            enzyme_names=overlay_enzyme_names,
        ),
        selected_left_enzyme=selected_left_enzyme,
        selected_right_enzyme=selected_right_enzyme,
        vector_map=CloningVisualSequenceMap(
            role="Vector",
            asset_name=vector_asset.name,
            sequence_length=len(vector_asset.sequence),
            is_circular=vector_asset.is_circular,
            map_shape=_visual_map_shape_for_topology(is_circular=vector_asset.is_circular),
            regions=vector_regions,
            restriction_sites=_build_visual_restriction_sites(
                role="vector",
                sequence=vector_asset.sequence,
                vector_sequence=vector_asset.sequence,
                insert_sequence=insert_asset.sequence,
                enzyme_names=overlay_enzyme_names,
                selected_left_enzyme=selected_left_enzyme,
                selected_right_enzyme=selected_right_enzyme,
            ),
            digest_fragments=_build_visual_digest_fragments(
                sequence=vector_asset.sequence,
                enzyme_names=overlay_enzyme_names,
                is_circular=vector_asset.is_circular,
            ),
        ),
        insert_map=CloningVisualSequenceMap(
            role="Insert",
            asset_name=insert_asset.name,
            sequence_length=len(insert_asset.sequence),
            is_circular=insert_asset.is_circular,
            map_shape=_visual_map_shape_for_topology(is_circular=insert_asset.is_circular),
            regions=insert_regions,
            restriction_sites=_build_visual_restriction_sites(
                role="insert",
                sequence=insert_asset.sequence,
                vector_sequence=vector_asset.sequence,
                insert_sequence=insert_asset.sequence,
                enzyme_names=overlay_enzyme_names,
                selected_left_enzyme=selected_left_enzyme,
                selected_right_enzyme=selected_right_enzyme,
            ),
            digest_fragments=_build_visual_digest_fragments(
                sequence=insert_asset.sequence,
                enzyme_names=overlay_enzyme_names,
                is_circular=insert_asset.is_circular,
            ),
        ),
    )


def _non_empty_segment(start: int, end: int, label: str) -> Optional[CloningDigestSegment]:
    if end <= start:
        return None
    return CloningDigestSegment(start=start, end=end, label=label)


def _compact_segments(segments: tuple[Optional[CloningDigestSegment], ...]) -> tuple[CloningDigestSegment, ...]:
    return tuple(segment for segment in segments if segment is not None)


def _infer_used_digest_segments(
    *,
    vector_sequence: str,
    insert_sequence: str,
    left_enzyme_name: str,
    right_enzyme_name: str,
    validation_result: CloningValidationResult,
    explicit_fragment_selection: bool = False,
) -> tuple[tuple[CloningDigestSegment, ...], tuple[CloningDigestSegment, ...]]:
    if not validation_result.is_valid:
        return tuple(), tuple()

    if validation_result.vector_fragment_segments or validation_result.insert_fragment_segments:
        vector_label = "Selected vector fragment" if explicit_fragment_selection else "Vector backbone"
        insert_label = "Selected insert fragment" if explicit_fragment_selection else "Insert fragment"
        return (
            tuple(
                CloningDigestSegment(
                    start=start,
                    end=end,
                    label=vector_label,
                )
                for start, end in validation_result.vector_fragment_segments
            ),
            tuple(
                CloningDigestSegment(
                    start=start,
                    end=end,
                    label=insert_label,
                )
                for start, end in validation_result.insert_fragment_segments
            ),
        )

    fragment_positions = (
        validation_result.vector_fragment_start,
        validation_result.vector_fragment_end,
        validation_result.insert_fragment_start,
        validation_result.insert_fragment_end,
    )
    if all(position is not None for position in fragment_positions):
        return (
            _compact_segments(
                (
                    _non_empty_segment(
                        validation_result.vector_fragment_start,
                        validation_result.vector_fragment_end,
                        "Selected vector fragment",
                    ),
                )
            ),
            _compact_segments(
                (
                    _non_empty_segment(
                        validation_result.insert_fragment_start,
                        validation_result.insert_fragment_end,
                        "Selected insert fragment",
                    ),
                )
            ),
        )

    if left_enzyme_name == right_enzyme_name:
        enzyme = _get_enzyme_by_name(left_enzyme_name)
        if enzyme is None:
            return tuple(), tuple()
        vector_cuts = sorted(_find_cut_positions(vector_sequence, enzyme))
        insert_cuts = sorted(_find_cut_positions(insert_sequence, enzyme))

        if len(vector_cuts) == 1 and len(insert_cuts) == 0:
            return (
                _compact_segments((_non_empty_segment(0, len(vector_sequence), "Linearized vector"),)),
                _compact_segments((_non_empty_segment(0, len(insert_sequence), "Full insert"),)),
            )
        if len(vector_cuts) == 2 and len(insert_cuts) == 2:
            vector_left_cut, vector_right_cut = vector_cuts
            insert_left_cut, insert_right_cut = insert_cuts
            return (
                _compact_segments(
                    (
                        _non_empty_segment(0, vector_left_cut, "Vector left flank"),
                        _non_empty_segment(vector_right_cut, len(vector_sequence), "Vector right flank"),
                    )
                ),
                _compact_segments(
                    (
                        _non_empty_segment(insert_left_cut, insert_right_cut, "Excised insert fragment"),
                    )
                ),
            )
        return tuple(), tuple()

    left_enzyme = _get_enzyme_by_name(left_enzyme_name)
    right_enzyme = _get_enzyme_by_name(right_enzyme_name)
    if left_enzyme is None or right_enzyme is None:
        return tuple(), tuple()
    left_hits = sorted(_find_cut_positions(vector_sequence, left_enzyme))
    right_hits = sorted(_find_cut_positions(vector_sequence, right_enzyme))
    if len(left_hits) != 1 or len(right_hits) != 1:
        return tuple(), tuple()
    left_cut = left_hits[0]
    right_cut = right_hits[0]
    if right_cut <= left_cut:
        return tuple(), tuple()
    return (
        _compact_segments(
            (
                _non_empty_segment(0, left_cut, "Vector left flank"),
                _non_empty_segment(right_cut, len(vector_sequence), "Vector right flank"),
            )
        ),
        _compact_segments((_non_empty_segment(0, len(insert_sequence), "Full insert"),)),
    )


def _range_overlaps_used_segment(
    *,
    start: int,
    end: int,
    used_segments: tuple[CloningDigestSegment, ...],
) -> bool:
    return any(segment.start < end and start < segment.end for segment in used_segments)


def _merge_preview_intervals(
    intervals: list[tuple[int, int]],
    *,
    sequence_length: int,
    merge_gap: int = 12,
) -> list[tuple[int, int]]:
    normalized = sorted(
        (max(0, start), min(sequence_length, end))
        for start, end in intervals
        if end > start
    )
    if not normalized:
        return [(0, min(sequence_length, 80))] if sequence_length else []

    merged = []
    for start, end in normalized:
        if not merged or start > merged[-1][1] + merge_gap:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def _build_sequence_preview_parts(
    *,
    sequence: str,
    used_segments: tuple[CloningDigestSegment, ...],
    cut_markers: tuple[CloningRestrictionCutMarker, ...],
    cut_context: int = 18,
    fragment_edge_context: int = 28,
    max_full_fragment_length: int = 96,
) -> tuple[CloningSequencePreviewPart, ...]:
    sequence = str(sequence or "")
    sequence_length = len(sequence)
    intervals = []
    for marker in cut_markers:
        intervals.append((marker.position - cut_context, marker.position + cut_context))
    for segment in used_segments:
        if segment.length <= max_full_fragment_length:
            intervals.append((segment.start, segment.end))
        else:
            intervals.append((segment.start, segment.start + fragment_edge_context))
            intervals.append((segment.end - fragment_edge_context, segment.end))

    merged_intervals = _merge_preview_intervals(intervals, sequence_length=sequence_length)
    markers_by_position = {}
    for marker in cut_markers:
        markers_by_position.setdefault(marker.position, []).append(marker)

    parts = []
    previous_end = None
    rendered_marker_positions = set()
    for start, end in merged_intervals:
        if previous_end is not None and start > previous_end:
            parts.append(CloningSequencePreviewPart(kind="omitted", text="..."))

        boundaries = {start, end}
        boundaries.update(
            marker.position
            for marker in cut_markers
            if start <= marker.position <= end
        )
        for segment in used_segments:
            if start < segment.start < end:
                boundaries.add(segment.start)
            if start < segment.end < end:
                boundaries.add(segment.end)
        ordered_boundaries = sorted(boundaries)

        for index, boundary_start in enumerate(ordered_boundaries[:-1]):
            for marker in markers_by_position.get(boundary_start, []):
                marker_key = (marker.position, marker.enzyme_name)
                if marker_key not in rendered_marker_positions:
                    parts.append(
                        CloningSequencePreviewPart(kind="cut", text="|", label=marker.label)
                    )
                    rendered_marker_positions.add(marker_key)

            boundary_end = ordered_boundaries[index + 1]
            if boundary_end <= boundary_start:
                continue
            parts.append(
                CloningSequencePreviewPart(
                    kind=(
                        "used"
                        if _range_overlaps_used_segment(
                            start=boundary_start,
                            end=boundary_end,
                            used_segments=used_segments,
                        )
                        else "flank"
                    ),
                    text=sequence[boundary_start:boundary_end],
                )
            )

        for marker in markers_by_position.get(end, []):
            marker_key = (marker.position, marker.enzyme_name)
            if marker_key not in rendered_marker_positions:
                parts.append(
                    CloningSequencePreviewPart(kind="cut", text="|", label=marker.label)
                )
                rendered_marker_positions.add(marker_key)
        previous_end = end

    return tuple(parts)


def _build_digest_sequence_views(
    *,
    vector_asset: ResolvedCloningAsset,
    insert_asset: ResolvedCloningAsset,
    left_enzyme_name: str,
    right_enzyme_name: str,
    validation_result: CloningValidationResult,
    vector_fragment_index: Optional[int],
    insert_fragment_index: Optional[int],
) -> tuple[CloningDigestSequenceView, ...]:
    enzyme_names = tuple(
        enzyme_name
        for index, enzyme_name in enumerate((left_enzyme_name, right_enzyme_name))
        if enzyme_name and enzyme_name not in (left_enzyme_name, right_enzyme_name)[:index]
    )
    vector_used_segments, insert_used_segments = _infer_used_digest_segments(
        vector_sequence=vector_asset.sequence,
        insert_sequence=insert_asset.sequence,
        left_enzyme_name=left_enzyme_name,
        right_enzyme_name=right_enzyme_name,
        validation_result=validation_result,
        explicit_fragment_selection=vector_fragment_index is not None or insert_fragment_index is not None,
    )
    vector_cut_markers = _build_cut_markers(
        sequence=vector_asset.sequence,
        enzyme_names=enzyme_names,
    )
    insert_cut_markers = _build_cut_markers(
        sequence=insert_asset.sequence,
        enzyme_names=enzyme_names,
    )
    same_enzyme_name = left_enzyme_name if left_enzyme_name == right_enzyme_name else None
    return (
        CloningDigestSequenceView(
            role="Plasmid / vector",
            asset_name=vector_asset.name,
            sequence_length=len(vector_asset.sequence),
            used_segments=vector_used_segments,
            cut_markers=vector_cut_markers,
            double_strand_cut_views=_build_double_strand_cut_views(
                sequence=vector_asset.sequence,
                enzyme_names=enzyme_names,
            ),
            fragment_options=(
                _build_fragment_options(
                    sequence=vector_asset.sequence,
                    enzyme_name=same_enzyme_name,
                    field_name="vector_fragment_index",
                    selected_fragment_index=vector_fragment_index,
                    is_circular=vector_asset.is_circular,
                )
                if same_enzyme_name
                else tuple()
            ),
            preview_parts=_build_sequence_preview_parts(
                sequence=vector_asset.sequence,
                used_segments=vector_used_segments,
                cut_markers=vector_cut_markers,
            ),
        ),
        CloningDigestSequenceView(
            role="Insert",
            asset_name=insert_asset.name,
            sequence_length=len(insert_asset.sequence),
            used_segments=insert_used_segments,
            cut_markers=insert_cut_markers,
            double_strand_cut_views=_build_double_strand_cut_views(
                sequence=insert_asset.sequence,
                enzyme_names=enzyme_names,
            ),
            fragment_options=(
                _build_fragment_options(
                    sequence=insert_asset.sequence,
                    enzyme_name=same_enzyme_name,
                    field_name="insert_fragment_index",
                    selected_fragment_index=insert_fragment_index,
                    is_circular=insert_asset.is_circular,
                )
                if same_enzyme_name
                else tuple()
            ),
            preview_parts=_build_sequence_preview_parts(
                sequence=insert_asset.sequence,
                used_segments=insert_used_segments,
                cut_markers=insert_cut_markers,
            ),
        ),
    )


@dataclass(frozen=True)
class DigestedFragment:
    index: int
    start: int
    end: int
    sequence: str
    wraps_origin: bool = False

    @property
    def length(self) -> int:
        return len(self.sequence)

    def source_segments(self, sequence_length: int) -> tuple[tuple[int, int], ...]:
        if self.wraps_origin:
            return _compact_source_segments(((self.start, sequence_length), (0, self.end)))
        return _compact_source_segments(((self.start, self.end),))


def _compact_source_segments(segments: tuple[tuple[int, int], ...]) -> tuple[tuple[int, int], ...]:
    return tuple((start, end) for start, end in segments if end > start)


def _sequence_from_segments(sequence: str, segments: tuple[tuple[int, int], ...]) -> str:
    return "".join(sequence[start:end] for start, end in segments)


def _offset_segments(
    segments: tuple[tuple[int, int], ...],
    *,
    assembled_offset: int,
) -> tuple[tuple[int, int, int], ...]:
    mapped_segments = []
    offset = assembled_offset
    for start, end in segments:
        if end <= start:
            continue
        mapped_segments.append((start, end, offset))
        offset += end - start
    return tuple(mapped_segments)


def _digest_sequence_fragments(sequence: str, enzyme, *, is_circular: bool = False) -> tuple[DigestedFragment, ...]:
    cut_positions = sorted(
        set(position for position in _find_cut_positions(sequence, enzyme) if 0 < position < len(sequence))
    )
    if is_circular and cut_positions:
        fragments = []
        for index, start in enumerate(cut_positions, start=1):
            end = cut_positions[index % len(cut_positions)]
            wraps_origin = end <= start
            segments = (
                ((start, len(sequence)), (0, end))
                if wraps_origin
                else ((start, end),)
            )
            fragment_sequence = _sequence_from_segments(sequence, _compact_source_segments(segments))
            if not fragment_sequence:
                continue
            fragments.append(
                DigestedFragment(
                    index=index,
                    start=start,
                    end=end,
                    sequence=fragment_sequence,
                    wraps_origin=wraps_origin,
                )
            )
        return tuple(fragments)

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
                wraps_origin=False,
            )
        )
    return tuple(fragments)


def _select_digest_fragment(
    *,
    sequence: str,
    enzyme,
    fragment_index: Optional[int],
    is_circular: bool = False,
    minimum_length: int = 1,
    label: str,
) -> DigestedFragment:
    fragments = _digest_sequence_fragments(sequence, enzyme, is_circular=is_circular)
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


def build_digest_fragment_choices(*, sequence: str, enzyme_name: str, minimum_length: int = 1, is_circular: bool = False):
    enzyme = _get_enzyme_by_name(enzyme_name)
    if enzyme is None:
        return []
    fragments = _digest_sequence_fragments(sequence, enzyme, is_circular=is_circular)
    eligible_fragments = sorted(
        (fragment for fragment in fragments if fragment.length >= minimum_length),
        key=lambda fragment: (-fragment.length, fragment.start, fragment.index),
    )
    return [
        (
            str(fragment.index),
            (
                f"Fragment {fragment.index} | {fragment.length} bp | bases "
                f"{fragment.start + 1}-{len(sequence)} + 1-{fragment.end}"
                if fragment.wraps_origin
                else f"Fragment {fragment.index} | {fragment.length} bp | bases {fragment.start + 1}-{fragment.end}"
            ),
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
        raise ValueError(f"Sequence file '{sequence_file.name}' could not be parsed for cloning.")
    return records


def _record_is_circular(record) -> bool:
    topology = str(getattr(record, "annotations", {}).get("topology", "") or "").strip().lower()
    return topology == "circular"


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
        is_circular=_record_is_circular(record),
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
        is_circular=_record_is_circular(record),
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
        is_circular=False,
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
        template_name=str(template_name or template_path.name).strip(),
        record_id=str(record.id),
        name=template_path.name,
        sequence=str(record.seq).upper(),
        is_circular=_record_is_circular(record),
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


@dataclass(frozen=True)
class RestrictionLigationStrategy:
    left_enzyme_name: str
    right_enzyme_name: str
    vector_fragment_index: Optional[int] = None
    insert_fragment_index: Optional[int] = None
    vector_is_circular: bool = False
    insert_is_circular: bool = False

    def validate(self, *, vector_sequence: str, insert_sequence: str) -> CloningValidationResult:
        return _validate_restriction_ligation(
            vector_sequence,
            insert_sequence,
            self.left_enzyme_name,
            self.right_enzyme_name,
            vector_fragment_index=self.vector_fragment_index,
            insert_fragment_index=self.insert_fragment_index,
            vector_is_circular=self.vector_is_circular,
            insert_is_circular=self.insert_is_circular,
        )


def _build_cloning_strategy(
    *,
    assembly_strategy: str,
    left_enzyme: str,
    right_enzyme: str,
    vector_fragment_index: Optional[int],
    insert_fragment_index: Optional[int],
    vector_is_circular: bool = False,
    insert_is_circular: bool = False,
) -> RestrictionLigationStrategy:
    if assembly_strategy != CloningConstruct.STRATEGY_RESTRICTION_LIGATION:
        raise ValueError("Only restriction ligation is currently supported.")

    return RestrictionLigationStrategy(
        left_enzyme_name=left_enzyme,
        right_enzyme_name=right_enzyme,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
        vector_is_circular=vector_is_circular,
        insert_is_circular=insert_is_circular,
    )


def _validate_same_enzyme_fragment_ligation(
    *,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_name: str,
    enzyme,
    vector_fragment_index: Optional[int],
    insert_fragment_index: Optional[int],
    vector_is_circular: bool = False,
    insert_is_circular: bool = False,
) -> CloningValidationResult:
    try:
        vector_fragment = _select_digest_fragment(
            sequence=vector_sequence,
            enzyme=enzyme,
            fragment_index=vector_fragment_index,
            is_circular=vector_is_circular,
            label="Vector",
        )
        insert_fragment = _select_digest_fragment(
            sequence=insert_sequence,
            enzyme=enzyme,
            fragment_index=insert_fragment_index,
            is_circular=insert_is_circular,
            label="Insert",
        )
    except ValueError as exc:
        return _invalid_cloning_result(str(exc))

    vector_fragment_segments = vector_fragment.source_segments(len(vector_sequence))
    insert_fragment_segments = insert_fragment.source_segments(len(insert_sequence))
    if vector_fragment.wraps_origin:
        left_vector_segments = _compact_source_segments(((0, vector_fragment.end),))
        right_vector_segments = _compact_source_segments(((vector_fragment.start, len(vector_sequence)),))
        insert_offset = sum(end - start for start, end in left_vector_segments)
        assembled_sequence = (
            _sequence_from_segments(vector_sequence, left_vector_segments)
            + insert_fragment.sequence
            + _sequence_from_segments(vector_sequence, right_vector_segments)
        )
        vector_assembly_segments = (
            *_offset_segments(left_vector_segments, assembled_offset=0),
            *_offset_segments(
                right_vector_segments,
                assembled_offset=insert_offset + len(insert_fragment.sequence),
            ),
        )
    else:
        insert_offset = len(vector_fragment.sequence)
        assembled_sequence = vector_fragment.sequence + insert_fragment.sequence
        vector_assembly_segments = _offset_segments(vector_fragment_segments, assembled_offset=0)

    return CloningValidationResult(
        is_valid=True,
        validation_messages=(
            f"Construct assembled successfully using {enzyme_name} fragment selection (vector fragment {vector_fragment.index} and insert fragment {insert_fragment.index}).",
        ),
        assembled_sequence=assembled_sequence,
        insertion_start=insert_offset,
        inserted_length=len(insert_fragment.sequence),
        vector_fragment_start=vector_fragment.start,
        vector_fragment_end=vector_fragment.end,
        insert_fragment_start=insert_fragment.start,
        insert_fragment_end=insert_fragment.end,
        vector_fragment_segments=vector_fragment_segments,
        insert_fragment_segments=insert_fragment_segments,
        vector_assembly_segments=vector_assembly_segments,
        insert_assembly_segments=_offset_segments(
            insert_fragment_segments,
            assembled_offset=insert_offset,
        ),
    )


def _validate_same_enzyme_ligation(
    *,
    vector_sequence: str,
    insert_sequence: str,
    enzyme_name: str,
    enzyme,
    vector_is_circular: bool = False,
    insert_is_circular: bool = False,
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
            vector_fragment_segments=((0, len(vector_sequence)),),
            insert_fragment_segments=((0, len(insert_sequence)),),
            vector_assembly_segments=(
                (0, insertion_point, 0),
                (insertion_point, len(vector_sequence), insertion_point + len(insert_sequence)),
            ),
            insert_assembly_segments=((0, len(insert_sequence), insertion_point),),
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
            vector_fragment_segments=_compact_source_segments(
                ((0, vector_left_cut), (vector_right_cut, len(vector_sequence)))
            ),
            insert_fragment_segments=((insert_left_cut, insert_right_cut),),
            vector_assembly_segments=(
                *_offset_segments(((0, vector_left_cut),), assembled_offset=0),
                *_offset_segments(
                    ((vector_right_cut, len(vector_sequence)),),
                    assembled_offset=vector_left_cut + len(insert_fragment),
                ),
            ),
            insert_assembly_segments=((insert_left_cut, insert_right_cut, vector_left_cut),),
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
    vector_is_circular: bool = False,
    insert_is_circular: bool = False,
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
        vector_fragment_segments=_compact_source_segments(
            ((0, left_cut), (right_cut, len(vector_sequence)))
        ),
        insert_fragment_segments=((0, len(insert_sequence)),),
        vector_assembly_segments=(
            *_offset_segments(((0, left_cut),), assembled_offset=0),
            *_offset_segments(
                ((right_cut, len(vector_sequence)),),
                assembled_offset=left_cut + len(insert_sequence),
            ),
        ),
        insert_assembly_segments=((0, len(insert_sequence), left_cut),),
    )


def _validate_restriction_ligation(
    vector_sequence: str,
    insert_sequence: str,
    left_enzyme_name: str,
    right_enzyme_name: str,
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
    vector_is_circular: bool = False,
    insert_is_circular: bool = False,
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
                vector_is_circular=vector_is_circular,
                insert_is_circular=insert_is_circular,
            )
        # Preserve the legacy same-enzyme topology for already-saved constructs that do not store fragment indexes.
        return _validate_same_enzyme_ligation(
            vector_sequence=vector_sequence,
            insert_sequence=insert_sequence,
            enzyme_name=left_enzyme_name,
            enzyme=left_enzyme,
            vector_is_circular=vector_is_circular,
            insert_is_circular=insert_is_circular,
        )

    return _validate_two_enzyme_ligation(
        vector_sequence=vector_sequence,
        insert_sequence=insert_sequence,
        left_enzyme_name=left_enzyme_name,
        right_enzyme_name=right_enzyme_name,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        vector_is_circular=vector_is_circular,
        insert_is_circular=insert_is_circular,
    )


def preview_cloning_construct(
    *,
    vector_asset: ResolvedCloningAsset,
    insert_asset: ResolvedCloningAsset,
    assembly_strategy: str,
    left_enzyme: str,
    right_enzyme: str,
    is_circular: Optional[bool] = None,
    vector_fragment_index: Optional[int] = None,
    insert_fragment_index: Optional[int] = None,
) -> CloningConstructPreview:
    vector_fragment_index = _normalize_fragment_index(vector_fragment_index)
    insert_fragment_index = _normalize_fragment_index(insert_fragment_index)
    result_is_circular = vector_asset.is_circular if is_circular is None else bool(is_circular)
    cloning_strategy = _build_cloning_strategy(
        assembly_strategy=assembly_strategy,
        left_enzyme=left_enzyme,
        right_enzyme=right_enzyme,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
        vector_is_circular=vector_asset.is_circular,
        insert_is_circular=insert_asset.is_circular,
    )

    validation_messages = []
    if vector_asset.message:
        validation_messages.append(vector_asset.message)
    if insert_asset.message:
        validation_messages.append(insert_asset.message)

    validation_result = cloning_strategy.validate(
        vector_sequence=vector_asset.sequence,
        insert_sequence=insert_asset.sequence,
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
        is_circular=result_is_circular,
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
        digest_sequence_views=_build_digest_sequence_views(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            left_enzyme_name=left_enzyme,
            right_enzyme_name=right_enzyme,
            validation_result=validation_result,
            vector_fragment_index=vector_fragment_index,
            insert_fragment_index=insert_fragment_index,
        ),
        vector_fragment_start=validation_result.vector_fragment_start,
        vector_fragment_end=validation_result.vector_fragment_end,
        insert_fragment_start=validation_result.insert_fragment_start,
        insert_fragment_end=validation_result.insert_fragment_end,
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
            vector_is_circular=preview_data.vector_asset.is_circular,
            insert_is_circular=preview_data.insert_asset.is_circular,
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
    detail_display_snapshot = getattr(construct, "detail_display_snapshot", None)
    if detail_display_snapshot:
        return _deserialize_detail_display(detail_display_snapshot)

    source_errors = []
    try:
        vector_asset = _resolve_construct_asset(
            source_type=construct.vector_source_type,
            sequence_file=construct.vector_sequence_file,
            pcr_product=construct.vector_pcr_product,
            template_name=getattr(construct, "vector_template_name", ""),
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
            template_name=getattr(construct, "insert_template_name", ""),
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
            vector_fragment_index=getattr(construct, "vector_fragment_index", None),
            insert_fragment_index=getattr(construct, "insert_fragment_index", None),
            vector_is_circular=vector_asset.is_circular,
            insert_is_circular=insert_asset.is_circular,
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
        is_circular=preview_data.is_circular,
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
    is_circular: Optional[bool] = None,
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
        is_circular=is_circular,
        vector_fragment_index=vector_fragment_index,
        insert_fragment_index=insert_fragment_index,
    )
    return save_cloning_construct(
        name=name,
        description=description,
        preview_data=preview_data,
        user=user,
    )
