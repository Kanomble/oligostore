from django.contrib.auth import authenticate, get_user_model
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.db import OperationalError, ProgrammingError
from django.contrib.auth.models import User
from django import forms
from .access import (
    accessible_pcr_products,
    accessible_primer_pairs,
    accessible_primers,
    accessible_sequence_files,
)
from .models import CloningConstruct, Primer, PrimerPair, Project, SequenceFile
import re
from .services.cloning import (
    build_digest_fragment_choices,
    build_pcr_product_asset_choice,
    build_sequence_file_asset_choice,
    get_detected_enzyme_choices,
    get_template_asset_choices,
    resolve_asset_choice,
)
from .services.sequence_loader import load_sequences
from .services.sequence_records import get_sequence_records

def apply_tailwind_classes(fields):
    for field in fields.values():
        if isinstance(field.widget, forms.Textarea):
            field.widget.attrs.update({"class": "textarea textarea-bordered w-full"})
        elif isinstance(field.widget, forms.Select):
            field.widget.attrs.update({"class": "select select-bordered w-full"})
        elif isinstance(field, (forms.CharField, forms.IntegerField, forms.FloatField)):
            field.widget.attrs.update({"class": "input input-bordered w-full"})

def clean_sequence_value(value, allow_n=True, max_length=None):
    seq = re.sub(r"\s+", "", value or "").upper()
    if max_length is not None and len(seq) > max_length:
        raise forms.ValidationError(
            f"Sequence must be {max_length} bases or fewer."
        )
    pattern = r"[ACGTN]+" if allow_n else r"[ACGT]+"
    if not re.fullmatch(pattern, seq):
        raise forms.ValidationError(
            "Sequence may only contain the characters A, C, G, T"
            + (" or N" if allow_n else "")
            + " (no spaces or numbers)."
        )
    return seq

def clean_optional_sequence_value(value, allow_n=True):
    seq = re.sub(r"\s+", "", value or "").upper()
    if not seq:
        return ""
    return clean_sequence_value(seq, allow_n=allow_n)

class ProjectForm(forms.ModelForm):
    class Meta:
        model = Project
        fields = ["name", "description"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        apply_tailwind_classes(self.fields)


def _build_cloning_asset_choices(user):
    if user is None:
        return []

    asset_choices = []
    for sequence_file in accessible_sequence_files(user).order_by("name"):
        try:
            records = get_sequence_records(sequence_file, load_sequences)
        except Exception:
            records = None
        if records:
            multi_record_label = (
                "multi-record file" if len(records) > 1 else "single-record file"
            )
            for record in records:
                asset_choice = build_sequence_file_asset_choice(
                    sequence_file_id=sequence_file.id,
                    record_id=str(record.id),
                )
                asset_choices.append(
                    (
                        asset_choice.encoded_value,
                        (
                            f"Sequence file | {sequence_file.name} | "
                            f"record {record.id} | {len(record.seq)} bp | "
                            f"{multi_record_label}"
                        ),
                    )
                )
        else:
            asset_choices.append(
                (
                    f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{sequence_file.id}",
                    f"Sequence file | {sequence_file.name} | record unavailable | length unavailable | parse required",
                )
            )
    for pcr_product in accessible_pcr_products(user).order_by("name"):
        asset_choice = build_pcr_product_asset_choice(pcr_product_id=pcr_product.id)
        asset_choices.append(
            (
                asset_choice.encoded_value,
                (
                    f"PCR product | {pcr_product.name} | "
                    f"record {pcr_product.record_id or 'n/a'} | "
                    f"{pcr_product.length} bp | single record"
                ),
            )
        )
    asset_choices.extend(get_template_asset_choices())
    return asset_choices


class BaseCloningConstructForm(forms.Form):
    def clean(self):
        cleaned_data = super().clean()
        vector_asset = cleaned_data.get("vector_asset")
        insert_asset = cleaned_data.get("insert_asset")
        if vector_asset and insert_asset and vector_asset == insert_asset:
            raise forms.ValidationError("Vector and insert must be different assets.")
        return cleaned_data


class CloningConstructAssetForm(BaseCloningConstructForm):
    name = forms.CharField(max_length=255)
    description = forms.CharField(required=False, widget=forms.Textarea)
    assembly_strategy = forms.ChoiceField(choices=CloningConstruct.STRATEGY_CHOICES)
    vector_asset = forms.ChoiceField(choices=())
    insert_asset = forms.ChoiceField(choices=())

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)
        asset_choices = _build_cloning_asset_choices(user)
        self.fields["vector_asset"].choices = asset_choices
        self.fields["insert_asset"].choices = list(asset_choices)
        self.fields["vector_asset"].help_text = "Choose an uploaded sequence-file record, template, or saved PCR product to use as the vector."
        self.fields["insert_asset"].help_text = "Choose an uploaded sequence-file record, template, or saved PCR product to insert into the vector."
        apply_tailwind_classes(self.fields)


def _configure_fragment_choice_field(field, *, choices, enzyme_name: str, asset_role: str):
    field.widget = forms.Select()
    if len(choices) > 1:
        field.choices = [("", f"Select a {asset_role} fragment")] + choices
    else:
        field.choices = choices or [("", f"No {enzyme_name} fragments available")]
    field.required = False
    field.help_text = f"Choose which {enzyme_name} fragment to keep from the {asset_role} asset."


def _ensure_choice(choices, value, label=None):
    normalized_value = str(value or "").strip()
    if not normalized_value:
        return choices
    if any(str(choice_value) == normalized_value for choice_value, _ in choices):
        return choices
    return list(choices) + [(normalized_value, label or normalized_value)]


def _format_circular_choice(value):
    if isinstance(value, bool):
        return "1" if value else "0"
    normalized = str(value or "").strip().lower()
    return "1" if normalized in {"1", "true", "circular", "yes", "on"} else "0"


class CloningConstructAssemblyForm(BaseCloningConstructForm):
    RESULT_TOPOLOGY_CHOICES = (
        ("1", "Circular"),
        ("0", "Linear"),
    )

    name = forms.CharField(max_length=255, widget=forms.HiddenInput)
    description = forms.CharField(required=False, widget=forms.HiddenInput)
    assembly_strategy = forms.ChoiceField(choices=CloningConstruct.STRATEGY_CHOICES, widget=forms.HiddenInput)
    vector_asset = forms.ChoiceField(choices=(), widget=forms.HiddenInput)
    insert_asset = forms.ChoiceField(choices=(), widget=forms.HiddenInput)
    is_circular = forms.ChoiceField(
        choices=RESULT_TOPOLOGY_CHOICES,
        required=False,
        label="Result topology",
        help_text="Choose circular when the assembled sequence should be treated as a plasmid.",
    )
    left_enzyme = forms.ChoiceField(choices=(), required=False)
    right_enzyme = forms.ChoiceField(choices=(), required=False)
    selected_enzymes = forms.MultipleChoiceField(
        choices=(),
        required=False,
        widget=forms.SelectMultiple,
        label="Restriction enzymes",
        help_text="Search and select enzymes to show cut sites and digest fragments on the vector and insert maps.",
    )
    vector_fragment_index = forms.ChoiceField(choices=(), required=False)
    insert_fragment_index = forms.ChoiceField(choices=(), required=False)

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        self.user = user
        super().__init__(*args, **kwargs)
        asset_choices = _build_cloning_asset_choices(user)
        if self.is_bound:
            asset_choices = _ensure_choice(
                asset_choices,
                self.data.get("vector_asset"),
                "Submitted vector asset",
            )
            asset_choices = _ensure_choice(
                asset_choices,
                self.data.get("insert_asset"),
                "Submitted insert asset",
            )
        self.fields["vector_asset"].choices = asset_choices
        self.fields["insert_asset"].choices = list(asset_choices)
        selected_vector_asset = (
            self.data.get("vector_asset")
            if self.is_bound
            else self.initial.get("vector_asset")
        )
        if not self.is_bound:
            if self.initial.get("is_circular") not in (None, ""):
                self.initial["is_circular"] = _format_circular_choice(self.initial.get("is_circular"))
            else:
                self.initial["is_circular"] = _format_circular_choice(
                    self._default_result_is_circular(selected_vector_asset)
                )
        enzyme_choices = (
            get_detected_enzyme_choices(
                user=user,
                selected_asset_choice=selected_vector_asset,
            )
            if user is not None
            else []
        )
        self.fields["left_enzyme"].choices = [("", "Select left enzyme")] + list(enzyme_choices)
        self.fields["right_enzyme"].choices = [("", "Select right enzyme")] + list(enzyme_choices)
        self.fields["selected_enzymes"].choices = list(enzyme_choices)
        if self.is_bound:
            self.fields["left_enzyme"].choices = _ensure_choice(
                self.fields["left_enzyme"].choices,
                self.data.get("left_enzyme"),
            )
            self.fields["right_enzyme"].choices = _ensure_choice(
                self.fields["right_enzyme"].choices,
                self.data.get("right_enzyme"),
            )
            selected_enzyme_values = self.data.getlist("selected_enzymes")
            for selected_enzyme in selected_enzyme_values:
                self.fields["selected_enzymes"].choices = _ensure_choice(
                    self.fields["selected_enzymes"].choices,
                    selected_enzyme,
                )
        elif not self.initial.get("selected_enzymes"):
            self.initial["selected_enzymes"] = [
                enzyme_name
                for enzyme_name in (
                    self.initial.get("left_enzyme"),
                    self.initial.get("right_enzyme"),
                )
                if enzyme_name
            ]
        self.fields["left_enzyme"].widget = forms.HiddenInput()
        self.fields["right_enzyme"].widget = forms.HiddenInput()

        selected_insert_asset = (
            self.data.get("insert_asset")
            if self.is_bound
            else self.initial.get("insert_asset")
        )
        selected_left_enzyme = (
            self.data.get("left_enzyme")
            if self.is_bound
            else self.initial.get("left_enzyme")
        )
        selected_right_enzyme = (
            self.data.get("right_enzyme")
            if self.is_bound
            else self.initial.get("right_enzyme")
        )
        selected_left_enzyme = str(selected_left_enzyme or "").strip()
        selected_right_enzyme = str(selected_right_enzyme or "").strip()
        same_enzyme_name = selected_left_enzyme if selected_left_enzyme and selected_left_enzyme == selected_right_enzyme else None
        if same_enzyme_name and user is not None and selected_vector_asset and selected_insert_asset:
            try:
                vector_asset = resolve_asset_choice(user=user, choice=selected_vector_asset)
                insert_asset = resolve_asset_choice(user=user, choice=selected_insert_asset)
                vector_choices = build_digest_fragment_choices(
                    sequence=vector_asset.sequence,
                    enzyme_name=same_enzyme_name,
                    is_circular=vector_asset.is_circular,
                )
                insert_choices = build_digest_fragment_choices(
                    sequence=insert_asset.sequence,
                    enzyme_name=same_enzyme_name,
                    is_circular=insert_asset.is_circular,
                )
            except ValueError:
                vector_choices = []
                insert_choices = []

            _configure_fragment_choice_field(
                self.fields["vector_fragment_index"],
                choices=vector_choices,
                enzyme_name=same_enzyme_name,
                asset_role="vector",
            )
            _configure_fragment_choice_field(
                self.fields["insert_fragment_index"],
                choices=insert_choices,
                enzyme_name=same_enzyme_name,
                asset_role="insert",
            )
        apply_tailwind_classes(self.fields)
        self.fields["selected_enzymes"].widget.attrs.update(
            {
                "class": "hidden",
                "data-cloning-selected-enzymes-control": "1",
                "aria-hidden": "true",
            }
        )
        self.fields["vector_fragment_index"].widget = forms.HiddenInput()
        self.fields["insert_fragment_index"].widget = forms.HiddenInput()

    def _default_result_is_circular(self, selected_vector_asset):
        if not selected_vector_asset:
            return False
        try:
            vector_asset = resolve_asset_choice(user=self.user, choice=selected_vector_asset)
        except ValueError:
            return False
        return vector_asset.is_circular

    def clean_is_circular(self):
        value = self.cleaned_data.get("is_circular")
        if value == "":
            return self._default_result_is_circular(self.cleaned_data.get("vector_asset"))
        return _format_circular_choice(value) == "1"

    def clean_selected_enzymes(self):
        selected_enzymes = [
            str(enzyme_name or "").strip()
            for enzyme_name in self.cleaned_data.get("selected_enzymes", [])
            if str(enzyme_name or "").strip()
        ]
        if selected_enzymes:
            return selected_enzymes
        return [
            enzyme_name
            for enzyme_name in (
                self.cleaned_data.get("left_enzyme"),
                self.cleaned_data.get("right_enzyme"),
            )
            if enzyme_name
        ]


class CloningConstructSequenceFileForm(forms.Form):
    name = forms.CharField(max_length=255)
    description = forms.CharField(required=False, widget=forms.Textarea)
    file_type = forms.ChoiceField(
        choices=[
            (SequenceFile.FILE_GENBANK, "GenBank"),
            (SequenceFile.FILE_FASTA, "FASTA"),
        ]
    )

    def __init__(self, *args, **kwargs):
        construct = kwargs.pop("construct", None)
        super().__init__(*args, **kwargs)
        if construct is not None and not self.is_bound:
            self.fields["name"].initial = construct.name
            self.fields["description"].initial = construct.description
        self.fields["file_type"].initial = SequenceFile.FILE_GENBANK
        self.fields["name"].help_text = "Name for the saved sequence file."
        self.fields["description"].help_text = "Optional description copied into the new sequence file."
        self.fields["file_type"].help_text = "GenBank keeps annotations; FASTA stores only the sequence."
        apply_tailwind_classes(self.fields)

class PrimerForm(forms.ModelForm):
    class Meta:
        model = Primer
        fields = ["primer_name", "overhang_sequence", "sequence"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        apply_tailwind_classes(self.fields)

    def clean_sequence(self):
        return clean_sequence_value(
            self.cleaned_data.get("sequence", ""),
            allow_n=False,
            max_length=60,
        )

    def clean_overhang_sequence(self):
        return clean_optional_sequence_value(
            self.cleaned_data.get("overhang_sequence", ""),
            allow_n=False,
        )


class PrimerExcelUploadForm(forms.Form):
    excel_file = forms.FileField(label="Excel file (.xlsx)")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        apply_tailwind_classes(self.fields)


class PrimerExcelColumnMapForm(forms.Form):
    name_column = forms.ChoiceField(label="Name column")
    sequence_column = forms.ChoiceField(label="Sequence column")
    overhang_column = forms.ChoiceField(
        label="Overhang column (optional)", required=False
    )

    def __init__(self, *args, **kwargs):
        columns = kwargs.pop("columns", [])
        super().__init__(*args, **kwargs)
        choices = [(col, col) for col in columns]
        self.fields["name_column"].choices = choices
        self.fields["sequence_column"].choices = choices
        self.fields["overhang_column"].choices = [("", "None")] + choices
        apply_tailwind_classes(self.fields)

class PrimerPairForm(forms.ModelForm ):
    class Meta:
        model = PrimerPair
        fields = ["name", "forward_primer", "reverse_primer"]

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)
        self.fields["forward_primer"].queryset = accessible_primers(user)
        self.fields["reverse_primer"].queryset = accessible_primers(user)
        apply_tailwind_classes(self.fields)

    def clean(self):
        cleaned_data = super().clean()
        forward = cleaned_data.get("forward_primer")
        reverse = cleaned_data.get("reverse_primer")
        if forward and reverse and forward == reverse:
            raise forms.ValidationError(
                "Forward and reverse primers must be different."
            )
        return cleaned_data


class PrimerPairCreateCombinedForm(forms.Form):
    pair_name = forms.CharField(max_length=100)

    # forward primer fields
    forward_name = forms.CharField(max_length=100)
    forward_sequence = forms.CharField(widget=forms.Textarea)
    forward_overhang = forms.CharField(required=False)

    # reverse primer fields
    reverse_name = forms.CharField(max_length=100)
    reverse_sequence = forms.CharField(widget=forms.Textarea)
    reverse_overhang = forms.CharField(required=False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        apply_tailwind_classes(self.fields)

    def clean_forward_sequence(self):
        return clean_sequence_value(
            self.cleaned_data.get("forward_sequence", ""),
            allow_n=True,
            max_length=60,
        )

    def clean_reverse_sequence(self):
        return clean_sequence_value(
            self.cleaned_data.get("reverse_sequence", ""),
            allow_n=True,
            max_length=60,
        )

    def clean_forward_overhang(self):
        return clean_optional_sequence_value(
            self.cleaned_data.get("forward_overhang", ""),
            allow_n=False,
        )

    def clean_reverse_overhang(self):
        return clean_optional_sequence_value(
            self.cleaned_data.get("reverse_overhang", ""),
            allow_n=False,
        )


class PCRProductDiscoveryForm(forms.Form):
    primer_pair = forms.ModelChoiceField(queryset=PrimerPair.objects.none())
    sequence_file = forms.ModelChoiceField(queryset=SequenceFile.objects.none())
    max_mismatches = forms.IntegerField(min_value=0, max_value=5, initial=0)
    block_3prime_mismatch = forms.BooleanField(required=False, initial=True)

    def __init__(self, *args, **kwargs):
        user = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)
        if user is not None:
            self.fields["primer_pair"].queryset = accessible_primer_pairs(user).order_by("name")
            self.fields["sequence_file"].queryset = accessible_sequence_files(user).order_by("name")
        apply_tailwind_classes(self.fields)

class CustomUserCreationForm(UserCreationForm):
    class Meta:
        model = User
        fields = ["username", "email"]

class CustomAuthenticationForm(AuthenticationForm):
    def clean(self):
        username = self.cleaned_data.get("username")
        password = self.cleaned_data.get("password")

        if username:
            user_model = get_user_model()
            try:
                has_user = user_model.objects.filter(username__iexact=username).exists()
            except (OperationalError, ProgrammingError):
                raise forms.ValidationError(
                    "Authentication is unavailable because the user database is not ready."
                )

            if not has_user:
                raise forms.ValidationError(
                    "No account found with that username. Please register first."
                )

        if username and password:
            try:
                self.user_cache = authenticate(
                    self.request, username=username, password=password
                )
            except (OperationalError, ProgrammingError):
                raise forms.ValidationError(
                    "Authentication is unavailable because the user database is not ready."
                )
            else:
                if self.user_cache is None:
                    raise forms.ValidationError(
                        self.error_messages["invalid_login"],
                        code="invalid_login",
                    )
                self.confirm_login_allowed(self.user_cache)

        return self.cleaned_data


class Primer3GlobalArgsForm(forms.Form):
    # PRODUCT SIZE
    PRIMER_PRODUCT_SIZE_RANGE = forms.CharField(
        label="Product Size Range",
        initial="100-300",
        help_text="Format: min-max or multiple ranges, e.g. 100-300 400-600"
    )

    def clean_PRIMER_PRODUCT_SIZE_RANGE(self):
        value = self.cleaned_data.get("PRIMER_PRODUCT_SIZE_RANGE", "")
        normalized = " ".join(value.split())
        pattern = r"^\d+-\d+( \d+-\d+)*$"

        if not re.fullmatch(pattern, normalized):
            raise forms.ValidationError(
                "Enter product size ranges as 'min-max' with optional additional ranges separated by spaces."
            )

        for size_range in normalized.split():
            min_size, max_size = (int(part) for part in size_range.split("-", 1))
            if min_size > max_size:
                raise forms.ValidationError(
                    "Each product size range must have a minimum less than or equal to the maximum."
                )

        return normalized

    # PRIMER LENGTH
    PRIMER_OPT_SIZE = forms.IntegerField(initial=20)
    PRIMER_MIN_SIZE = forms.IntegerField(initial=18)
    PRIMER_MAX_SIZE = forms.IntegerField(initial=27)

    # MELTING TEMPERATURES
    PRIMER_OPT_TM = forms.FloatField(initial=60.0)
    PRIMER_MIN_TM = forms.FloatField(initial=57.0)
    PRIMER_MAX_TM = forms.FloatField(initial=63.0)

    # GC CONTENT
    PRIMER_MIN_GC = forms.FloatField(initial=40.0)
    PRIMER_MAX_GC = forms.FloatField(initial=60.0)
    PRIMER_OPT_GC_PERCENT = forms.FloatField(initial=50.0)

    # SELF-COMPLEMENTARITY
    PRIMER_MAX_SELF_ANY = forms.FloatField(initial=8.0)
    PRIMER_MAX_SELF_END = forms.FloatField(initial=3.0)

    # HAIRPINS / DIMERS (melting temperatures accepted)
    PRIMER_MAX_HAIRPIN_TH = forms.FloatField(initial=47.0, required=False)
    PRIMER_MAX_END_STABILITY = forms.FloatField(initial=9.0, required=False)

    # POLY-X (runs of same base)
    PRIMER_MAX_POLY_X = forms.IntegerField(initial=4)

    # MISPRIMING
    PRIMER_MAX_LIBRARY_MISPRIMING = forms.FloatField(initial=12.00, required=False)
    PRIMER_PAIR_MAX_LIBRARY_MISPRIMING = forms.FloatField(initial=20.00, required=False)

    # PRODUCT MELTING TEMPERATURE
    PRIMER_PRODUCT_OPT_TM = forms.FloatField(required=False)
    PRIMER_PRODUCT_MIN_TM = forms.FloatField(required=False)
    PRIMER_PRODUCT_MAX_TM = forms.FloatField(required=False)

    # CLAMPING
    PRIMER_GC_CLAMP = forms.IntegerField(initial=1)

    # SALT / CHEMISTRY CONDITIONS
    PRIMER_SALT_MONOVALENT = forms.FloatField(initial=50.0)
    PRIMER_SALT_DIVALENT = forms.FloatField(initial=1.5)
    PRIMER_DNTP_CONC = forms.FloatField(initial=0.2)
    PRIMER_DNA_CONC = forms.FloatField(initial=50.0)

    # OTHER SETTINGS
    PRIMER_EXPLAIN_FLAG = forms.IntegerField(initial=1)
    PRIMER_NUM_RETURN = forms.IntegerField(initial=5)
    PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = forms.BooleanField(
        required=False, initial=True
    )
    PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = forms.BooleanField(
        required=False, initial=True
    )

    PRIMER_SIDE_CHOICES = (
        ("LEFT", "Forward (left) primers"),
        ("RIGHT", "Reverse (right) primers"),
    )

    PRIMER_SIDES = forms.MultipleChoiceField(
        choices=PRIMER_SIDE_CHOICES,
        widget=forms.CheckboxSelectMultiple,
        required=False,
        initial=["LEFT", "RIGHT"],
        label="Generate primers",
    )
