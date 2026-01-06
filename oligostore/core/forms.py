from django.contrib.auth import authenticate, get_user_model
from django.contrib.auth.forms import AuthenticationForm, UserCreationForm
from django.db import OperationalError, ProgrammingError
from django.contrib.auth.models import User
from django import forms
from .models import Primer, PrimerPair, Project
import re

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
        self.fields["forward_primer"].queryset = Primer.objects.filter(users=user)
        self.fields["reverse_primer"].queryset = Primer.objects.filter(users=user)
        apply_tailwind_classes(self.fields)

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