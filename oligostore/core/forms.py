from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
from django import forms
from .models import Primer, PrimerPair, Project
import re

class ProjectForm(forms.ModelForm):
    class Meta:
        model = Project
        fields = ["name", "description"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["name"].widget.attrs.update({"class": "input input-bordered w-full"})
        self.fields["description"].widget.attrs.update({"class": "textarea textarea-bordered w-full"})

class PrimerForm(forms.ModelForm):
    class Meta:
        model = Primer
        fields = ["primer_name", "sequence"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields["primer_name"].widget.attrs.update({"class": "input input-bordered w-full"})
        self.fields["sequence"].widget.attrs.update({"class": "textarea textarea-bordered w-full"})

    def clean_sequence(self):
        seq = self.cleaned_data.get("sequence", "").upper()

        # Allow only A, C, G, T
        if not re.fullmatch(r"[ACGT]+", seq):
            raise forms.ValidationError(
                "Sequence may only contain the characters A, C, G, T (no spaces or numbers)."
            )

        return seq

class PrimerPairForm(forms.ModelForm):
    class Meta:
        model = PrimerPair
        fields = ["name", "forward_primer", "reverse_primer"]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Show ALL primers, since access should NOT be limited
        self.fields["forward_primer"].queryset = Primer.objects.all()
        self.fields["reverse_primer"].queryset = Primer.objects.all()

        # Add DaisyUI/Tailwind classes for styling
        self.fields["name"].widget.attrs.update({"class": "input input-bordered w-full"})
        self.fields["forward_primer"].widget.attrs.update({"class": "select select-bordered w-full"})
        self.fields["reverse_primer"].widget.attrs.update({"class": "select select-bordered w-full"})

class PrimerPairCreateCombinedForm(forms.Form):
    pair_name = forms.CharField(max_length=100)

    # forward primer fields
    forward_name = forms.CharField(max_length=100)
    forward_sequence = forms.CharField(widget=forms.Textarea)

    # reverse primer fields
    reverse_name = forms.CharField(max_length=100)
    reverse_sequence = forms.CharField(widget=forms.Textarea)

    # styling
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for k, field in self.fields.items():
            if k.endswith("sequence"):
                field.widget.attrs.update({"class": "textarea textarea-bordered w-full"})
            else:
                field.widget.attrs.update({"class": "input input-bordered w-full"})


class CustomUserCreationForm(UserCreationForm):
    class Meta:
        model = User
        fields = ["username", "email"]


class Primer3GlobalArgsForm(forms.Form):
    # PRODUCT SIZE
    PRIMER_PRODUCT_SIZE_RANGE = forms.CharField(
        label="Product Size Range",
        initial="100-300",
        help_text="Format: min-max or multiple ranges, e.g. 100-300 400-600"
    )

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