from django.db import models
from django.contrib.auth.models import User
from functools import lru_cache
from Bio.Seq import Seq
from Bio.Restriction import CommOnly


@lru_cache(maxsize=1024)
def _find_overhang_restriction_sites(sequence: str):
    sequence = (sequence or "").upper().strip()
    if not sequence:
        return tuple()

    hits = []
    try:
        restriction_results = CommOnly.search(Seq(sequence), linear=True)
        for enzyme, cut_positions in restriction_results.items():
            site = str(getattr(enzyme, "site", "") or "")
            site_length = len(site)
            if site_length <= 0:
                continue

            cut_offset = int(getattr(enzyme, "fst5", 0))
            for cut_position in cut_positions:
                start = int(cut_position) - cut_offset
                end = start + site_length - 1
                if start < 1 or end > len(sequence):
                    continue
                hits.append(
                    {
                        "enzyme": str(enzyme),
                        "site": site,
                        "cut_offset": cut_offset,
                        "start": start,
                        "end": end,
                    }
                )
    except Exception:
        return tuple()

    return tuple(sorted(hits, key=lambda hit: (hit["start"], hit["enzyme"])))

class AccessControllModel(models.Model):
    creator = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="%(class)s_created",
    )

    users = models.ManyToManyField(
        User,
        related_name="%(class)s_access",
        blank=True,
    )

    class Meta:
        abstract = True

class SequenceFile(models.Model):
    FILE_FASTA = "fasta"
    FILE_GENBANK = "genbank"

    FILE_TYPE_CHOICES = [
        (FILE_FASTA, "FASTA"),
        (FILE_GENBANK, "GenBank"),
    ]

    name = models.CharField(max_length=255)
    file = models.FileField(upload_to="sequence_files/")
    file_type = models.CharField(
        max_length=20,
        choices=FILE_TYPE_CHOICES,
    )

    uploaded_by = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="sequence_files",
    )

    uploaded_at = models.DateTimeField(auto_now_add=True)

    description = models.TextField(blank=True)

    class Meta:
        ordering = ["-uploaded_at"]

    def __str__(self):
        return self.name


class SequenceFeature(models.Model):
    TYPE_PRIMER_BIND = "primer_bind"
    TYPE_CUSTOM = "custom"

    TYPE_CHOICES = [
        (TYPE_PRIMER_BIND, "Primer binding site"),
        (TYPE_CUSTOM, "Custom"),
    ]

    STRAND_FORWARD = 1
    STRAND_REVERSE = -1
    STRAND_CHOICES = [
        (STRAND_FORWARD, "Forward (+)"),
        (STRAND_REVERSE, "Reverse (-)"),
    ]

    sequence_file = models.ForeignKey(
        SequenceFile,
        on_delete=models.CASCADE,
        related_name="user_features",
    )
    primer = models.ForeignKey(
        "Primer",
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="sequence_features",
    )
    record_id = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.IntegerField(choices=STRAND_CHOICES, default=STRAND_FORWARD)
    feature_type = models.CharField(
        max_length=50,
        choices=TYPE_CHOICES,
        default=TYPE_PRIMER_BIND,
    )
    label = models.CharField(max_length=255)
    created_by = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="created_sequence_features",
    )
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["record_id", "start", "end", "label"]

    def __str__(self):
        return f"{self.sequence_file_id}:{self.record_id}:{self.label}"


class PrimerBindingResult(models.Model):
    primer = models.ForeignKey("Primer", on_delete=models.CASCADE)
    sequence_file = models.ForeignKey(SequenceFile, on_delete=models.CASCADE)

    record_id = models.CharField(max_length=255)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.CharField(max_length=1)
    mismatches = models.IntegerField()

    created_at = models.DateTimeField(auto_now_add=True)

class Project(AccessControllModel):
    name = models.CharField(max_length=100)
    description = models.TextField(blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    primerpairs = models.ManyToManyField(
        "PrimerPair",
        related_name="projects",
        blank=True
    )
    sequence_files = models.ManyToManyField(
        SequenceFile,
        related_name="projects",
        blank=True,
    )

    def __str__(self):
        return self.name

class Primer(AccessControllModel):
    primer_name = models.CharField(max_length=100)
    sequence = models.TextField()
    overhang_sequence = models.TextField(null=True, blank=True, default="")
    length = models.IntegerField(null=True,blank=True)

    # --- Analysis fields ---
    gc_content = models.FloatField(null=True, blank=True)
    tm = models.FloatField(null=True, blank=True)
    hairpin_dg = models.FloatField(null=True, blank=True)
    self_dimer_dg = models.FloatField(null=True, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.primer_name

    @property
    def full_sequence(self):
        return f"{self.overhang_sequence}{self.sequence}"

    @property
    def overhang_restriction_sites(self):
        return list(_find_overhang_restriction_sites(self.overhang_sequence or ""))

    @property
    def restriction_site_summary(self):
        hits = self.overhang_restriction_sites
        if not hits:
            return ""
        return ", ".join(
            f"{hit['enzyme']} ({hit['site']})" for hit in hits
        )

    @classmethod
    def create_with_analysis(
        cls, *, primer_name, sequence, user, overhang_sequence=""
    ):
        from .services.primer_analysis import analyze_primer

        analysis = analyze_primer(sequence)
        primer = cls(
            primer_name=primer_name,
            sequence=sequence,
            overhang_sequence=(overhang_sequence or "").strip().upper(),
            length=len(sequence),
            gc_content=analysis["gc_content"],
            tm=analysis["tm"],
            hairpin_dg=analysis["hairpin_dg"],
            self_dimer_dg=analysis["self_dimer_dg"],
            creator=user,
        )
        primer.save()
        primer.users.add(user)
        return primer

class PrimerPair(AccessControllModel):
    name = models.CharField(max_length=100)
    description = models.TextField(null=True,
                                   blank=True)

    forward_primer = models.ForeignKey(
        Primer,
        on_delete=models.CASCADE,
        related_name="as_forward_primer",
    )
    reverse_primer = models.ForeignKey(
        Primer,
        on_delete=models.CASCADE,
        related_name="as_reverse_primer",
    )

    def __str__(self):
        return self.name
