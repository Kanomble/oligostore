from django.db import models
from django.contrib.auth.models import User

RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "XhoI": "CTCGAG",
    "XbaI": "TCTAGA",
    "NheI": "GCTAGC",
    "SpeI": "ACTAGT",
    "NotI": "GCGGCCGC",
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
}

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
    overhang_sequence = models.TextField(blank=True, default="")
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
        sequence = (self.overhang_sequence or "").upper()
        if not sequence:
            return []

        hits = []
        for enzyme, site in RESTRICTION_ENZYMES.items():
            start = 0
            while True:
                idx = sequence.find(site, start)
                if idx == -1:
                    break
                hits.append(
                    {
                        "enzyme": enzyme,
                        "site": site,
                        "start": idx + 1,
                    }
                )
                start = idx + 1
        return hits

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