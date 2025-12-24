from django.db import models
from django.contrib.auth.models import User

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

class Project(models.Model):
    name = models.CharField(max_length=100)
    description = models.TextField(blank=True, null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    creator = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="created_projects"
    )

    users = models.ManyToManyField(
        User,
        related_name="projects",
        blank=True
    )

    primerpairs = models.ManyToManyField(
        "PrimerPair",
        related_name="projects",
        blank=True
    )

    def __str__(self):
        return self.name

class Primer(models.Model):
    primer_name = models.CharField(max_length=100)
    sequence = models.TextField()
    length = models.IntegerField(null=True,blank=True)

    creator = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="created_primers"
    )

    users = models.ManyToManyField(
        User,
        related_name="accessible_primers",
        blank=True
    )

    # --- Analysis fields ---
    gc_content = models.FloatField(null=True, blank=True)
    tm = models.FloatField(null=True, blank=True)
    hairpin_dg = models.FloatField(null=True, blank=True)
    self_dimer_dg = models.FloatField(null=True, blank=True)

    created_at = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.primer_name

class PrimerPair(models.Model):
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

    creator = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="created_primerpairs"
    )

    users = models.ManyToManyField(
        User,
        related_name="accessible_primerpairs",
        blank=True
    )

    def __str__(self):
        return self.name