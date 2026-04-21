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
    FILE_SNAPGENE = "snapgene"

    FILE_TYPE_CHOICES = [
        (FILE_FASTA, "FASTA"),
        (FILE_GENBANK, "GenBank"),
        (FILE_SNAPGENE, "SnapGene (.dna)"),
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
    users = models.ManyToManyField(
        User,
        related_name="sequencefile_access",
        blank=True,
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
    pcr_products = models.ManyToManyField(
        "PCRProduct",
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
    def binding_sequence(self):
        return (self.sequence or "").strip().upper()

    @property
    def full_sequence(self):
        return f"{self.overhang_sequence}{self.binding_sequence}"

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


class PCRProduct(AccessControllModel):
    name = models.CharField(max_length=255)
    sequence_file = models.ForeignKey(
        SequenceFile,
        on_delete=models.CASCADE,
        related_name="pcr_products",
    )
    record_id = models.CharField(max_length=255)
    forward_primer = models.ForeignKey(
        Primer,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="pcr_products_as_forward",
    )
    reverse_primer = models.ForeignKey(
        Primer,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="pcr_products_as_reverse",
    )
    forward_feature = models.ForeignKey(
        SequenceFeature,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="pcr_products_as_forward_feature",
    )
    reverse_feature = models.ForeignKey(
        SequenceFeature,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="pcr_products_as_reverse_feature",
    )
    forward_primer_label = models.CharField(max_length=255, blank=True, default="")
    reverse_primer_label = models.CharField(max_length=255, blank=True, default="")
    start = models.IntegerField()
    end = models.IntegerField()
    sequence = models.TextField()
    length = models.IntegerField(default=0)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at", "name"]

    def save(self, *args, **kwargs):
        self.sequence = (self.sequence or "").strip().upper()
        self.length = len(self.sequence)
        super().save(*args, **kwargs)

    def __str__(self):
        return self.name


class AnalysisJob(models.Model):
    TYPE_PRIMER_BINDING = "primer_binding"
    TYPE_PCR_PRODUCT_DISCOVERY = "pcr_product_discovery"
    TYPE_CHOICES = [
        (TYPE_PRIMER_BINDING, "Primer binding"),
        (TYPE_PCR_PRODUCT_DISCOVERY, "PCR product discovery"),
    ]

    STATUS_PENDING = "pending"
    STATUS_RUNNING = "running"
    STATUS_SUCCESS = "success"
    STATUS_FAILURE = "failure"
    STATUS_CHOICES = [
        (STATUS_PENDING, "Pending"),
        (STATUS_RUNNING, "Running"),
        (STATUS_SUCCESS, "Success"),
        (STATUS_FAILURE, "Failure"),
    ]

    owner = models.ForeignKey(
        User,
        on_delete=models.CASCADE,
        related_name="analysis_jobs",
    )
    job_type = models.CharField(max_length=50, choices=TYPE_CHOICES)
    celery_task_id = models.CharField(max_length=255, blank=True, default="")
    status = models.CharField(
        max_length=20,
        choices=STATUS_CHOICES,
        default=STATUS_PENDING,
    )
    error_message = models.TextField(blank=True, default="")
    result_payload = models.JSONField(null=True, blank=True)
    target_object_id = models.PositiveBigIntegerField(null=True, blank=True)
    sequence_file = models.ForeignKey(
        SequenceFile,
        on_delete=models.CASCADE,
        related_name="analysis_jobs",
        null=True,
        blank=True,
    )
    primer = models.ForeignKey(
        Primer,
        on_delete=models.CASCADE,
        related_name="analysis_jobs",
        null=True,
        blank=True,
    )
    primer_pair = models.ForeignKey(
        PrimerPair,
        on_delete=models.CASCADE,
        related_name="analysis_jobs",
        null=True,
        blank=True,
    )
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ["-created_at"]

    def __str__(self):
        return f"{self.job_type}:{self.owner_id}:{self.status}"


class CloningConstruct(AccessControllModel):
    STRATEGY_RESTRICTION_LIGATION = "restriction_ligation"
    STRATEGY_CHOICES = [
        (STRATEGY_RESTRICTION_LIGATION, "Restriction ligation"),
    ]
    SOURCE_SEQUENCE_FILE = "sequence_file"
    SOURCE_PCR_PRODUCT = "pcr_product"
    SOURCE_CHOICES = [
        (SOURCE_SEQUENCE_FILE, "Sequence file"),
        (SOURCE_PCR_PRODUCT, "PCR product"),
    ]

    name = models.CharField(max_length=255)
    description = models.TextField(blank=True, default="")
    vector_source_type = models.CharField(max_length=20, choices=SOURCE_CHOICES)
    vector_sequence_file = models.ForeignKey(
        SequenceFile,
        on_delete=models.CASCADE,
        related_name="cloning_constructs_as_vector",
        null=True,
        blank=True,
    )
    vector_pcr_product = models.ForeignKey(
        "PCRProduct",
        on_delete=models.CASCADE,
        related_name="cloning_constructs_as_vector",
        null=True,
        blank=True,
    )
    vector_record_id = models.CharField(max_length=255, blank=True, default="")
    insert_source_type = models.CharField(max_length=20, choices=SOURCE_CHOICES)
    insert_sequence_file = models.ForeignKey(
        SequenceFile,
        on_delete=models.CASCADE,
        related_name="cloning_constructs_as_insert",
        null=True,
        blank=True,
    )
    insert_pcr_product = models.ForeignKey(
        "PCRProduct",
        on_delete=models.CASCADE,
        related_name="cloning_constructs_as_insert",
        null=True,
        blank=True,
    )
    insert_record_id = models.CharField(max_length=255, blank=True, default="")
    assembly_strategy = models.CharField(
        max_length=50,
        choices=STRATEGY_CHOICES,
        default=STRATEGY_RESTRICTION_LIGATION,
    )
    left_enzyme = models.CharField(max_length=100)
    right_enzyme = models.CharField(max_length=100)
    assembled_sequence = models.TextField(blank=True, default="")
    assembled_length = models.IntegerField(default=0)
    is_valid = models.BooleanField(default=False)
    validation_messages = models.JSONField(default=list, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)

    class Meta:
        ordering = ["-created_at", "name"]

    def save(self, *args, **kwargs):
        self.assembled_sequence = (self.assembled_sequence or "").strip().upper()
        self.assembled_length = len(self.assembled_sequence)
        super().save(*args, **kwargs)

    def __str__(self):
        return self.name

    @property
    def vector_name(self):
        if self.vector_sequence_file_id:
            return self.vector_sequence_file.name
        if self.vector_pcr_product_id:
            return self.vector_pcr_product.name
        return "Unknown vector"

    @property
    def insert_name(self):
        if self.insert_sequence_file_id:
            return self.insert_sequence_file.name
        if self.insert_pcr_product_id:
            return self.insert_pcr_product.name
        return "Unknown insert"
