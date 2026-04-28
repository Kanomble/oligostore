from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings

from core.access import accessible_pcr_products, accessible_sequence_files
from core.models import (
    AnalysisJob,
    CloningConstruct,
    PCRProduct,
    Primer,
    PrimerBindingResult,
    PrimerPair,
    Project,
    SequenceFeature,
    SequenceFile,
)

import shutil
import tempfile


class SequenceFileModelTests(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls._temp_media_root = tempfile.mkdtemp(prefix="oligostore-tests-")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._temp_media_root)
        super().tearDownClass()

    def setUp(self):
        super().setUp()
        self._media_override = override_settings(MEDIA_ROOT=self._temp_media_root)
        self._media_override.enable()
        self.user = User.objects.create_user(
            username="uploader",
            email="uploader@example.com",
            password="testpass123",
        )

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_sequence_file_str_and_ordering(self):
        first_file = SimpleUploadedFile("first.fasta", b">seq1\nATCG")
        second_file = SimpleUploadedFile("second.gb", b"LOCUS")

        first = SequenceFile.objects.create(
            name="First",
            file=first_file,
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
            description="First upload",
        )
        second = SequenceFile.objects.create(
            name="Second",
            file=second_file,
            file_type=SequenceFile.FILE_GENBANK,
            uploaded_by=self.user,
            description="Second upload",
        )

        self.assertEqual(str(first), "First")
        self.assertEqual(str(second), "Second")

        ordered = list(SequenceFile.objects.all())
        self.assertEqual(ordered[0], second)
        self.assertEqual(ordered[1], first)

    def test_sequence_feature_str_and_defaults(self):
        sequence = SequenceFile.objects.create(
            name="Feature Seq",
            file=SimpleUploadedFile("feature.fasta", b">r1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        primer = Primer.objects.create(
            primer_name="FeaturePrimer",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        feature = SequenceFeature.objects.create(
            sequence_file=sequence,
            primer=primer,
            record_id="r1",
            start=2,
            end=9,
            strand=SequenceFeature.STRAND_FORWARD,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="FeaturePrimer",
            created_by=self.user,
        )
        self.assertEqual(str(feature), f"{sequence.id}:r1:FeaturePrimer")
        self.assertEqual(feature.feature_type, SequenceFeature.TYPE_PRIMER_BIND)


class PrimerAndProjectModelTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        cls.creator = User.objects.create_user(
            username="creator",
            email="creator@example.com",
            password="testpass123",
        )
        cls.collaborator = User.objects.create_user(
            username="collab",
            email="collab@example.com",
            password="testpass123",
        )

        cls.forward_primer = Primer.objects.create(
            primer_name="Forward",
            sequence="ATCGATCG",
            length=8,
            creator=cls.creator,
        )
        cls.reverse_primer = Primer.objects.create(
            primer_name="Reverse",
            sequence="CGATCGAT",
            length=8,
            creator=cls.creator,
        )

        cls.primer_pair = PrimerPair.objects.create(
            name="Pair A",
            description="Example pair",
            forward_primer=cls.forward_primer,
            reverse_primer=cls.reverse_primer,
            creator=cls.creator,
        )

    def test_primer_str(self):
        self.assertEqual(str(self.forward_primer), "Forward")

    def test_create_with_analysis_sets_analysis_fields(self):
        primer = Primer.create_with_analysis(
            primer_name="Analyzed Primer",
            sequence="ATGCATGC",
            overhang_sequence="aagcttggatcc",
            user=self.creator,
        )

        self.assertEqual(primer.creator, self.creator)
        self.assertIn(self.creator, primer.users.all())
        self.assertEqual(primer.length, 8)
        self.assertIsNotNone(primer.gc_content)
        self.assertIsNotNone(primer.tm)
        self.assertIsNotNone(primer.hairpin_dg)
        self.assertIsNotNone(primer.self_dimer_dg)
        self.assertEqual(primer.overhang_sequence, "AAGCTTGGATCC")
        self.assertIn("HindIII", primer.restriction_site_summary)
        self.assertIn("BamHI", primer.restriction_site_summary)

    def test_binding_sequence_excludes_overhang(self):
        primer = Primer.create_with_analysis(
            primer_name="Binding Primer",
            sequence=" atgcatgc ",
            overhang_sequence="GGATCC",
            user=self.creator,
        )

        self.assertEqual(primer.binding_sequence, "ATGCATGC")
        self.assertEqual(primer.full_sequence, "GGATCCATGCATGC")

    def test_overhang_detects_restriction_sites_beyond_legacy_set(self):
        primer = Primer.create_with_analysis(
            primer_name="Expanded Restriction Primer",
            sequence="ATGCATGC",
            overhang_sequence="GGGCCC",  # ApaI site, not part of legacy hardcoded set
            user=self.creator,
        )

        self.assertIn("ApaI", primer.restriction_site_summary)

    def test_project_access_and_pairs(self):
        project = Project.objects.create(
            name="Project Alpha",
            description="Demo",
            creator=self.creator,
        )
        project.users.add(self.collaborator)
        project.primerpairs.add(self.primer_pair)

        self.assertEqual(str(project), "Project Alpha")
        self.assertIn(self.collaborator, project.users.all())
        self.assertIn(self.primer_pair, project.primerpairs.all())

    def test_pcr_product_save_normalizes_sequence_and_length(self):
        sequence_file = SequenceFile.objects.create(
            name="PCR Source",
            file=SimpleUploadedFile("pcr_source.fasta", b">r1\nATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.creator,
        )
        product = PCRProduct.objects.create(
            name="Amplicon 1",
            sequence_file=sequence_file,
            record_id="r1",
            forward_primer=self.forward_primer,
            reverse_primer=self.reverse_primer,
            forward_primer_label="Forward",
            reverse_primer_label="Reverse",
            start=1,
            end=4,
            sequence=" atcg ",
            creator=self.creator,
        )
        product.users.add(self.creator)

        self.assertEqual(str(product), "Amplicon 1")
        self.assertEqual(product.sequence, "ATCG")
        self.assertEqual(product.length, 4)

    def test_cloning_construct_normalizes_assembled_sequence_and_asset_names(self):
        vector = SequenceFile.objects.create(
            name="Vector A",
            description="Backbone",
            file=SimpleUploadedFile("vector_a.fasta", b">vecA\nGAATTCAAAAGGATCC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.creator,
        )
        vector.users.add(self.creator)
        insert_source = SequenceFile.objects.create(
            name="Insert Source",
            description="Insert",
            file=SimpleUploadedFile("insert_a.fasta", b">insA\nATGC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.creator,
        )
        insert_source.users.add(self.creator)
        construct = CloningConstruct.objects.create(
            name="Construct A",
            vector_source_type=CloningConstruct.SOURCE_SEQUENCE_FILE,
            vector_sequence_file=vector,
            vector_record_id="vecA",
            insert_source_type=CloningConstruct.SOURCE_SEQUENCE_FILE,
            insert_sequence_file=insert_source,
            insert_record_id="insA",
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
            assembled_sequence=" gaattcatgcggatcc ",
            creator=self.creator,
        )
        construct.users.add(self.creator)

        self.assertEqual(construct.assembled_sequence, "GAATTCATGCGGATCC")
        self.assertEqual(construct.assembled_length, 16)
        self.assertEqual(construct.vector_name, "Vector A")
        self.assertEqual(construct.insert_name, "Insert Source")
        self.assertEqual(construct.vector_asset_label, "Vector A / vecA")
        self.assertEqual(construct.insert_asset_label, "Insert Source / insA")


class PrimerBindingResultModelTests(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls._temp_media_root = tempfile.mkdtemp(prefix="oligostore-tests-")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._temp_media_root)
        super().tearDownClass()

    def setUp(self):
        super().setUp()
        self._media_override = override_settings(MEDIA_ROOT=self._temp_media_root)
        self._media_override.enable()
        self.user = User.objects.create_user(
            username="binding_user",
            email="binding@example.com",
            password="testpass123",
        )
        self.primer = Primer.objects.create(
            primer_name="Primer X",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        sequence_file = SimpleUploadedFile("binding.fasta", b">seq\nATCG")
        self.sequence = SequenceFile.objects.create(
            name="Binding",
            file=sequence_file,
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
            description="Binding test",
        )

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_primer_binding_result_fields(self):
        result = PrimerBindingResult.objects.create(
            primer=self.primer,
            sequence_file=self.sequence,
            record_id="record-1",
            start=10,
            end=20,
            strand="+",
            mismatches=1,
        )

        self.assertEqual(result.primer, self.primer)
        self.assertEqual(result.sequence_file, self.sequence)
        self.assertEqual(result.record_id, "record-1")
        self.assertEqual(result.start, 10)
        self.assertEqual(result.end, 20)
        self.assertEqual(result.strand, "+")
        self.assertEqual(result.mismatches, 1)


class AccessHelperTests(TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls._temp_media_root = tempfile.mkdtemp(prefix="oligostore-tests-")

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls._temp_media_root)
        super().tearDownClass()

    def setUp(self):
        super().setUp()
        self._media_override = override_settings(MEDIA_ROOT=self._temp_media_root)
        self._media_override.enable()
        self.user = User.objects.create_user(
            username="access_user",
            email="access@example.com",
            password="testpass123",
        )
        self.other_user = User.objects.create_user(
            username="access_other",
            email="access_other@example.com",
            password="testpass123",
        )

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_accessible_sequence_files_includes_shared_files(self):
        owned = SequenceFile.objects.create(
            name="Owned File",
            file=SimpleUploadedFile("owned.fasta", b">r1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        shared = SequenceFile.objects.create(
            name="Other File",
            file=SimpleUploadedFile("other.fasta", b">r1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.other_user,
        )
        shared.users.add(self.user)

        self.assertEqual(
            {item.name for item in accessible_sequence_files(self.user)},
            {owned.name, shared.name},
        )

    def test_accessible_pcr_products_includes_shared_records(self):
        sequence_file = SequenceFile.objects.create(
            name="PCR Access File",
            file=SimpleUploadedFile("pcr_access.fasta", b">r1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        owned = PCRProduct.objects.create(
            name="Owned Product",
            sequence_file=sequence_file,
            record_id="r1",
            start=1,
            end=4,
            sequence="ATCG",
            creator=self.user,
        )
        shared = PCRProduct.objects.create(
            name="Shared Product",
            sequence_file=sequence_file,
            record_id="r1",
            start=1,
            end=4,
            sequence="ATCG",
            creator=self.other_user,
        )
        shared.users.add(self.user)

        self.assertEqual(
            {product.name for product in accessible_pcr_products(self.user)},
            {owned.name, shared.name},
        )

    def test_analysis_job_str(self):
        sequence_file = SequenceFile.objects.create(
            name="Job File",
            file=SimpleUploadedFile("job.fasta", b">r1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        job = AnalysisJob.objects.create(
            owner=self.user,
            job_type=AnalysisJob.TYPE_PRIMER_BINDING,
            sequence_file=sequence_file,
            status=AnalysisJob.STATUS_PENDING,
        )

        self.assertIn("primer_binding", str(job))
