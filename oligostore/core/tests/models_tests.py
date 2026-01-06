from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings

from core.models import (
    Primer,
    PrimerBindingResult,
    PrimerPair,
    Project,
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