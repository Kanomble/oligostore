from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse
from types import SimpleNamespace
from unittest.mock import patch
from core.models import (
    Primer,
    PrimerBindingResult,
    PrimerPair,
    Project,
    SequenceFile,
)

import shutil
import tempfile

class SequenceFileViewTests(TestCase):
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
            username="sequence_user",
            email="sequence@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_sequencefile_upload_get(self):
        response = self.client.get(reverse("sequencefile_upload"))
        self.assertEqual(response.status_code, 200)

    def test_sequencefile_upload_post_validation(self):
        response = self.client.post(reverse("sequencefile_upload"), {"name": "Missing"})
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "All required fields must be provided.")

    def test_sequencefile_upload_post_success(self):
        file_data = SimpleUploadedFile("upload.fasta", b">seq\nATCG")
        response = self.client.post(
            reverse("sequencefile_upload"),
            {
                "name": "Upload",
                "file": file_data,
                "file_type": SequenceFile.FILE_FASTA,
                "description": "Test upload",
            },
        )
        self.assertRedirects(response, reverse("sequencefile_list"))
        self.assertTrue(SequenceFile.objects.filter(name="Upload").exists())

    def test_sequencefile_list_filters(self):
        other_user = User.objects.create_user(
            username="other_user",
            email="other@example.com",
            password="testpass123",
        )
        SequenceFile.objects.create(
            name="Alpha",
            file=SimpleUploadedFile("alpha.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        SequenceFile.objects.create(
            name="Beta",
            file=SimpleUploadedFile("beta.gb", b"LOCUS"),
            file_type=SequenceFile.FILE_GENBANK,
            uploaded_by=self.user,
            description="Genbank file",
        )
        SequenceFile.objects.create(
            name="Gamma",
            file=SimpleUploadedFile("gamma.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=other_user,
        )

        response = self.client.get(reverse("sequencefile_list"), {"q": "beta"})
        self.assertEqual(response.status_code, 200)
        listed = list(response.context["sequence_files"])
        self.assertEqual(len(listed), 1)
        self.assertEqual(listed[0].name, "Beta")

        response = self.client.get(reverse("sequencefile_list"), {"type": "fasta"})
        listed = list(response.context["sequence_files"])
        self.assertTrue(all(item.file_type == SequenceFile.FILE_FASTA for item in listed))
        self.assertEqual({item.name for item in listed}, {"Alpha"})


class PrimerBindingViewTests(TestCase):
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
            username="binding_view_user",
            email="binding_view@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)
        self.primer = Primer.objects.create(
            primer_name="Primer B",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        self.primer.users.add(self.user)
        self.sequence_file = SequenceFile.objects.create(
            name="Sequence",
            file=SimpleUploadedFile("seq.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_primer_binding_analysis_get(self):
        response = self.client.get(reverse("primer_binding"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("primers", response.context)

    @patch("core.views.analyze_primer_binding")
    def test_primer_binding_analysis_post(self, analyze_primer_binding):
        analyze_primer_binding.return_value = [{"record_id": "rec1"}]
        response = self.client.post(
            reverse("primer_binding"),
            {
                "primer_id": self.primer.id,
                "sequence_file_id": self.sequence_file.id,
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["hits"], [{"record_id": "rec1"}])


class DownloadProductSequenceViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="download_user",
            email="download@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def test_download_product_sequence_success(self):
        response = self.client.post(
            reverse("download_product_sequence"),
            {"product_sequence": "ATCG", "pair_index": "3"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertIn(b">PCR_product_pair_3", response.content)
        self.assertEqual(
            response["Content-Disposition"],
            "attachment; filename=pcr_product_pair_3.fasta",
        )

    def test_download_product_sequence_missing_payload(self):
        response = self.client.post(reverse("download_product_sequence"), {})
        self.assertEqual(response.status_code, 400)


class AnalyzePrimerViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="analysis_user",
            email="analysis@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def test_analyze_primer_requires_post(self):
        response = self.client.get(reverse("analyze_primer"))
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json()["error"], "POST only")

    def test_analyze_primer_invalid_sequence(self):
        response = self.client.post(
            reverse("analyze_primer"),
            {"sequence": "ATCG123"},
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn("invalid characters", response.json()["error"])

    @patch("core.views.analyze_primer")
    def test_analyze_primer_valid_sequence(self, analyze_primer):
        analyze_primer.return_value = {"tm": 60.0}
        response = self.client.post(
            reverse("analyze_primer"),
            {"sequence": "ATCG"},
        )
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.json(), {"tm": 60.0})


class AnalyzePrimerPairViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="pair_user",
            email="pair@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def test_analyze_primerpair_requires_post(self):
        response = self.client.get(reverse("analyze_primerpair"))
        self.assertEqual(response.status_code, 400)
        self.assertEqual(response.json()["error"], "POST required")

    def test_analyze_primerpair_invalid_forward(self):
        response = self.client.post(
            reverse("analyze_primerpair"),
            {"forward_sequence": "ATCGN1", "reverse_sequence": "ATCG"},
        )
        self.assertEqual(response.status_code, 400)
        self.assertIn("Forward sequence contains invalid characters", response.json()["error"])

    @patch("core.views.analyze_cross_dimer")
    @patch("core.views.analyze_primer")
    def test_analyze_primerpair_success(self, analyze_primer, analyze_cross_dimer):
        analyze_primer.side_effect = [
            {"tm": 60.0, "gc_content": 0.5, "hairpin_dg": -1.0, "self_dimer_dg": -2.0},
            {"tm": 61.0, "gc_content": 0.4, "hairpin_dg": -1.5, "self_dimer_dg": -2.5},
        ]
        analyze_cross_dimer.return_value = SimpleNamespace(dg=-3.2)

        response = self.client.post(
            reverse("analyze_primerpair"),
            {"forward_sequence": "ATCG", "reverse_sequence": "CGAT"},
        )

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertEqual(payload["forward"]["tm"], 60.0)
        self.assertEqual(payload["reverse"]["tm"], 61.0)
        self.assertEqual(payload["hetero_dimer_dg"], -3.2)


class ProjectViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="project_user",
            email="project@example.com",
            password="testpass123",
        )
        self.other_user = User.objects.create_user(
            username="project_other",
            email="project_other@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)
        self.project = Project.objects.create(
            name="Project One",
            creator=self.user,
        )
        self.project.users.add(self.user)
        self.forward = Primer.objects.create(
            primer_name="Forward",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        self.reverse = Primer.objects.create(
            primer_name="Reverse",
            sequence="CGATCGAT",
            length=8,
            creator=self.user,
        )
        self.pair = PrimerPair.objects.create(
            name="Pair One",
            forward_primer=self.forward,
            reverse_primer=self.reverse,
            creator=self.user,
        )
        self.pair.users.add(self.user)

    def test_project_dashboard_forbidden_for_non_member(self):
        self.client.force_login(self.other_user)
        response = self.client.get(reverse("project_dashboard", args=[self.project.id]))
        self.assertEqual(response.status_code, 403)

    def test_project_add_and_remove_primerpair(self):
        response = self.client.get(
            reverse("project_add_primerpair", args=[self.project.id, self.pair.id])
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertIn(self.pair, self.project.primerpairs.all())

        response = self.client.get(
            reverse("project_remove_primerpair", args=[self.project.id, self.pair.id])
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertNotIn(self.pair, self.project.primerpairs.all())


class PrimerViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="primer_user",
            email="primer@example.com",
            password="testpass123",
        )
        self.other_user = User.objects.create_user(
            username="primer_other",
            email="primer_other@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    @patch("core.views.analyze_primer")
    def test_primer_create_and_list(self, analyze_primer):
        analyze_primer.return_value = {
            "gc_content": 0.5,
            "tm": 60.0,
            "hairpin_dg": -1.0,
            "self_dimer_dg": -2.0,
        }
        response = self.client.post(
            reverse("primer_create"),
            {"primer_name": "Primer Z", "sequence": "ATCG"},
        )
        self.assertRedirects(response, reverse("primer_list"))
        response = self.client.get(reverse("primer_list"))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(list(response.context["primers"])[0].primer_name, "Primer Z")

    def test_primer_delete_permission(self):
        primer = Primer.objects.create(
            primer_name="To Delete",
            sequence="ATCG",
            length=4,
            creator=self.user,
        )
        primer.users.add(self.user)
        response = self.client.get(reverse("primer_delete", args=[primer.id]))
        self.assertRedirects(response, reverse("primer_list"))
        self.assertFalse(Primer.objects.filter(id=primer.id).exists())

        primer = Primer.objects.create(
            primer_name="Other Primer",
            sequence="ATCG",
            length=4,
            creator=self.other_user,
        )
        self.client.force_login(self.user)
        response = self.client.get(reverse("primer_delete", args=[primer.id]))
        self.assertEqual(response.status_code, 403)


class PrimerPairViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="pair_creator",
            email="pair_creator@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)
        self.forward = Primer.objects.create(
            primer_name="Forward",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        self.reverse = Primer.objects.create(
            primer_name="Reverse",
            sequence="CGATCGAT",
            length=8,
            creator=self.user,
        )
        self.forward.users.add(self.user)
        self.reverse.users.add(self.user)

    @patch("core.views.PrimerPairForm")
    def test_primerpair_create_uses_form(self, mock_form):
        pair = PrimerPair(
            name="Pair Test",
            forward_primer=self.forward,
            reverse_primer=self.reverse,
        )
        form_instance = mock_form.return_value
        form_instance.is_valid.return_value = True
        form_instance.save.return_value = pair

        response = self.client.post(reverse("primerpair_create"), {"name": "Pair Test"})

        self.assertRedirects(response, reverse("primerpair_list"))
        self.assertTrue(PrimerPair.objects.filter(name="Pair Test").exists())

    @patch("core.views.analyze_primer")
    def test_primerpair_combined_create(self, analyze_primer):
        analyze_primer.return_value = {
            "gc_content": 0.5,
            "tm": 60.0,
            "hairpin_dg": -1.0,
            "self_dimer_dg": -2.0,
        }
        response = self.client.post(
            reverse("primerpair_combined_create"),
            {
                "pair_name": "Combined Pair",
                "forward_name": "Forward",
                "forward_sequence": "ATCG",
                "reverse_name": "Reverse",
                "reverse_sequence": "CGAT",
            },
        )
        self.assertRedirects(response, reverse("primerpair_list"))
        self.assertTrue(PrimerPair.objects.filter(name="Combined Pair").exists())

    def test_primerpair_delete_permission(self):
        pair = PrimerPair.objects.create(
            name="Delete Pair",
            forward_primer=self.forward,
            reverse_primer=self.reverse,
            creator=self.user,
        )
        response = self.client.get(reverse("primerpair_delete", args=[pair.id]))
        self.assertRedirects(response, reverse("primerpair_list"))
        self.assertFalse(PrimerPair.objects.filter(id=pair.id).exists())


class MiscViewTests(TestCase):
    def setUp(self):
        super().setUp()
        self.user = User.objects.create_user(
            username="misc_user",
            email="misc@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def test_home_view(self):
        response = self.client.get(reverse("home"))
        self.assertEqual(response.status_code, 200)

    def test_register_view_get(self):
        response = self.client.get(reverse("register"))
        self.assertEqual(response.status_code, 200)

    def test_project_list_view(self):
        Project.objects.create(name="List Project", creator=self.user)
        response = self.client.get(reverse("project_list"))
        self.assertEqual(response.status_code, 200)

    def test_analyze_sequence_get(self):
        response = self.client.get(reverse("analyze_sequence"))
        self.assertEqual(response.status_code, 200)