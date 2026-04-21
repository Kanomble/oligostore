from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase, override_settings
from django.urls import reverse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from types import SimpleNamespace
from unittest.mock import patch
from io import StringIO
import json
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
from core.tasks import analyze_primerpair_products_task

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

    def test_sequencefile_upload_post_success_for_snapgene(self):
        file_data = SimpleUploadedFile("plasmid.dna", b"SnapGene", content_type="application/octet-stream")
        response = self.client.post(
            reverse("sequencefile_upload"),
            {
                "name": "SnapGene Upload",
                "file": file_data,
                "file_type": SequenceFile.FILE_SNAPGENE,
                "description": "SnapGene plasmid",
            },
        )

        self.assertRedirects(response, reverse("sequencefile_list"))
        uploaded = SequenceFile.objects.get(name="SnapGene Upload")
        self.assertEqual(uploaded.file_type, SequenceFile.FILE_SNAPGENE)

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
        SequenceFile.objects.create(
            name="Snap",
            file=SimpleUploadedFile("snap.dna", b"SnapGene"),
            file_type=SequenceFile.FILE_SNAPGENE,
            uploaded_by=self.user,
            description="SnapGene file",
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

        response = self.client.get(reverse("sequencefile_list"), {"type": "snapgene"})
        listed = list(response.context["sequence_files"])
        self.assertEqual(len(listed), 1)
        self.assertEqual(listed[0].file_type, SequenceFile.FILE_SNAPGENE)
        self.assertEqual(listed[0].name, "Snap")

    def test_sequencefile_update_type_changes_type_and_redirects_back_to_listing(self):
        sequence_file = SequenceFile.objects.create(
            name="Mutable Type",
            file=SimpleUploadedFile("mutable.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_update_type", args=[sequence_file.id]),
            {
                "file_type": SequenceFile.FILE_GENBANK,
                "next": f'{reverse("sequencefile_list")}?type=fasta&order=name',
            },
        )

        self.assertRedirects(
            response,
            f'{reverse("sequencefile_list")}?type=fasta&order=name',
        )
        sequence_file.refresh_from_db()
        self.assertEqual(sequence_file.file_type, SequenceFile.FILE_GENBANK)

    def test_sequencefile_update_type_rejects_invalid_type(self):
        sequence_file = SequenceFile.objects.create(
            name="Immutable Type",
            file=SimpleUploadedFile("immutable.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_update_type", args=[sequence_file.id]),
            {
                "file_type": "invalid",
                "next": reverse("sequencefile_list"),
            },
            follow=True,
        )

        self.assertEqual(response.status_code, 200)
        sequence_file.refresh_from_db()
        self.assertEqual(sequence_file.file_type, SequenceFile.FILE_FASTA)
        self.assertContains(response, "Invalid sequence file type.")

    def test_sequencefile_linear_view_payload(self):
        sequence_file = SequenceFile.objects.create(
            name="Linear",
            file=SimpleUploadedFile("linear.fasta", b">record1\nATCGATCGAA"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.get(
            reverse("sequencefile_linear_view", args=[sequence_file.id])
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["sequence_file"].id, sequence_file.id)
        self.assertEqual(len(response.context["records_payload"]), 1)
        self.assertEqual(response.context["records_payload"][0]["id"], "record1")
        self.assertEqual(response.context["records_payload"][0]["length"], 10)

    def test_sequencefile_linear_view_includes_requested_pcr_product(self):
        sequence_file = SequenceFile.objects.create(
            name="Linear With PCR",
            file=SimpleUploadedFile("linear_with_pcr.fasta", b">record1\nATCGATCGAA"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        product = PCRProduct.objects.create(
            name="record1_1_8",
            sequence_file=sequence_file,
            record_id="record1",
            start=1,
            end=8,
            sequence="ATCGATCG",
            creator=self.user,
        )
        product.users.add(self.user)

        response = self.client.get(
            reverse("sequencefile_linear_view", args=[sequence_file.id]),
            {"pcr_product": product.id},
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["initial_pcr_product"]["id"], product.id)
        self.assertEqual(response.context["initial_pcr_product"]["record_id"], "record1")
        self.assertEqual(response.context["initial_pcr_product"]["start"], 1)
        self.assertEqual(response.context["initial_pcr_product"]["end"], 8)

    def test_sequencefile_linear_view_ignores_inaccessible_pcr_product(self):
        own_sequence_file = SequenceFile.objects.create(
            name="Own Viewer File",
            file=SimpleUploadedFile("own_viewer.fasta", b">record1\nATCGATCGAA"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        other_user = User.objects.create_user(
            username="pcr_other",
            email="pcr_other@example.com",
            password="testpass123",
        )
        other_sequence_file = SequenceFile.objects.create(
            name="Other Viewer File",
            file=SimpleUploadedFile("other_viewer.fasta", b">record1\nATCGATCGAA"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=other_user,
        )
        hidden_product = PCRProduct.objects.create(
            name="Hidden Product",
            sequence_file=other_sequence_file,
            record_id="record1",
            start=1,
            end=8,
            sequence="ATCGATCG",
            creator=other_user,
        )
        hidden_product.users.add(other_user)

        response = self.client.get(
            reverse("sequencefile_linear_view", args=[own_sequence_file.id]),
            {"pcr_product": hidden_product.id},
        )

        self.assertEqual(response.status_code, 200)
        self.assertIsNone(response.context["initial_pcr_product"])

    def test_sequencefile_linear_view_forbidden_for_non_owner(self):
        other_user = User.objects.create_user(
            username="viewer_other",
            email="viewer_other@example.com",
            password="testpass123",
        )
        sequence_file = SequenceFile.objects.create(
            name="Other private file",
            file=SimpleUploadedFile("other.fasta", b">record1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=other_user,
        )

        response = self.client.get(
            reverse("sequencefile_linear_view", args=[sequence_file.id])
        )
        self.assertEqual(response.status_code, 404)

    def test_sequencefile_circular_view_payload(self):
        sequence_file = SequenceFile.objects.create(
            name="Circular",
            file=SimpleUploadedFile("circular.fasta", b">record1\nATCGATCGAA"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.get(
            reverse("sequencefile_circular_view", args=[sequence_file.id])
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["sequence_file"].id, sequence_file.id)
        self.assertEqual(len(response.context["records_payload"]), 1)
        self.assertEqual(response.context["records_payload"][0]["id"], "record1")
        self.assertContains(response, reverse("sequencefile_linear_view", args=[sequence_file.id]))

    def test_sequencefile_linear_record_data_includes_user_features(self):
        sequence_file = SequenceFile.objects.create(
            name="Feature Merge",
            file=SimpleUploadedFile("feature_merge.fasta", b">record1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        primer = Primer.objects.create(
            primer_name="MergedPrimer",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        primer.users.add(self.user)
        SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id="record1",
            start=2,
            end=9,
            strand=1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="MergedPrimer",
            created_by=self.user,
        )

        response = self.client.get(
            reverse("sequencefile_linear_record_data", args=[sequence_file.id]),
            {"record_index": 0, "start": 1, "end": 12},
        )

        self.assertEqual(response.status_code, 200)
        payload = response.json()
        self.assertIn("features", payload)
        self.assertTrue(any(f.get("label") == "MergedPrimer" for f in payload["features"]))
        merged = next(f for f in payload["features"] if f.get("label") == "MergedPrimer")
        self.assertEqual(merged["source"], "user")
        self.assertIn("feature_id", merged)
        self.assertEqual(merged["primer_id"], primer.id)

    @patch("core.views.sequence_files.Primer.create_with_analysis")
    def test_sequencefile_linear_create_primer_with_feature_attach(self, create_with_analysis):
        sequence_file = SequenceFile.objects.create(
            name="Attach Seq",
            file=SimpleUploadedFile("attach.fasta", b">record1\nATCGATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        def _create_primer(**kwargs):
            primer = Primer.objects.create(
                primer_name=kwargs["primer_name"],
                sequence=kwargs["sequence"],
                overhang_sequence=kwargs.get("overhang_sequence", ""),
                length=len(kwargs["sequence"]),
                gc_content=50.0,
                tm=60.0,
                hairpin_dg=-1.0,
                self_dimer_dg=-2.0,
                creator=kwargs["user"],
            )
            primer.users.add(kwargs["user"])
            return primer

        create_with_analysis.side_effect = _create_primer
        response = self.client.post(
            reverse("sequencefile_linear_create_primer", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "primer_name": "AttachPrimer",
                    "sequence": "ATCGATCG",
                    "attach_feature": True,
                    "record_id": "record1",
                    "feature_start": 3,
                    "feature_end": 10,
                    "feature_strand": 1,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        payload = response.json()
        self.assertTrue(payload["ok"])
        self.assertIsNotNone(payload["attached_feature"])
        feature = SequenceFeature.objects.get(id=payload["attached_feature"]["id"])
        self.assertEqual(feature.record_id, "record1")
        self.assertEqual(feature.start, 3)
        self.assertEqual(feature.end, 10)

    @patch("core.views.sequence_files.Primer.create_with_analysis")
    def test_sequencefile_linear_create_primer_can_attach_feature_without_saving_primer(self, create_with_analysis):
        sequence_file = SequenceFile.objects.create(
            name="Attach Only",
            file=SimpleUploadedFile("attach_only.fasta", b">record1\nATCGATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_create_primer", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "primer_name": "AttachOnlyPrimer",
                    "sequence": "ATCGATCG",
                    "save_to_primers": False,
                    "attach_feature": True,
                    "record_id": "record1",
                    "feature_start": 2,
                    "feature_end": 9,
                    "feature_strand": -1,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        payload = response.json()
        self.assertTrue(payload["ok"])
        self.assertIsNone(payload["primer"])
        feature = SequenceFeature.objects.get(id=payload["attached_feature"]["id"])
        self.assertIsNone(feature.primer)
        self.assertEqual(feature.label, "AttachOnlyPrimer")
        self.assertEqual(feature.strand, -1)
        create_with_analysis.assert_not_called()

    @patch("core.views.sequence_files.Primer.create_with_analysis")
    def test_sequencefile_linear_create_primer_requires_destination(self, create_with_analysis):
        sequence_file = SequenceFile.objects.create(
            name="No Destination",
            file=SimpleUploadedFile("no_destination.fasta", b">record1\nATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_create_primer", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "primer_name": "NowherePrimer",
                    "sequence": "ATCGATCG",
                    "save_to_primers": False,
                    "attach_feature": False,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn("Select at least one destination", response.json()["error"])
        create_with_analysis.assert_not_called()

    def test_sequencefile_linear_delete_primer_feature_only(self):
        sequence_file = SequenceFile.objects.create(
            name="Delete Feature Only",
            file=SimpleUploadedFile("delete_feature_only.fasta", b">record1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        primer = Primer.objects.create(
            primer_name="FeatureOnlyPrimer",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        primer.users.add(self.user)
        feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id="record1",
            start=1,
            end=8,
            strand=1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="FeatureOnlyPrimer",
            created_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_delete_primer", args=[sequence_file.id]),
            data=json.dumps({"feature_id": feature.id, "delete_primer": False}),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 200)
        self.assertFalse(SequenceFeature.objects.filter(id=feature.id).exists())
        self.assertTrue(Primer.objects.filter(id=primer.id).exists())

    def test_sequencefile_linear_delete_primer_feature_and_oligostore_primer(self):
        sequence_file = SequenceFile.objects.create(
            name="Delete Primer Everywhere",
            file=SimpleUploadedFile("delete_everywhere.fasta", b">record1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        primer = Primer.objects.create(
            primer_name="EverywherePrimer",
            sequence="ATCGATCG",
            length=8,
            creator=self.user,
        )
        primer.users.add(self.user)
        feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id="record1",
            start=2,
            end=9,
            strand=-1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="EverywherePrimer",
            created_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_delete_primer", args=[sequence_file.id]),
            data=json.dumps({"feature_id": feature.id, "delete_primer": True}),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 200)
        self.assertFalse(SequenceFeature.objects.filter(id=feature.id).exists())
        self.assertFalse(Primer.objects.filter(id=primer.id).exists())

    def test_sequencefile_linear_delete_primer_forbidden_for_unowned_primer(self):
        other_user = User.objects.create_user(
            username="primer_owner_other",
            email="primer_owner_other@example.com",
            password="testpass123",
        )
        sequence_file = SequenceFile.objects.create(
            name="Delete Primer Forbidden",
            file=SimpleUploadedFile("delete_forbidden.fasta", b">record1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        primer = Primer.objects.create(
            primer_name="SharedPrimer",
            sequence="ATCGATCG",
            length=8,
            creator=other_user,
        )
        primer.users.add(self.user)
        feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=primer,
            record_id="record1",
            start=3,
            end=10,
            strand=1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="SharedPrimer",
            created_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_delete_primer", args=[sequence_file.id]),
            data=json.dumps({"feature_id": feature.id, "delete_primer": True}),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 403)
        self.assertTrue(SequenceFeature.objects.filter(id=feature.id).exists())
        self.assertTrue(Primer.objects.filter(id=primer.id).exists())

    @patch("core.views.sequence_files.Primer.create_with_analysis")
    def test_sequencefile_linear_create_primer_rejects_invalid_attachment(self, create_with_analysis):
        sequence_file = SequenceFile.objects.create(
            name="Attach Invalid",
            file=SimpleUploadedFile("attach_invalid.fasta", b">record1\nATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_create_primer", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "primer_name": "AttachPrimer",
                    "sequence": "ATCGATCG",
                    "attach_feature": True,
                    "record_id": "record1",
                    "feature_start": 1,
                    "feature_end": 9999,
                    "feature_strand": 1,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn("Feature coordinates exceed record length", response.json()["error"])
        create_with_analysis.assert_not_called()

    def test_sequencefile_linear_save_pcr_product(self):
        sequence_file = SequenceFile.objects.create(
            name="PCR Save",
            file=SimpleUploadedFile("pcr_save.fasta", b">record1\nATCGATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        forward = Primer.objects.create(
            primer_name="Forward Save",
            sequence="ATCG",
            length=4,
            creator=self.user,
        )
        reverse = Primer.objects.create(
            primer_name="Reverse Save",
            sequence="GATC",
            length=4,
            creator=self.user,
        )
        forward.users.add(self.user)
        reverse.users.add(self.user)
        forward_feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=forward,
            record_id="record1",
            start=1,
            end=4,
            strand=1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="Forward Save",
            created_by=self.user,
        )
        reverse_feature = SequenceFeature.objects.create(
            sequence_file=sequence_file,
            primer=reverse,
            record_id="record1",
            start=9,
            end=12,
            strand=-1,
            feature_type=SequenceFeature.TYPE_PRIMER_BIND,
            label="Reverse Save",
            created_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_save_pcr_product", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "name": "Saved Amplicon",
                    "record_id": "record1",
                    "start": 1,
                    "end": 12,
                    "sequence": "ATCGATCGATCG",
                    "forward_primer_label": "Forward Save",
                    "reverse_primer_label": "Reverse Save",
                    "forward_primer_id": forward.id,
                    "reverse_primer_id": reverse.id,
                    "forward_feature_id": forward_feature.id,
                    "reverse_feature_id": reverse_feature.id,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        payload = response.json()
        self.assertTrue(payload["ok"])
        product = PCRProduct.objects.get(id=payload["pcr_product"]["id"])
        self.assertEqual(product.name, "Saved Amplicon")
        self.assertEqual(product.sequence_file, sequence_file)
        self.assertEqual(product.forward_primer, forward)
        self.assertEqual(product.reverse_primer, reverse)
        self.assertEqual(product.forward_feature, forward_feature)
        self.assertEqual(product.reverse_feature, reverse_feature)
        self.assertEqual(product.length, 12)

    def test_sequencefile_linear_save_pcr_product_rejects_mismatched_sequence(self):
        sequence_file = SequenceFile.objects.create(
            name="PCR Save Invalid",
            file=SimpleUploadedFile("pcr_save_invalid.fasta", b">record1\nATCGATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

        response = self.client.post(
            reverse("sequencefile_linear_save_pcr_product", args=[sequence_file.id]),
            data=json.dumps(
                {
                    "record_id": "record1",
                    "start": 1,
                    "end": 4,
                    "sequence": "TTTT",
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 400)
        self.assertIn("does not match", response.json()["error"])

    def test_pcrproduct_list_shows_only_current_user_products(self):
        own_sequence_file = SequenceFile.objects.create(
            name="Own PCR File",
            file=SimpleUploadedFile("own_pcr.fasta", b">record1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        other_user = User.objects.create_user(
            username="pcr_products_other",
            email="pcr_products_other@example.com",
            password="testpass123",
        )
        other_sequence_file = SequenceFile.objects.create(
            name="Other PCR File",
            file=SimpleUploadedFile("other_pcr.fasta", b">record1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=other_user,
        )
        own_product = PCRProduct.objects.create(
            name="Own Product",
            sequence_file=own_sequence_file,
            record_id="record1",
            start=1,
            end=4,
            sequence="ATCG",
            creator=self.user,
        )
        own_product.users.add(self.user)
        other_product = PCRProduct.objects.create(
            name="Other Product",
            sequence_file=other_sequence_file,
            record_id="record1",
            start=1,
            end=4,
            sequence="ATCG",
            creator=other_user,
        )
        other_product.users.add(other_user)

        response = self.client.get(reverse("pcrproduct_list"))

        self.assertEqual(response.status_code, 200)
        listed = list(response.context["pcr_products"])
        self.assertEqual([item.name for item in listed], ["Own Product"])

    def test_pcrproduct_list_open_viewer_link_targets_saved_product(self):
        sequence_file = SequenceFile.objects.create(
            name="PCR Link File",
            file=SimpleUploadedFile("pcr_link.fasta", b">record1\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        product = PCRProduct.objects.create(
            name="Linked Product",
            sequence_file=sequence_file,
            record_id="record1",
            start=1,
            end=4,
            sequence="ATCG",
            creator=self.user,
        )
        product.users.add(self.user)

        response = self.client.get(reverse("pcrproduct_list"))

        self.assertContains(
            response,
            f'{reverse("sequencefile_linear_view", args=[sequence_file.id])}?pcr_product={product.id}',
        )


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

    @patch("core.views.sequence_files.analyze_primer_binding")
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
            username="download_user",
            email="download@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

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

    def test_download_product_sequence_genbank_includes_source_features(self):
        record = SeqRecord(Seq("AAATTTCC"), id="seq1", name="seq1", description="template")
        record.annotations["molecule_type"] = "DNA"
        record.features.append(
            SeqFeature(
                FeatureLocation(0, 6, strand=1),
                type="promoter",
                qualifiers={"label": ["source_feature"]},
            )
        )
        handle = StringIO()
        SeqIO.write(record, handle, "genbank")
        sequence_file = SequenceFile.objects.create(
            name="Template",
            file=SimpleUploadedFile("template.gb", handle.getvalue().encode("utf-8")),
            file_type=SequenceFile.FILE_GENBANK,
            uploaded_by=self.user,
        )
        sequence_file.users.add(self.user)
        SequenceFeature.objects.create(
            sequence_file=sequence_file,
            record_id="seq1",
            start=2,
            end=5,
            strand=1,
            feature_type=SequenceFeature.TYPE_CUSTOM,
            label="user_feature",
            created_by=self.user,
        )

        response = self.client.post(
            reverse("download_product_sequence"),
            {
                "export_format": "genbank",
                "pair_index": "7",
                "sequence_file_id": sequence_file.id,
                "record_id": "seq1",
                "product_start": 1,
                "product_end": 6,
                "wraps_origin": "false",
            },
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(
            response["Content-Disposition"],
            "attachment; filename=pcr_product_pair_7.gb",
        )
        exported_record = SeqIO.read(StringIO(response.content.decode("utf-8")), "genbank")
        self.assertEqual(str(exported_record.seq), "AAATTT")
        labels = {
            (feature.qualifiers.get("label") or [""])[0]
            for feature in exported_record.features
        }
        self.assertIn("source_feature", labels)
        self.assertIn("user_feature", labels)


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

    @patch("core.views.analysis.analyze_primer")
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

    @patch("core.views.analysis.analyze_cross_dimer")
    @patch("core.views.analysis.analyze_primer")
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
        self.sequence_file = SequenceFile.objects.create(
            name="Project Sequence",
            file=SimpleUploadedFile("project.fasta", b">seq\nATCG"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
            description="Project file description",
        )
        self.pcr_product = PCRProduct.objects.create(
            name="Project Product",
            sequence_file=self.sequence_file,
            record_id="seq",
            start=1,
            end=4,
            sequence="ATCG",
            creator=self.user,
        )
        self.pcr_product.users.add(self.user)

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

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

    def test_project_add_and_remove_sequencefile(self):
        response = self.client.get(
            reverse(
                "project_add_sequencefile",
                args=[self.project.id, self.sequence_file.id],
            )
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertIn(self.sequence_file, self.project.sequence_files.all())

        response = self.client.get(
            reverse(
                "project_remove_sequencefile",
                args=[self.project.id, self.sequence_file.id],
            )
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertNotIn(self.sequence_file, self.project.sequence_files.all())

    def test_project_add_and_remove_pcr_product(self):
        response = self.client.get(
            reverse(
                "project_add_pcr_product",
                args=[self.project.id, self.pcr_product.id],
            )
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertIn(self.pcr_product, self.project.pcr_products.all())

        response = self.client.get(
            reverse(
                "project_remove_pcr_product",
                args=[self.project.id, self.pcr_product.id],
            )
        )
        self.assertRedirects(
            response,
            reverse("project_dashboard", args=[self.project.id]),
        )
        self.assertNotIn(self.pcr_product, self.project.pcr_products.all())

    def test_project_dashboard_contains_sequence_view_links(self):
        self.project.sequence_files.add(self.sequence_file)
        self.project.pcr_products.add(self.pcr_product)

        response = self.client.get(reverse("project_dashboard", args=[self.project.id]))

        self.assertContains(
            response,
            reverse("sequencefile_linear_view", args=[self.sequence_file.id]),
        )
        self.assertContains(
            response,
            f'{reverse("sequencefile_linear_view", args=[self.sequence_file.id])}?pcr_product={self.pcr_product.id}',
        )
        self.assertContains(
            response,
            reverse("sequencefile_circular_view", args=[self.sequence_file.id]),
        )

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

    @patch("core.views.primers.Primer.create_with_analysis")
    def test_primer_create_and_list(self, create_with_analysis):
        def _create_primer(**kwargs):
            primer = Primer.objects.create(
                primer_name=kwargs["primer_name"],
                sequence=kwargs["sequence"],
                overhang_sequence=kwargs.get("overhang_sequence", ""),
                length=len(kwargs["sequence"]),
                gc_content=50.0,
                tm=60.0,
                hairpin_dg=-1.0,
                self_dimer_dg=-2.0,
                creator=kwargs["user"],
            )
            primer.users.add(kwargs["user"])
            return primer

        create_with_analysis.side_effect = _create_primer
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

    @patch("core.views.primerpairs.PrimerPairForm")
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

        self.assertEqual(response.status_code, 200)
        self.assertTrue(PrimerPair.objects.filter(name="Pair Test").exists())

    @patch("core.views.primerpairs.Primer.create_with_analysis")
    def test_primerpair_combined_create(self, create_with_analysis):
        created_primers = []

        def _create_primer(**kwargs):
            primer = Primer.objects.create(
                primer_name=kwargs["primer_name"],
                sequence=kwargs["sequence"],
                overhang_sequence=kwargs.get("overhang_sequence", ""),
                length=len(kwargs["sequence"]),
                gc_content=50.0,
                tm=60.0,
                hairpin_dg=-1.0,
                self_dimer_dg=-2.0,
                creator=kwargs["user"],
            )
            primer.users.add(kwargs["user"])
            created_primers.append(primer)
            return primer

        create_with_analysis.side_effect = _create_primer
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
        self.assertEqual(response.status_code, 200)
        self.assertEqual(create_with_analysis.call_count, 2)
        pair = PrimerPair.objects.get(name="Combined Pair")
        self.assertEqual(pair.forward_primer, created_primers[0])
        self.assertEqual(pair.reverse_primer, created_primers[1])

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


class PCRProductDiscoveryViewTests(TestCase):
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
            username="pcr_products_user",
            email="pcr_products@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)
        self.forward = Primer.objects.create(
            primer_name="Forward PCR",
            sequence="AAA",
            length=3,
            creator=self.user,
        )
        self.reverse = Primer.objects.create(
            primer_name="Reverse PCR",
            sequence="AAA",
            length=3,
            creator=self.user,
        )
        self.forward.users.add(self.user)
        self.reverse.users.add(self.user)
        self.pair = PrimerPair.objects.create(
            name="PCR Pair",
            forward_primer=self.forward,
            reverse_primer=self.reverse,
            creator=self.user,
        )
        self.pair.users.add(self.user)
        self.sequence_file = SequenceFile.objects.create(
            name="PCR Sequence",
            file=SimpleUploadedFile("pcr.fasta", b">seq1\nAAATTT"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_pcr_product_discovery_get(self):
        response = self.client.get(reverse("primerpair_products"))
        self.assertEqual(response.status_code, 200)
        self.assertIn("form", response.context)

    @patch("core.views.primerpairs.analyze_primerpair_products_task.delay")
    def test_pcr_product_discovery_async_post_starts_task(self, delay_mock):
        delay_mock.return_value = SimpleNamespace(id="task-123")
        response = self.client.post(
            reverse("primerpair_products_async"),
            {
                "primer_pair": self.pair.id,
                "sequence_file": self.sequence_file.id,
                "max_mismatches": 0,
                "block_3prime_mismatch": "on",
            },
        )
        self.assertEqual(response.status_code, 202)
        self.assertEqual(response.json()["task_id"], "task-123")
        self.assertIn("job_id", response.json())
        job = AnalysisJob.objects.get(id=response.json()["job_id"])
        self.assertEqual(job.owner, self.user)
        self.assertEqual(job.status, AnalysisJob.STATUS_RUNNING)
        delay_mock.assert_called_once_with(
            analysis_job_id=job.id,
            primer_pair_id=self.pair.id,
            sequence_file_id=self.sequence_file.id,
            max_mismatches=0,
            block_3prime_mismatch=True,
        )

    @patch("core.tasks.analyze_primerpair_products")
    @patch("core.tasks.mark_job_success")
    def test_pcr_product_discovery_task_uses_binding_sequences_without_overhangs(
        self,
        mark_job_success_mock,
        analyze_products_mock,
    ):
        self.forward.overhang_sequence = "GGATCC"
        self.forward.save(update_fields=["overhang_sequence"])
        self.reverse.overhang_sequence = "GACT"
        self.reverse.save(update_fields=["overhang_sequence"])
        analyze_products_mock.return_value = []

        payload = analyze_primerpair_products_task.run(
            analysis_job_id=123,
            primer_pair_id=self.pair.id,
            sequence_file_id=self.sequence_file.id,
            max_mismatches=0,
            block_3prime_mismatch=True,
        )

        analyze_products_mock.assert_called_once_with(
            forward_primer_sequence="AAA",
            reverse_primer_sequence="AAA",
            forward_overhang_sequence="GGATCC",
            reverse_overhang_sequence="GACT",
            sequence_file=self.sequence_file,
            max_mismatches=0,
            block_3prime_mismatch=True,
        )
        self.assertEqual(payload["products"], [])
        mark_job_success_mock.assert_called_once()

    def test_pcr_product_discovery_save_creates_pcr_product(self):
        response = self.client.post(
            reverse("primerpair_products_save"),
            data=json.dumps(
                {
                    "primer_pair_id": self.pair.id,
                    "sequence_file_id": self.sequence_file.id,
                    "name": "Saved Finder Product",
                    "record_id": "seq1",
                    "product_start": 1,
                    "product_end": 6,
                    "product_sequence": "AAATTT",
                    "wraps_origin": False,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        saved = PCRProduct.objects.get(name="Saved Finder Product")
        self.assertEqual(saved.sequence_file, self.sequence_file)
        self.assertEqual(saved.record_id, "seq1")
        self.assertEqual(saved.start, 1)
        self.assertEqual(saved.end, 6)
        self.assertEqual(saved.sequence, "AAATTT")
        self.assertEqual(saved.forward_primer, self.forward)
        self.assertEqual(saved.reverse_primer, self.reverse)

    def test_pcr_product_discovery_save_accepts_overhang_enriched_sequence(self):
        self.forward.overhang_sequence = "GGATCC"
        self.forward.save(update_fields=["overhang_sequence"])
        self.reverse.overhang_sequence = "GACT"
        self.reverse.save(update_fields=["overhang_sequence"])

        response = self.client.post(
            reverse("primerpair_products_save"),
            data=json.dumps(
                {
                    "primer_pair_id": self.pair.id,
                    "sequence_file_id": self.sequence_file.id,
                    "record_id": "seq1",
                    "product_start": 1,
                    "product_end": 6,
                    "product_sequence": "GGATCCAAATTTAGTC",
                    "wraps_origin": False,
                }
            ),
            content_type="application/json",
        )

        self.assertEqual(response.status_code, 201)
        saved = PCRProduct.objects.get(id=response.json()["pcr_product"]["id"])
        self.assertEqual(saved.sequence, "GGATCCAAATTTAGTC")
        self.assertEqual(saved.length, 16)


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
        self.assertContains(response, reverse("pcrproduct_list"))

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


class CloningViewTests(TestCase):
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
            username="cloning_user",
            email="cloning@example.com",
            password="testpass123",
        )
        self.client.force_login(self.user)

    def tearDown(self):
        self._media_override.disable()
        super().tearDown()

    def test_cloning_asset_list_shows_sequence_files_and_pcr_products(self):
        sequence_file = SequenceFile.objects.create(
            name="Vector Source",
            file=SimpleUploadedFile("vector.fasta", b">vec1\nGAATTCACCCCGGATCC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        sequence_file.users.add(self.user)
        pcr_product = PCRProduct.objects.create(
            name="Insert Product",
            sequence_file=sequence_file,
            record_id="vec1",
            start=1,
            end=8,
            sequence="ATGCATGC",
            creator=self.user,
        )
        pcr_product.users.add(self.user)

        response = self.client.get(reverse("cloning_asset_list"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Vector Source")
        self.assertContains(response, "Insert Product")

    def test_cloning_construct_create_and_detail(self):
        vector = SequenceFile.objects.create(
            name="Vector Backbone",
            file=SimpleUploadedFile("vector_backbone.fasta", b">vec1\nGAATTCACCCCGGATCC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        vector.users.add(self.user)
        insert_source = SequenceFile.objects.create(
            name="Insert Fragment",
            file=SimpleUploadedFile("insert_fragment.fasta", b">ins1\nATGCATGC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        insert_source.users.add(self.user)
        insert = PCRProduct.objects.create(
            name="Insert Product",
            sequence_file=insert_source,
            record_id="ins1",
            start=1,
            end=8,
            sequence="ATGCATGC",
            creator=self.user,
        )
        insert.users.add(self.user)

        response = self.client.post(
            reverse("cloning_construct_create"),
            {
                "name": "Construct One",
                "description": "Restriction ligation construct",
                "vector_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{vector.id}",
                "insert_asset": f"{CloningConstruct.SOURCE_PCR_PRODUCT}:{insert.id}",
                "assembly_strategy": CloningConstruct.STRATEGY_RESTRICTION_LIGATION,
                "left_enzyme": "EcoRI",
                "right_enzyme": "BamHI",
            },
        )

        construct = CloningConstruct.objects.get(name="Construct One")
        self.assertRedirects(
            response,
            reverse("cloning_construct_detail", args=[construct.id]),
        )
        self.assertTrue(construct.is_valid)
        self.assertIn("ATGCATGC", construct.assembled_sequence)
        self.assertEqual(construct.vector_sequence_file, vector)
        self.assertEqual(construct.insert_pcr_product, insert)

    def test_cloning_construct_rejects_same_asset_for_vector_and_insert(self):
        vector = SequenceFile.objects.create(
            name="Shared Asset",
            file=SimpleUploadedFile("shared_asset.fasta", b">vec1\nGAATTCACCCCGGATCC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        vector.users.add(self.user)

        response = self.client.post(
            reverse("cloning_construct_create"),
            {
                "name": "Construct Invalid",
                "description": "Invalid restriction ligation construct",
                "vector_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{vector.id}",
                "insert_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{vector.id}",
                "assembly_strategy": CloningConstruct.STRATEGY_RESTRICTION_LIGATION,
                "left_enzyme": "EcoRI",
                "right_enzyme": "BamHI",
            },
        )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Vector and insert must be different assets.")
        self.assertFalse(CloningConstruct.objects.filter(name="Construct Invalid").exists())

    def test_cloning_construct_shows_parse_error_for_bad_sequence_file(self):
        vector = SequenceFile.objects.create(
            name="Bad Vector",
            file=SimpleUploadedFile("bad_vector.fasta", b"not-a-valid-sequence-file"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        vector.users.add(self.user)
        insert_source = SequenceFile.objects.create(
            name="Insert Source",
            file=SimpleUploadedFile("insert_ok.fasta", b">ins1\nATGCATGC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        insert_source.users.add(self.user)

        response = self.client.post(
            reverse("cloning_construct_create"),
            {
                "name": "Construct Parse Error",
                "description": "Should fail on parse",
                "vector_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{vector.id}",
                "insert_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{insert_source.id}",
                "assembly_strategy": CloningConstruct.STRATEGY_RESTRICTION_LIGATION,
                "left_enzyme": "EcoRI",
                "right_enzyme": "BamHI",
            },
        )

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "could not be parsed for cloning")
        self.assertFalse(CloningConstruct.objects.filter(name="Construct Parse Error").exists())

    def test_cloning_construct_allows_single_enzyme_mode(self):
        vector = SequenceFile.objects.create(
            name="Single Enzyme Vector",
            file=SimpleUploadedFile("single_enzyme_vector.fasta", b">vec1\nAAAAGAATTCTTTT"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        vector.users.add(self.user)
        insert_source = SequenceFile.objects.create(
            name="Insert Source",
            file=SimpleUploadedFile("single_enzyme_insert.fasta", b">ins1\nATGCATGC"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        insert_source.users.add(self.user)

        response = self.client.post(
            reverse("cloning_construct_create"),
            {
                "name": "Single Enzyme Construct",
                "description": "Single enzyme ligation",
                "vector_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{vector.id}",
                "insert_asset": f"{CloningConstruct.SOURCE_SEQUENCE_FILE}:{insert_source.id}",
                "assembly_strategy": CloningConstruct.STRATEGY_RESTRICTION_LIGATION,
                "left_enzyme": "EcoRI",
                "right_enzyme": "EcoRI",
            },
        )

        construct = CloningConstruct.objects.get(name="Single Enzyme Construct")
        self.assertRedirects(
            response,
            reverse("cloning_construct_detail", args=[construct.id]),
        )
        self.assertTrue(construct.is_valid)
        self.assertIn("ATGCATGC", construct.assembled_sequence)

    def test_cloning_construct_form_filters_enzymes_to_detected_sites(self):
        sequence_file = SequenceFile.objects.create(
            name="EcoRI Only",
            file=SimpleUploadedFile("ecori_only.fasta", b">vec1\nAAAAGAATTCTTTT"),
            file_type=SequenceFile.FILE_FASTA,
            uploaded_by=self.user,
        )
        sequence_file.users.add(self.user)

        response = self.client.get(reverse("cloning_construct_create"))

        self.assertEqual(response.status_code, 200)
        form = response.context["form"]
        left_choices = {value for value, _ in form.fields["left_enzyme"].choices}
        self.assertIn("EcoRI", left_choices)
