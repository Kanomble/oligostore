import tempfile
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from django.test import SimpleTestCase

from core.services import primer_analysis
from core.services import primer_binding
from core.services import sequence_loader
from core.services import user_assignment
from core.services import cloning
from core.services import cloning_exports


class PrimerAnalysisTests(SimpleTestCase):
    def test_sanitize_sequence_normalizes_and_validates(self):
        cleaned = primer_analysis.sanitize_sequence(" aCg t\nN ")
        self.assertEqual(cleaned, "ACGTN")

        with self.assertRaises(ValueError) as exc:
            primer_analysis.sanitize_sequence("ACGTX")
        self.assertIn("X", str(exc.exception))

    def test_reverse_complement(self):
        self.assertEqual(primer_analysis.reverse_complement("ATGC"), "GCAT")

    def test_find_binding_site(self):
        self.assertEqual(primer_analysis.find_binding_site("AAACCC", "ACC"), 2)
        self.assertIsNone(primer_analysis.find_binding_site("AAACCC", "GGG"))

    def test_render_binding_line(self):
        self.assertEqual(
            primer_analysis.render_binding_line("AAACCC", "ACC", 2),
            "..ACC.",
        )
        self.assertEqual(
            primer_analysis.render_binding_line("AAACCC", "ACC", None),
            "No binding site found.",
        )

    def test_window_sequence_and_render_windowed_line(self):
        window, primer_start_in_window, primer_len = primer_analysis.window_sequence(
            "ACGTACGTAC", pos=3, primer_len=3, flank=2
        )
        self.assertEqual(window, "...CGTACGT...")
        self.assertEqual(primer_start_in_window, 5)
        self.assertEqual(primer_len, 3)
        self.assertEqual(
            primer_analysis.render_windowed_line(window, primer_start_in_window, primer_len),
            ".....TAC.....",
        )

        window, primer_start, primer_len = primer_analysis.window_sequence(
            "ACGTACGTAC", pos=0, primer_len=3, flank=2
        )
        self.assertEqual(window, "ACGTA...")
        self.assertEqual(primer_start, 0)
        self.assertEqual(primer_len, 3)

    def test_highlight_binding(self):
        result = primer_analysis.highlight_binding("AACCGG", 2, 2)
        self.assertIn("<span", result)
        self.assertIn("CC", result)

    def test_analyze_sequence_pair_left_right_none(self):
        primer3_result = {
            "PRIMER_LEFT_NUM_RETURNED": 1,
            "PRIMER_RIGHT_NUM_RETURNED": 1,
            "PRIMER_PAIR_NUM_RETURNED": 1,
            "PRIMER_LEFT_0_SEQUENCE": "AAA",
            "PRIMER_RIGHT_0_SEQUENCE": "TTT",
            "PRIMER_LEFT_0_TM": 60.0,
            "PRIMER_RIGHT_0_TM": 61.0,
            "PRIMER_PAIR_0_PRODUCT_SIZE": 100,
            "PRIMER_PAIR_0_PENALTY": 0.1,
        }
        with mock.patch.object(
            primer_analysis.primer3.bindings,
            "design_primers",
            return_value=primer3_result,
        ):
            primer_list, results, mode = primer_analysis.analyze_sequence(
                "ACGT", {"PRIMER_NUM_RETURN": 1}
            )
        self.assertEqual(mode, "PAIR")
        self.assertEqual(len(primer_list), 1)
        self.assertEqual(primer_list[0]["left_seq"], "AAA")
        self.assertEqual(results, primer3_result)

        left_only = {
            "PRIMER_LEFT_NUM_RETURNED": 1,
            "PRIMER_RIGHT_NUM_RETURNED": 0,
            "PRIMER_PAIR_NUM_RETURNED": 0,
            "PRIMER_LEFT_0_SEQUENCE": "AAA",
            "PRIMER_LEFT_0_TM": 59.0,
            "PRIMER_LEFT_0_PENALTY": 1.0,
        }
        with mock.patch.object(
            primer_analysis.primer3.bindings,
            "design_primers",
            return_value=left_only,
        ):
            primer_list, _, mode = primer_analysis.analyze_sequence(
                "ACGT", {"PRIMER_NUM_RETURN": 1}
            )
        self.assertEqual(mode, "LEFT")
        self.assertEqual(primer_list[0]["seq"], "AAA")

        right_only = {
            "PRIMER_LEFT_NUM_RETURNED": 0,
            "PRIMER_RIGHT_NUM_RETURNED": 1,
            "PRIMER_PAIR_NUM_RETURNED": 0,
            "PRIMER_RIGHT_0_SEQUENCE": "TTT",
            "PRIMER_RIGHT_0_TM": 58.0,
            "PRIMER_RIGHT_0_PENALTY": 2.0,
        }
        with mock.patch.object(
            primer_analysis.primer3.bindings,
            "design_primers",
            return_value=right_only,
        ):
            primer_list, _, mode = primer_analysis.analyze_sequence(
                "ACGT", {"PRIMER_NUM_RETURN": 1}
            )
        self.assertEqual(mode, "RIGHT")
        self.assertEqual(primer_list[0]["seq"], "TTT")

        none_result = {
            "PRIMER_LEFT_NUM_RETURNED": 0,
            "PRIMER_RIGHT_NUM_RETURNED": 0,
            "PRIMER_PAIR_NUM_RETURNED": 0,
        }
        with mock.patch.object(
            primer_analysis.primer3.bindings,
            "design_primers",
            return_value=none_result,
        ):
            primer_list, _, mode = primer_analysis.analyze_sequence(
                "ACGT", {"PRIMER_NUM_RETURN": 1}
            )
        self.assertEqual(mode, "NONE")
        self.assertEqual(primer_list, [])

    def test_analyze_primer_and_cross_dimer(self):
        hairpin = SimpleNamespace(dg=-5000, structure_found=True)
        dimer = SimpleNamespace(dg=-2500, structure_found=False)
        with mock.patch.object(primer_analysis.primer3, "calc_hairpin", return_value=hairpin), mock.patch.object(
            primer_analysis.primer3, "calc_homodimer", return_value=dimer
        ), mock.patch.object(primer_analysis.primer3, "calcTm", return_value=60.1234):
            result = primer_analysis.analyze_primer("acgt")
        self.assertEqual(result["gc_content"], 0.5)
        self.assertEqual(result["tm"], 60.12)
        self.assertEqual(result["hairpin_dg"], -5.0)
        self.assertEqual(result["self_dimer_dg"], -2.5)
        self.assertTrue(result["hairpin"])
        self.assertFalse(result["self_dimer"])

        with mock.patch.object(
            primer_analysis.primer3,
            "calcHeterodimer",
            return_value="hetero",
        ):
            self.assertEqual(
                primer_analysis.analyze_cross_dimer("AAA", "TTT"),
                "hetero",
            )

    def test_analyze_primer_truncates_long_sequence_for_thermo(self):
        hairpin = SimpleNamespace(dg=-5000, structure_found=True)
        dimer = SimpleNamespace(dg=-2500, structure_found=False)
        long_seq = "A" * 61
        with mock.patch.object(
            primer_analysis.primer3,
            "calc_hairpin",
            return_value=hairpin,
        ) as hairpin_mock, mock.patch.object(
            primer_analysis.primer3,
            "calc_homodimer",
            return_value=dimer,
        ) as dimer_mock, mock.patch.object(
            primer_analysis.primer3,
            "calcTm",
            return_value=60.1234,
        ):
            primer_analysis.analyze_primer(long_seq)

        hairpin_mock.assert_called_once_with(long_seq[-60:])
        dimer_mock.assert_called_once_with(long_seq[-60:])

class PrimerBindingTests(SimpleTestCase):
    def test_reverse_complement(self):
        self.assertEqual(primer_binding.reverse_complement("ATGC"), "GCAT")

    def test_count_mismatches(self):
        iter_mismatches = primer_binding.iter_mismatch_counts("AAAA", "AAAT")
        mismatch = 0
        for i, mismatch in enumerate(iter_mismatches):
            mismatch = mismatch
        self.assertEqual(mismatch, 1)

    def test_scan_sequence_respects_3prime_block(self):
        hits = primer_binding.scan_sequence(
            "AAAC", "AAT", strand="+", max_mismatches=1, block_3prime_mismatch=True
        )
        self.assertEqual(hits, [])

        hits = primer_binding.scan_sequence(
            "AAAC", "AAT", strand="+", max_mismatches=1, block_3prime_mismatch=False
        )
        self.assertEqual(len(hits), 2)
        self.assertEqual(hits[0].start, 0)

    def test_analyze_primer_binding(self):
        record = SeqRecord(Seq("AAACCC"), id="seq1", description="")
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta") as handle:
            SeqIO.write(record, handle, "fasta")
            handle.flush()
            sequence_file = SimpleNamespace(
                file=SimpleNamespace(path=handle.name),
                file_type="fasta",
            )
            results = primer_binding.analyze_primer_binding("AAA", sequence_file)

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].record_id, "seq1")
        self.assertEqual(results[0].strand, "+")
        self.assertEqual(results[0].start, 0)
        self.assertEqual(results[0].end, 3)

    def test_analyze_primerpair_products(self):
        record = SeqRecord(Seq("AAATTT"), id="seq1", description="")
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta") as handle:
            SeqIO.write(record, handle, "fasta")
            handle.flush()
            sequence_file = SimpleNamespace(
                file=SimpleNamespace(path=handle.name),
                file_type="fasta",
            )
            results = primer_binding.analyze_primerpair_products(
                forward_primer_sequence="AAA",
                reverse_primer_sequence="AAA",
                sequence_file=sequence_file,
            )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].record_id, "seq1")
        self.assertEqual(results[0].forward_start, 1)
        self.assertEqual(results[0].reverse_end, 6)
        self.assertEqual(results[0].product_sequence, "AAATTT")
        self.assertFalse(results[0].wraps_origin)
        self.assertFalse(results[0].is_circular_record)

    def test_analyze_primerpair_products_includes_overhangs_in_product_sequence(self):
        record = SeqRecord(Seq("AAATTT"), id="seq1", description="")
        sequence_file = SimpleNamespace(
            file=SimpleNamespace(path="ignored"),
            file_type="fasta",
        )

        with mock.patch.object(
            primer_binding,
            "load_sequences_from_sequence_file",
            return_value=iter([record]),
        ):
            results = primer_binding.analyze_primerpair_products(
                forward_primer_sequence="AAA",
                reverse_primer_sequence="AAA",
                forward_overhang_sequence="GGATCC",
                reverse_overhang_sequence="GACT",
                sequence_file=sequence_file,
            )

        self.assertEqual(len(results), 1)
        self.assertEqual(results[0].product_sequence, "GGATCCAAATTTAGTC")
        self.assertEqual(results[0].product_length, 16)

    def test_analyze_primerpair_products_skips_wraparound_product_on_linear_record(self):
        record = SeqRecord(Seq("TTTAAA"), id="seq1", description="")
        sequence_file = SimpleNamespace(
            file=SimpleNamespace(path="ignored"),
            file_type="fasta",
        )

        with mock.patch.object(
            primer_binding,
            "load_sequences_from_sequence_file",
            return_value=iter([record]),
        ):
            results = primer_binding.analyze_primerpair_products(
                forward_primer_sequence="AAA",
                reverse_primer_sequence="AAA",
                sequence_file=sequence_file,
            )

        self.assertEqual(results, [])

    def test_analyze_primerpair_products_allows_wraparound_product_on_circular_record(self):
        record = SeqRecord(Seq("TTTAAA"), id="plasmid", description="")
        record.annotations["topology"] = "circular"
        sequence_file = SimpleNamespace(
            file=SimpleNamespace(path="ignored"),
            file_type="snapgene",
        )

        with mock.patch.object(
            primer_binding,
            "load_sequences_from_sequence_file",
            return_value=iter([record]),
        ):
            results = primer_binding.analyze_primerpair_products(
                forward_primer_sequence="AAA",
                reverse_primer_sequence="AAA",
                sequence_file=sequence_file,
            )

        self.assertEqual(len(results), 1)
        self.assertTrue(results[0].is_circular_record)
        self.assertTrue(results[0].wraps_origin)
        self.assertEqual(results[0].forward_start, 4)
        self.assertEqual(results[0].reverse_end, 3)
        self.assertEqual(results[0].product_length, 6)
        self.assertEqual(results[0].product_sequence, "AAATTT")

    def test_binding_services_use_sequence_file_loader_helper(self):
        record = SeqRecord(Seq("AAATTT"), id="seq1", description="")
        sequence_file = SimpleNamespace(
            file=SimpleNamespace(path="ignored"),
            file_type="snapgene",
        )

        with mock.patch.object(
            primer_binding,
            "load_sequences_from_sequence_file",
            return_value=iter([record]),
        ) as loader_mock:
            binding_hits = primer_binding.analyze_primer_binding("AAA", sequence_file)
            product_hits = primer_binding.analyze_primerpair_products(
                forward_primer_sequence="AAA",
                reverse_primer_sequence="AAA",
                sequence_file=sequence_file,
            )

        self.assertEqual(len(binding_hits), 1)
        self.assertEqual(len(product_hits), 1)
        self.assertEqual(loader_mock.call_count, 2)
        loader_mock.assert_any_call(sequence_file)


class SequenceLoaderTests(SimpleTestCase):
    def test_load_sequences_fasta_and_genbank(self):
        record = SeqRecord(Seq("AAACCC"), id="seq1", description="")
        record.annotations["molecule_type"] = "DNA"

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".fasta") as fasta_handle:
            SeqIO.write(record, fasta_handle, "fasta")
            fasta_handle.flush()
            records = list(sequence_loader.load_sequences(fasta_handle.name, "fasta"))
            self.assertEqual(records[0].id, "seq1")

        with tempfile.NamedTemporaryFile(mode="w+", suffix=".gb") as gb_handle:
            SeqIO.write(record, gb_handle, "genbank")
            gb_handle.flush()
            records = list(sequence_loader.load_sequences(gb_handle.name, "genbank"))
            self.assertEqual(records[0].id, "seq1")

        with self.assertRaises(ValueError):
            list(sequence_loader.load_sequences("fake.txt", "txt"))

    def test_load_sequences_snapgene_uses_biopython_snapgene_parser(self):
        with mock.patch.object(sequence_loader.SeqIO, "parse", return_value=iter(())) as parse_mock:
            records = list(sequence_loader.load_sequences("example.dna", "snapgene"))

        self.assertEqual(records, [])
        parse_mock.assert_called_once_with("example.dna", "snapgene")

    def test_load_sequences_normalizes_file_type_before_parser_lookup(self):
        with mock.patch.object(sequence_loader.SeqIO, "parse", return_value=iter(())) as parse_mock:
            records = list(sequence_loader.load_sequences("example.dna", " SnapGene "))

        self.assertEqual(records, [])
        parse_mock.assert_called_once_with("example.dna", "snapgene")

    def test_load_sequences_from_sequence_file_uses_sequence_file_metadata(self):
        sequence_file = SimpleNamespace(
            file=SimpleNamespace(path="example.dna"),
            file_type="snapgene",
        )

        with mock.patch.object(sequence_loader.SeqIO, "parse", return_value=iter(())) as parse_mock:
            records = list(sequence_loader.load_sequences_from_sequence_file(sequence_file))

        self.assertEqual(records, [])
        parse_mock.assert_called_once_with("example.dna", "snapgene")


class UserAssignmentTests(SimpleTestCase):
    def test_assign_creator_updates_object(self):
        class Users:
            def __init__(self):
                self.added = None

            def add(self, user):
                self.added = user

        class Obj:
            def __init__(self):
                self.creator = None
                self.saved = False
                self.users = Users()

            def save(self):
                self.saved = True

        obj = Obj()
        user = "user-1"
        returned = user_assignment.assign_creator(obj, user)

        self.assertIs(returned, obj)
        self.assertEqual(obj.creator, user)
        self.assertFalse(obj.saved)
        self.assertIsNone(obj.users.added)

    def test_grant_user_access_adds_user(self):
        class Users:
            def __init__(self):
                self.added = None

            def add(self, user):
                self.added = user

        class Obj:
            def __init__(self):
                self.users = Users()

        obj = Obj()
        user = "user-1"
        returned = user_assignment.grant_user_access(obj, user)

        self.assertIs(returned, obj)
        self.assertEqual(obj.users.added, user)


class CloningServiceTests(SimpleTestCase):
    def test_asset_choice_helpers_encode_record_aware_values(self):
        sequence_choice = cloning.build_sequence_file_asset_choice(
            sequence_file_id=7,
            record_id="plasmidA",
        )
        product_choice = cloning.build_pcr_product_asset_choice(pcr_product_id=12)

        self.assertEqual(sequence_choice.encoded_value, "sequence_file:7:plasmidA")
        self.assertEqual(product_choice.encoded_value, "pcr_product:12")

    def test_template_asset_choices_support_bbsi_fragment_selection(self):
        template_choices = cloning.get_template_asset_choices()
        lvl25_choice = None
        tprom_choice = None
        vector_asset = None
        insert_asset = None
        vector_fragments = []
        insert_fragments = []

        for value, _ in template_choices:
            if lvl25_choice is None and value.startswith("template:LvL25_without_J23100.gb:"):
                candidate = cloning.resolve_asset_choice(user=None, choice=value)
                candidate_fragments = cloning.build_digest_fragment_choices(
                    sequence=candidate.sequence,
                    enzyme_name="BbsI",
                )
                if candidate_fragments:
                    lvl25_choice = value
                    vector_asset = candidate
                    vector_fragments = candidate_fragments
            if tprom_choice is None and value.startswith("template:tProm_leader_sequence_CRISPR_1.gb:"):
                candidate = cloning.resolve_asset_choice(user=None, choice=value)
                candidate_fragments = cloning.build_digest_fragment_choices(
                    sequence=candidate.sequence,
                    enzyme_name="BbsI",
                )
                if candidate_fragments:
                    tprom_choice = value
                    insert_asset = candidate
                    insert_fragments = candidate_fragments
            if lvl25_choice and tprom_choice:
                break

        self.assertIsNotNone(lvl25_choice)
        self.assertIsNotNone(tprom_choice)
        self.assertIsNotNone(vector_asset)
        self.assertIsNotNone(insert_asset)

        self.assertTrue(vector_fragments)
        self.assertTrue(insert_fragments)

        preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="BbsI",
            right_enzyme="BbsI",
            vector_fragment_index=vector_fragments[-1][0],
            insert_fragment_index=insert_fragments[-1][0],
        )

        self.assertTrue(preview_data.is_valid)
        self.assertLess(preview_data.assembled_length, 312)

    def test_preview_cloning_construct_supports_fragment_selection_for_other_same_enzyme(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Vector",
            sequence="AAAAGAATTCTTTT",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="TTTTGAATTCAAAA",
        )

        enzyme = cloning._get_enzyme_by_name("EcoRI")
        self.assertIsNotNone(enzyme)
        vector_fragments = cloning._digest_sequence_fragments(vector_asset.sequence, enzyme)
        insert_fragments = cloning._digest_sequence_fragments(insert_asset.sequence, enzyme)
        self.assertGreater(len(vector_fragments), 1)
        self.assertGreater(len(insert_fragments), 1)

        selected_vector_fragment = vector_fragments[-1]
        selected_insert_fragment = insert_fragments[-1]

        preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="EcoRI",
            right_enzyme="EcoRI",
            vector_fragment_index=selected_vector_fragment.index,
            insert_fragment_index=selected_insert_fragment.index,
        )

        self.assertTrue(preview_data.is_valid)
        self.assertEqual(preview_data.vector_fragment_index, selected_vector_fragment.index)
        self.assertEqual(preview_data.insert_fragment_index, selected_insert_fragment.index)
        self.assertEqual(preview_data.vector_fragment_start, selected_vector_fragment.start)
        self.assertEqual(preview_data.insert_fragment_start, selected_insert_fragment.start)
        self.assertEqual(
            preview_data.assembled_sequence,
            selected_vector_fragment.sequence + selected_insert_fragment.sequence,
        )
        self.assertEqual(len(preview_data.digest_sequence_views), 2)
        vector_digest_view = preview_data.digest_sequence_views[0]
        insert_digest_view = preview_data.digest_sequence_views[1]
        self.assertEqual(vector_digest_view.role, "Plasmid / vector")
        self.assertEqual(vector_digest_view.used_segments[0].start, selected_vector_fragment.start)
        self.assertEqual(vector_digest_view.used_segments[0].end, selected_vector_fragment.end)
        self.assertEqual(insert_digest_view.used_segments[0].start, selected_insert_fragment.start)
        self.assertTrue(any(marker.enzyme_name == "EcoRI" for marker in vector_digest_view.cut_markers))
        self.assertTrue(any(part.kind == "used" for part in vector_digest_view.preview_parts))
        self.assertTrue(vector_digest_view.double_strand_cut_views)
        self.assertEqual(vector_digest_view.double_strand_cut_views[0].overhang_sequence, "AATT")
        self.assertIn("|", vector_digest_view.double_strand_cut_views[0].top_line)
        self.assertIn("|", vector_digest_view.double_strand_cut_views[0].bottom_line)
        self.assertTrue(vector_digest_view.fragment_options)
        self.assertTrue(any(option.selected for option in vector_digest_view.fragment_options))

    def test_same_enzyme_fragment_selection_supports_circular_vector_backbone(self):
        sequence_files_dir = next(
            (
                parent / "media" / "sequence_files"
                for parent in Path(__file__).resolve().parents
                if (parent / "media" / "sequence_files").exists()
            )
        )
        expected_path = sequence_files_dir / "LvL25_leader_DR_metH_spacer_DR.dna"
        if not expected_path.exists():
            expected_path = sequence_files_dir / "LvL25_leader_DR_metH_sapcer_DR.dna"

        vector_record = next(
            SeqIO.parse(
                sequence_files_dir / "LvL25_leader_sequence_DR_placeholder_DR.dna",
                "snapgene",
            )
        )
        insert_record = next(SeqIO.parse(sequence_files_dir / "spacer_metH.fasta", "fasta"))
        expected_record = next(SeqIO.parse(expected_path, "snapgene"))

        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id=str(vector_record.id),
            name="LvL25 leader DR placeholder DR",
            sequence=str(vector_record.seq).upper(),
            is_circular=True,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id=str(insert_record.id),
            name="spacer_metH",
            sequence=str(insert_record.seq).upper(),
        )
        enzyme = cloning._get_enzyme_by_name("BsmBI")
        vector_backbone_fragment = max(
            cloning._digest_sequence_fragments(
                vector_asset.sequence,
                enzyme,
                is_circular=True,
            ),
            key=lambda fragment: fragment.length,
        )
        insert_fragment = max(
            cloning._digest_sequence_fragments(insert_asset.sequence, enzyme),
            key=lambda fragment: fragment.length,
        )

        preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="BsmBI",
            right_enzyme="BsmBI",
            vector_fragment_index=vector_backbone_fragment.index,
            insert_fragment_index=insert_fragment.index,
        )

        self.assertTrue(preview_data.is_valid)
        self.assertTrue(vector_backbone_fragment.wraps_origin)
        self.assertEqual(preview_data.assembled_sequence, str(expected_record.seq).upper())
        self.assertEqual(preview_data.assembled_length, len(expected_record.seq))
        self.assertTrue(preview_data.digest_sequence_views[0].fragment_options)
        self.assertTrue(
            any(
                option.selected and option.wraps_origin
                for option in preview_data.digest_sequence_views[0].fragment_options
            )
        )

    def test_resolve_cloning_assets_returns_vector_and_insert_assets(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Vector",
            sequence="GAATTC",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGC",
        )

        with mock.patch.object(cloning, "resolve_asset_choice", side_effect=[vector_asset, insert_asset]) as resolve_mock:
            assets = cloning.resolve_cloning_assets(
                user="user-1",
                vector_asset_choice="sequence_file:1:vec1",
                insert_asset_choice="pcr_product:2",
            )

        self.assertEqual(assets.vector_asset, vector_asset)
        self.assertEqual(assets.insert_asset, insert_asset)
        self.assertEqual(resolve_mock.call_count, 2)

    def test_preview_cloning_construct_returns_persistable_payload(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file="vector-file",
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Vector",
            sequence="GAATTCACCCCGGATCC",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product="insert-product",
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )

        with mock.patch.object(
            cloning,
            "_validate_restriction_ligation",
            return_value=cloning.CloningValidationResult(
                is_valid=True,
                validation_messages=("ok",),
                assembled_sequence="assembled-seq",
            ),
        ) as validate_mock:
            preview_data = cloning.preview_cloning_construct(
                vector_asset=vector_asset,
                insert_asset=insert_asset,
                assembly_strategy="restriction_ligation",
                left_enzyme="EcoRI",
                right_enzyme="BamHI",
            )

        self.assertEqual(preview_data.vector_asset.sequence_file, "vector-file")
        self.assertEqual(preview_data.insert_asset.pcr_product, "insert-product")
        self.assertEqual(preview_data.assembled_sequence, "assembled-seq")
        self.assertTrue(preview_data.is_valid)
        self.assertEqual(preview_data.validation_messages, ("ok",))
        validate_mock.assert_called_once_with(
            "GAATTCACCCCGGATCC",
            "ATGCATGC",
            "EcoRI",
            "BamHI",
            vector_fragment_index=None,
            insert_fragment_index=None,
            vector_is_circular=False,
            insert_is_circular=False,
        )

    def test_preview_cloning_construct_defaults_result_topology_from_vector(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Circular vector",
            sequence="GAATTCACCCCGGATCC",
            is_circular=True,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )

        preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
        )
        linear_preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
            is_circular=False,
        )

        self.assertTrue(preview_data.is_valid)
        self.assertTrue(preview_data.is_circular)
        self.assertFalse(linear_preview_data.is_circular)

    def test_assembly_visual_preview_supports_multiple_selected_enzymes_on_linear_sequence(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Linear vector",
            sequence="AAAAGAATTCTTTTGGATCCAAAA",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )

        visual_preview = cloning.build_cloning_assembly_visual_preview(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            map_enzyme_names=("EcoRI", "BamHI"),
        )

        enzyme_names = {site.enzyme_name for site in visual_preview.vector_map.restriction_sites}
        self.assertEqual(visual_preview.selected_enzyme_names, ("EcoRI", "BamHI"))
        self.assertEqual(enzyme_names, {"EcoRI", "BamHI"})
        self.assertTrue(
            all(site.stable_id.startswith("vector-") for site in visual_preview.vector_map.restriction_sites)
        )
        self.assertTrue(
            all(site.position >= 0 for site in visual_preview.vector_map.restriction_sites)
        )
        self.assertGreaterEqual(len(visual_preview.vector_map.digest_fragments), 3)
        self.assertEqual(visual_preview.vector_map.digest_fragments[0].start, 0)
        self.assertEqual(
            visual_preview.vector_map.digest_fragments[-1].end,
            len(vector_asset.sequence),
        )

    def test_assembly_visual_preview_renders_full_maps_before_enzyme_selection(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Circular vector",
            sequence="AAAAGAATTCTTTT",
            is_circular=True,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Linear insert",
            sequence="ATGCATGC",
        )

        visual_preview = cloning.build_cloning_assembly_visual_preview(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
        )

        self.assertTrue(visual_preview.vector_map.is_circular)
        self.assertEqual(visual_preview.vector_map.map_shape, "circular")
        self.assertTrue(visual_preview.vector_map.is_circular_map)
        self.assertEqual(visual_preview.vector_map.map_shape_label, "circular map")
        self.assertEqual(visual_preview.vector_map.source_topology_label, "circular source")
        self.assertFalse(visual_preview.insert_map.is_circular)
        self.assertEqual(visual_preview.insert_map.map_shape, "linear")
        self.assertFalse(visual_preview.insert_map.is_circular_map)
        self.assertEqual(visual_preview.insert_map.map_shape_label, "linear map")
        self.assertEqual(visual_preview.insert_map.source_topology_label, "linear source")
        self.assertEqual(visual_preview.selected_enzyme_names, ())
        self.assertEqual(visual_preview.vector_map.restriction_sites, ())
        self.assertEqual(visual_preview.insert_map.restriction_sites, ())
        self.assertEqual(visual_preview.vector_map.regions[0].length, len(vector_asset.sequence))
        self.assertEqual(visual_preview.insert_map.regions[0].length, len(insert_asset.sequence))

    def test_assembly_visual_preview_uses_topology_based_map_shapes(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Linear vector source",
            sequence="AAAAGAATTCTTTT",
            is_circular=False,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Circular insert source",
            sequence="ATGCATGC",
            is_circular=True,
        )

        visual_preview = cloning.build_cloning_assembly_visual_preview(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
        )

        self.assertFalse(visual_preview.vector_map.is_circular)
        self.assertEqual(visual_preview.vector_map.map_shape, "linear")
        self.assertFalse(visual_preview.vector_map.is_circular_map)
        self.assertTrue(visual_preview.insert_map.is_circular)
        self.assertEqual(visual_preview.insert_map.map_shape, "circular")
        self.assertTrue(visual_preview.insert_map.is_circular_map)

    def test_assembly_visual_preview_supports_circular_digest_wraparound(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Circular vector",
            sequence="AAAAGAATTCTTTT",
            is_circular=True,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )

        visual_preview = cloning.build_cloning_assembly_visual_preview(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            map_enzyme_names=("EcoRI",),
        )

        self.assertTrue(visual_preview.vector_map.is_circular)
        self.assertEqual(len(visual_preview.vector_map.digest_fragments), 1)
        self.assertTrue(visual_preview.vector_map.digest_fragments[0].wraps_origin)
        self.assertEqual(
            visual_preview.vector_map.digest_fragments[0].length,
            len(vector_asset.sequence),
        )

    def test_assembly_visual_preview_handles_selected_enzyme_with_no_cut_sites(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="No-cut vector",
            sequence="AAAATTTTCCCC",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="No-cut insert",
            sequence="TTTTCCCCAAAA",
        )

        visual_preview = cloning.build_cloning_assembly_visual_preview(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            map_enzyme_names=("EcoRI",),
        )

        self.assertEqual(visual_preview.vector_map.restriction_sites, ())
        self.assertEqual(visual_preview.insert_map.restriction_sites, ())
        self.assertEqual(visual_preview.vector_map.digest_fragments, ())
        self.assertEqual(visual_preview.enzyme_summaries[0].vector_cut_count, 0)
        self.assertEqual(visual_preview.enzyme_summaries[0].insert_cut_count, 0)

    def test_build_cloning_construct_record_uses_saved_result_topology(self):
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Vector",
            sequence="GAATTCACCCCGGATCC",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )
        construct = SimpleNamespace(
            name="Circular Construct",
            description="",
            vector_source_type="pcr_product",
            vector_sequence_file=None,
            vector_template_name="",
            vector_pcr_product=None,
            vector_record_id="vec1",
            vector_fragment_index=None,
            insert_source_type="pcr_product",
            insert_sequence_file=None,
            insert_template_name="",
            insert_pcr_product=None,
            insert_record_id="ins1",
            insert_fragment_index=None,
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
            assembled_sequence="GATGCATGCGATCC",
            vector_asset_label="Vector",
            insert_asset_label="Insert",
            is_circular=True,
        )

        with mock.patch.object(
            cloning_exports,
            "_resolve_construct_asset",
            side_effect=[vector_asset, insert_asset],
        ):
            record = cloning_exports.build_cloning_construct_record(construct)

        self.assertEqual(record.annotations["topology"], "circular")

    def test_export_snapgene_fragment_construct_preserves_features_and_round_trips(self):
        sequence_files_dir = next(
            (
                parent / "media" / "sequence_files"
                for parent in Path(__file__).resolve().parents
                if (parent / "media" / "sequence_files").exists()
            )
        )
        vector_record = next(
            SeqIO.parse(
                sequence_files_dir / "LvL25_leader_sequence_DR_placeholder_DR.dna",
                "snapgene",
            )
        )
        insert_record = next(SeqIO.parse(sequence_files_dir / "spacer_metH.fasta", "fasta"))
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id=str(vector_record.id),
            name="LvL25 leader DR placeholder DR",
            sequence=str(vector_record.seq).upper(),
            is_circular=True,
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file=None,
            pcr_product=None,
            template_name=None,
            record_id=str(insert_record.id),
            name="spacer_metH",
            sequence=str(insert_record.seq).upper(),
        )
        enzyme = cloning._get_enzyme_by_name("BsmBI")
        vector_fragment = max(
            cloning._digest_sequence_fragments(
                vector_asset.sequence,
                enzyme,
                is_circular=True,
            ),
            key=lambda fragment: fragment.length,
        )
        insert_fragment = max(
            cloning._digest_sequence_fragments(insert_asset.sequence, enzyme),
            key=lambda fragment: fragment.length,
        )
        preview_data = cloning.preview_cloning_construct(
            vector_asset=vector_asset,
            insert_asset=insert_asset,
            assembly_strategy="restriction_ligation",
            left_enzyme="BsmBI",
            right_enzyme="BsmBI",
            vector_fragment_index=vector_fragment.index,
            insert_fragment_index=insert_fragment.index,
        )
        construct = SimpleNamespace(
            name="LvL25 leader DR metH spacer DR",
            description="Feature-rich circular construct",
            vector_source_type="sequence_file",
            vector_sequence_file=None,
            vector_template_name="",
            vector_pcr_product=None,
            vector_record_id=str(vector_record.id),
            vector_fragment_index=vector_fragment.index,
            insert_source_type="sequence_file",
            insert_sequence_file=None,
            insert_template_name="",
            insert_pcr_product=None,
            insert_record_id=str(insert_record.id),
            insert_fragment_index=insert_fragment.index,
            left_enzyme="BsmBI",
            right_enzyme="BsmBI",
            assembled_sequence=preview_data.assembled_sequence,
            vector_asset_label="LvL25 leader DR placeholder DR",
            insert_asset_label="spacer_metH",
            is_circular=True,
        )

        with mock.patch.object(
            cloning_exports,
            "_resolve_construct_asset",
            side_effect=[vector_asset, insert_asset],
        ), mock.patch.object(
            cloning_exports,
            "_build_asset_bundle",
            side_effect=[
                {
                    "sequence": vector_asset.sequence,
                    "features": list(vector_record.features),
                    "annotations": dict(getattr(vector_record, "annotations", {}) or {}),
                },
                {
                    "sequence": insert_asset.sequence,
                    "features": list(getattr(insert_record, "features", [])),
                    "annotations": dict(getattr(insert_record, "annotations", {}) or {}),
                },
            ],
        ):
            content = cloning_exports.export_cloning_construct_sequence(
                construct,
                "genbank",
            )
            records = cloning_exports._validate_exported_sequence_content(
                content=content,
                file_type="genbank",
            )

        labels = {
            (feature.qualifiers.get("label") or [""])[0]
            for feature in records[0].features
        }
        self.assertEqual(records[0].annotations["topology"], "circular")
        self.assertGreaterEqual(len(records[0].features), 30)
        self.assertIn("leader_sequence_CRISPR", labels)
        self.assertIn("KanR", labels)

    def test_validate_exported_sequence_content_rejects_unparseable_generated_assembly(self):
        with mock.patch.object(cloning_exports.SeqIO, "parse", side_effect=ValueError("bad export")):
            with self.assertRaises(ValueError) as exc:
                cloning_exports._validate_exported_sequence_content(
                    content="LOCUS       bad\n",
                    file_type="genbank",
                )

        self.assertIn("exported sequence file could not be parsed", str(exc.exception))

    def test_validate_exported_sequence_content_accepts_generated_fasta(self):
        records = cloning_exports._validate_exported_sequence_content(
            content=">assembled\nATGCATGC\n",
            file_type="fasta",
        )

        self.assertEqual(len(records), 1)
        self.assertEqual(str(records[0].seq), "ATGCATGC")

    def test_build_cloning_construct_detail_display_returns_read_model(self):
        construct = SimpleNamespace(
            vector_source_type="sequence_file",
            vector_sequence_file="vector-file",
            vector_pcr_product=None,
            vector_record_id="vec1",
            insert_source_type="pcr_product",
            insert_sequence_file=None,
            insert_pcr_product="insert-product",
            insert_record_id="ins1",
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
            assembly_strategy="restriction_ligation",
            assembled_sequence="GATGCATGCGATCC",
        )
        vector_asset = cloning.ResolvedCloningAsset(
            source_type="sequence_file",
            sequence_file="vector-file",
            pcr_product=None,
            template_name=None,
            record_id="vec1",
            name="Vector",
            sequence="GAATTCACCCCGGATCC",
        )
        insert_asset = cloning.ResolvedCloningAsset(
            source_type="pcr_product",
            sequence_file=None,
            pcr_product="insert-product",
            template_name=None,
            record_id="ins1",
            name="Insert",
            sequence="ATGCATGC",
        )

        with mock.patch.object(cloning, "_resolve_construct_asset", side_effect=[vector_asset, insert_asset]), \
             mock.patch.object(
                 cloning,
                 "_build_cut_site_preview",
                 side_effect=[
                     cloning.CloningCutSitePreview(
                         enzyme_name="EcoRI",
                         site_sequence="GAATTC",
                         vector_cut_positions=(1,),
                         insert_recognition_site_positions=(),
                     ),
                     cloning.CloningCutSitePreview(
                         enzyme_name="BamHI",
                         site_sequence="GGATCC",
                         vector_cut_positions=(11,),
                         insert_recognition_site_positions=(),
                     ),
                 ],
             ), \
             mock.patch.object(cloning, "_get_enzyme_by_name", side_effect=["eco", "bam"]), \
             mock.patch.object(cloning, "_find_cut_positions", side_effect=[[1], [9]]):
            detail_display = cloning.build_cloning_construct_detail_display(construct)

        self.assertEqual(detail_display.junction_context_window, 12)
        self.assertEqual(len(detail_display.cut_site_previews), 2)
        self.assertEqual(detail_display.cut_site_previews[0].enzyme_name, "EcoRI")
        self.assertEqual(detail_display.cut_site_previews[1].enzyme_name, "BamHI")
        self.assertEqual(detail_display.junction_contexts[0].display, "G|ATGCATGCGATC")
        self.assertEqual(detail_display.junction_contexts[1].display, "GATGCATGC|GATCC")
        self.assertEqual(detail_display.source_errors, ())

    def test_build_cloning_construct_detail_display_reports_source_errors(self):
        construct = SimpleNamespace(
            vector_source_type="sequence_file",
            vector_sequence_file="vector-file",
            vector_pcr_product=None,
            vector_record_id="vec1",
            insert_source_type="pcr_product",
            insert_sequence_file=None,
            insert_pcr_product="insert-product",
            insert_record_id="ins1",
            left_enzyme="EcoRI",
            right_enzyme="BamHI",
            assembly_strategy="restriction_ligation",
            assembled_sequence="GATGCATGCGATCC",
        )

        with mock.patch.object(
            cloning,
            "_resolve_construct_asset",
            side_effect=[ValueError("Vector source is unavailable."), cloning.ResolvedCloningAsset(
                source_type="pcr_product",
                sequence_file=None,
                pcr_product="insert-product",
                template_name=None,
                record_id="ins1",
                name="Insert",
                sequence="ATGCATGC",
            )],
        ):
            detail_display = cloning.build_cloning_construct_detail_display(construct)

        self.assertEqual(detail_display.cut_site_previews, ())
        self.assertEqual(detail_display.junction_contexts, ())
        self.assertEqual(detail_display.source_errors, ("Vector source is unavailable.",))

    def test_build_cloning_construct_detail_display_uses_persisted_snapshot(self):
        construct = SimpleNamespace(
            detail_display_snapshot={
                "junction_context_window": 12,
                "cut_site_previews": [
                    {
                        "enzyme_name": "EcoRI",
                        "site_sequence": "GAATTC",
                        "vector_cut_positions": [2],
                        "insert_recognition_site_positions": [],
                    }
                ],
                "junction_contexts": [
                    {
                        "label": "Vector -> insert junction",
                        "left_context": "G",
                        "right_context": "ATGCATGCGATC",
                    }
                ],
                "source_errors": [],
            },
        )

        detail_display = cloning.build_cloning_construct_detail_display(construct)

        self.assertEqual(detail_display.junction_context_window, 12)
        self.assertEqual(len(detail_display.cut_site_previews), 1)
        self.assertEqual(detail_display.cut_site_previews[0].enzyme_name, "EcoRI")
        self.assertEqual(detail_display.junction_contexts[0].display, "G|ATGCATGCGATC")
        self.assertEqual(detail_display.source_errors, ())

    def test_build_cloning_construct_detail_display_reads_legacy_insert_site_snapshot_key(self):
        construct = SimpleNamespace(
            detail_display_snapshot={
                "junction_context_window": 12,
                "cut_site_previews": [
                    {
                        "enzyme_name": "EcoRI",
                        "site_sequence": "GAATTC",
                        "vector_cut_positions": [2],
                        "insert_site_positions": [5],
                    }
                ],
                "junction_contexts": [],
                "source_errors": [],
            },
        )

        detail_display = cloning.build_cloning_construct_detail_display(construct)

        self.assertEqual(
            detail_display.cut_site_previews[0].insert_recognition_site_positions,
            (5,),
        )
