import tempfile
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
        self.assertTrue(obj.saved)
        self.assertEqual(obj.users.added, user)