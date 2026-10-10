import csv
import gzip
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

import plasmidprobe as checker


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


class SequenceInputTests(unittest.TestCase):
    def test_reads_synthetic_fasta_and_fastq(self):
        name, sequence = checker.read_single_fasta(str(FIXTURES / "demo_reference.fasta"))
        self.assertEqual(name, "synthetic_candidate")
        self.assertEqual(len(sequence), 320)

        reads = list(checker.read_fastq(str(FIXTURES / "demo_reads.fastq")))
        self.assertEqual(len(reads), 3)
        self.assertEqual(reads[1][0], "demo_feature_2_reverse")
        self.assertTrue(all(len(read) == 80 for _, read in reads))

    def test_reads_gzipped_fastq(self):
        fixture = (FIXTURES / "demo_reads.fastq").read_text()
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "reads.fastq.gz"
            with gzip.open(path, "wt") as handle:
                handle.write(fixture)
            reads = list(checker.read_fastq(str(path)))

        self.assertEqual(len(reads), 3)

    def test_cleans_sequences_and_reverse_complements(self):
        self.assertEqual(checker.clean_seq("acgU-Nx"), "ACGN")
        self.assertEqual(checker.reverse_complement("AACGTN"), "NACGTT")


class FeatureInputTests(unittest.TestCase):
    def test_reads_csv_coordinates_as_zero_based_intervals(self):
        name, sequence = checker.read_single_fasta(str(FIXTURES / "demo_reference.fasta"))
        features = checker.read_features_csv(str(FIXTURES / "demo_features.csv"), len(sequence))

        self.assertEqual(name, "synthetic_candidate")
        self.assertEqual(features[0], checker.Feature("demo_feature_1", 0, 80, "+", "demo_feature_1"))
        self.assertEqual(features[-1].start, 240)
        self.assertEqual(features[-1].end, 320)

    def test_reads_text_index_with_joined_intervals(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "features.index.txt"
            path.write_text(">joined_feature\nstart {1, 11}, end {4, 14}\n")
            features = checker.read_features(str(path), reference_length=20)

        self.assertEqual([feature.name for feature in features], ["joined_feature_part1", "joined_feature_part2"])
        self.assertEqual([(feature.start, feature.end) for feature in features], [(0, 4), (10, 14)])

    def test_reads_genbank_annotations(self):
        record = SeqRecord(Seq("ACGT" * 10), id="demo", name="demo", description="synthetic record")
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(
                FeatureLocation(4, 16, strand=-1),
                type="CDS",
                qualifiers={"gene": ["synthetic_gene"]},
            )
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "demo.gb"
            SeqIO.write(record, path, "genbank")
            features = checker.read_features(str(path), reference_length=len(record.seq))

        self.assertEqual(len(features), 1)
        self.assertEqual(features[0].name, "synthetic_gene")
        self.assertEqual((features[0].start, features[0].end, features[0].strand), (4, 16, "-"))


class AnalysisTests(unittest.TestCase):
    def setUp(self):
        self.reference_name, self.reference = checker.read_single_fasta(
            str(FIXTURES / "demo_reference.fasta")
        )
        self.features = checker.read_features_csv(
            str(FIXTURES / "demo_features.csv"), len(self.reference)
        )

    def analyse_demo(self):
        return checker.analyse_features(
            reads_path=str(FIXTURES / "demo_reads.fastq"),
            ref_name=self.reference_name,
            ref_seq=self.reference,
            features=self.features,
            aligner="smithwaterman",
            threads=1,
            match_score=2,
            mismatch_penalty=3,
            gap_penalty=5,
            min_alignment_score=40,
            min_alignment_score_density=0.0,
            min_read_identity=0.75,
            min_feature_breadth=0.80,
            min_feature_depth=1.0,
            min_evidence_reads=1,
        )

    def test_scores_present_and_uncovered_features(self):
        results = self.analyse_demo()
        status_by_name = {result.feature: result.status for result in results}

        self.assertEqual(status_by_name["demo_feature_1"], "EVIDENCE_DETECTED")
        self.assertEqual(status_by_name["demo_feature_2"], "EVIDENCE_DETECTED")
        self.assertEqual(status_by_name["demo_feature_3_unobserved"], "NOT_DETECTED")
        self.assertEqual(status_by_name["demo_feature_4"], "EVIDENCE_DETECTED")
        self.assertEqual(results[1].evidence_reads, 1)
        self.assertGreater(results[1].mean_identity, 0.95)

    def test_probe_mode_finds_observed_features(self):
        probes = checker.generate_probes(self.reference, self.features, probe_length=50, probe_step=25)
        hits = checker.probe_hits_from_reads(str(FIXTURES / "demo_reads.fastq"), probes)

        self.assertIn("demo_feature_1", hits)
        self.assertIn("demo_feature_2", hits)
        self.assertNotIn("demo_feature_3_unobserved", hits)

    def test_writes_results_and_creates_parent_directory(self):
        results = self.analyse_demo()
        with tempfile.TemporaryDirectory() as temp_dir:
            output = Path(temp_dir) / "nested" / "results.csv"
            checker.write_results(results, str(output))
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))

        self.assertEqual(len(rows), 4)
        self.assertEqual(rows[0]["status"], "EVIDENCE_DETECTED")
        self.assertEqual(rows[0]["evidence_reads"], "1")
        self.assertEqual(rows[2]["status"], "NOT_DETECTED")


if __name__ == "__main__":
    unittest.main()
