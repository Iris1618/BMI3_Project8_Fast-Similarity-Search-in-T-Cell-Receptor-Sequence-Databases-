"""Small portable tests; the real VDJdb package is not required."""

import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Inverted_Index.inverted_index import InvertedIndex

from .baseline import exhaustive_matches, search_one
from .compare import benchmark, load_truth, match_differences, memory_comparison, validation_rows


class ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.references = {"exact": "CASS", "sub": "CATS", "ins": "CASSS",
                           "del": "CAS", "none": "YYYY"}
        self.queries = {"q1": "CASS", "q2": "WWWW", "short": "C"}
        self.index = InvertedIndex(self.references, k=2)

    def test_baseline_boundaries_and_no_match(self):
        self.assertEqual(exhaustive_matches("CASS", self.references, 0), {"exact": 0})
        self.assertEqual(exhaustive_matches("CASS", self.references, 1),
                         {"exact": 0, "sub": 1, "ins": 1, "del": 1})
        self.assertEqual(search_one("q", "WWWW", self.references, 0),
                         {"query_id": "q", "matches": []})

    def test_validation_both_strategies_and_fallback(self):
        for strategy in ("spread", "prefix"):
            for d in (0, 1, 2, 3):
                with self.subTest(strategy=strategy, d=d):
                    rows = validation_rows(self.references, self.queries, self.index, d, strategy)
                    self.assertEqual(len(rows), 6)
                    self.assertTrue(all(row["match"] for row in rows))
        rows = validation_rows(self.references, {"short": "C"}, self.index, 1, "spread")
        self.assertTrue(all(row["index_mode"] == "fallback" for row in rows))

    def test_wrong_distance_is_visible(self):
        missing, unexpected, wrong = match_differences({"a": 0, "b": 1}, {"a": 1, "c": 1})
        self.assertEqual(missing, ["b"])
        self.assertEqual(unexpected, ["c"])
        self.assertEqual(wrong, [("a", 0, 1)])

    def test_benchmark_detects_a_broken_verifier(self):
        with patch("Exhaustive_Benchmark.compare.verify", return_value=None):
            with self.assertRaisesRegex(ValueError, "Correctness check failed"):
                benchmark(self.references, self.queries, 1, 2, "spread", 2)

    def test_benchmark_includes_separate_costs_and_samples(self):
        metadata, rows = benchmark(self.references, self.queries, 1, 2, "spread", 2)
        self.assertEqual(metadata["query_ids"], list(self.queries))
        self.assertEqual(len(metadata["index_build_ms_samples"]), 2)
        self.assertEqual({row["method"] for row in rows}, {
            "exhaustive_basic", "indexed_basic", "indexed_optimized",
            "candidate_basic", "candidate_optimized"})
        self.assertTrue(all(len(json.loads(row["samples"])) == 2 for row in rows))

    def test_zero_candidates_has_no_per_pair_timing(self):
        metadata, rows = benchmark({"r": "AAAA"}, {"q": "YYYY"}, 0, 2, "spread", 2)
        self.assertEqual(metadata["candidate_count_total"], 0)
        self.assertTrue(all(row["mean"] is None for row in rows if row["method"].startswith("candidate_")))

    def test_memory_results_have_explicit_scope(self):
        rows = memory_comparison(self.references, self.queries, 1, 2, "spread")
        self.assertEqual({row["method"] for row in rows},
                         {"exhaustive_basic", "indexed_basic", "indexed_optimized"})
        self.assertTrue(all(row["incremental_peak_python_bytes"] > 0 for row in rows))
        self.assertTrue(all("not RSS" in row["scope"] for row in rows))

    def test_truth_summary_includes_zero_matches_and_rejects_missing_rows(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            summary = root / "summary.tsv"
            truth = root / "truth.tsv"
            summary.write_text("reference_set\tthreshold\tquery_id\tn_matches\tn_reference\n"
                               "tiny\t0\tq1\t1\t5\n"
                               "tiny\t0\tq2\t0\t5\n", encoding="utf-8")
            truth.write_text("reference_set\tthreshold\tquery_id\treference_id\tdistance\n"
                             "tiny\t0\tq1\texact\t0\n", encoding="utf-8")
            loaded = load_truth(truth, summary, "tiny", {"q1": "CASS", "q2": "WWWW"}, 0, 5)
            self.assertEqual(loaded["q1"], {"exact": 0})
            self.assertEqual(loaded["q2"], {})
            with self.assertRaisesRegex(ValueError, "does not cover"):
                load_truth(truth, summary, "tiny", {"missing": "A"}, 0, 5)
            with self.assertRaisesRegex(ValueError, "reference count"):
                load_truth(truth, summary, "tiny", {"q1": "CASS"}, 0, 99)


if __name__ == "__main__":
    unittest.main()
