import itertools
import json
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from inverted_index import InvertedIndex, read_sequences, select_seed_positions
from integration_example import search_one


def oracle_distance(a, b):
    """Small full-matrix test oracle; not the group's production DP verifier."""
    matrix = [[0] * (len(b) + 1) for _ in range(len(a) + 1)]
    for i in range(len(a) + 1):
        matrix[i][0] = i
    for j in range(len(b) + 1):
        matrix[0][j] = j
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            matrix[i][j] = min(matrix[i-1][j] + 1, matrix[i][j-1] + 1,
                               matrix[i-1][j-1] + (a[i-1] != b[j-1]))
    return matrix[-1][-1]


class IndexTests(unittest.TestCase):
    def test_positions_and_disjointness(self):
        self.assertEqual(select_seed_positions(12, 3, 2), (0, 4, 9))
        self.assertEqual(select_seed_positions(12, 3, 2, "prefix"), (0, 3, 6))
        self.assertEqual(select_seed_positions(8, 3, 2), ())
        for m in range(41):
            for k in range(1, 7):
                for d in range(7):
                    for strategy in ("spread", "prefix"):
                        starts = select_seed_positions(m, k, d, strategy)
                        if m < k * (d + 1):
                            self.assertEqual(starts, ())
                        else:
                            self.assertEqual(len(starts), d + 1)
                            self.assertGreaterEqual(starts[0], 0)
                            self.assertLessEqual(starts[-1] + k, m)
                            self.assertTrue(all(b - a >= k for a, b in zip(starts, starts[1:])))

    def test_all_reference_windows_and_distinct_ids(self):
        index = InvertedIndex({"r1": "AAAA", "r2": "AAAA", "r3": "CASSLG", "short": "A"}, 2)
        self.assertEqual(index.postings["AA"], {"r1", "r2"})
        for seed in ("CA", "AS", "SS", "SL", "LG"):
            self.assertIn("r3", index.postings[seed])
        self.assertEqual(index.generate_candidates("A", 0).candidate_ids, {"short"})

    def test_snapshot_cannot_be_accidentally_modified(self):
        refs = {"a": "CASS"}
        index = InvertedIndex(refs, 2)
        refs["a"] = "YYYY"
        self.assertEqual(index.references["a"], "CASS")
        with self.assertRaises(TypeError):
            index.references["a"] = "YYYY"

    def test_tiny_exhaustive_completeness_and_candidate_definition(self):
        # All 31 strings of length 0..4 over {A,C}, including an empty string.
        strings = ["".join(p) for n in range(5) for p in itertools.product("AC", repeat=n)]
        refs = {str(i): s for i, s in enumerate(strings)}
        truth_distances = {q: {rid: oracle_distance(q, s) for rid, s in refs.items()} for q in strings}
        for k in range(1, 6):
            index = InvertedIndex(refs, k)
            for query in strings:
                for d in (0, 1, 2, 3, 5):
                    for strategy in ("spread", "prefix"):
                        result = index.generate_candidates(query, d, strategy)
                        truth = {rid for rid, distance in truth_distances[query].items() if distance <= d}
                        self.assertLessEqual(truth, result.candidate_ids)
                        # Independent string-membership scan also checks extra candidates.
                        expected = {rid for rid, s in refs.items()
                                    if abs(len(s)-len(query)) <= d
                                    and (result.mode == "fallback" or any(seed in s for seed in result.seeds))}
                        self.assertEqual(result.candidate_ids, expected)

    def test_indel_shift_and_extra_candidate(self):
        refs = {"R05": "CASSLG", "R06": "CDSSLG", "R07": "CGASSLG",
                "R08": "ASSLG", "R09": "CASWWW", "R10": "CASSLGGG"}
        result = InvertedIndex(refs, 2).generate_candidates("CASSLG", 1)
        self.assertEqual(result.seeds, ("CA", "LG"))
        self.assertEqual(result.candidate_ids, {"R05", "R06", "R07", "R08", "R09"})
        self.assertEqual(result.union_count, 6)

    def test_empty_safe_union_is_not_fallback(self):
        result = InvertedIndex({"a": "CASSLG"}, 2).generate_candidates("YYYYYY", 1)
        self.assertEqual(result.mode, "indexed")
        self.assertEqual(result.candidate_ids, frozenset())

    def test_large_threshold_and_empty_database(self):
        result = InvertedIndex({"a": "C", "b": "CASSLG"}, 3).generate_candidates("C", 1000000)
        self.assertEqual(result.mode, "fallback")
        self.assertEqual(result.candidate_ids, {"a", "b"})
        self.assertFalse(InvertedIndex({}, 2).generate_candidates("CA", 0).candidate_ids)

    def test_invalid_parameters(self):
        for k in (0, -1, True, 2.5, "3"):
            with self.assertRaises(ValueError):
                InvertedIndex({}, k)
        index = InvertedIndex({}, 2)
        for d in (-1, True, 0.5, "2"):
            with self.assertRaises(ValueError):
                index.generate_candidates("CA", d)
        for q in ("ca", "CX", " CA", None):
            with self.assertRaises(ValueError):
                index.generate_candidates(q, 1)
        with self.assertRaises(ValueError):
            index.generate_candidates("CA", 1, "rare")

    def test_tsv_and_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            refs = root / "references.tsv"
            queries = root / "queries.tsv"
            refs.write_text("reference_id\tsequence\nr1\tCASSLG\nr2\tCDSSLG\n", encoding="utf-8")
            queries.write_text("query_id\tsequence\tcategory\nq1\tCASSLG\tignored\nq2\tYYYYYY\tignored\n", encoding="utf-8")
            self.assertEqual(read_sequences(refs, "reference_id")["r1"], "CASSLG")
            command = [sys.executable, str(Path(__file__).with_name("inverted_index.py")),
                       "--reference", str(refs), "--queries", str(queries), "--k", "2", "--d", "1"]
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
            results = [json.loads(line) for line in completed.stdout.splitlines()]
            self.assertEqual(results[0]["candidate_ids"], ["r1", "r2"])
            self.assertEqual(results[1]["candidate_ids"], [])
            self.assertEqual(json.loads(completed.stderr)["query_count"], 2)
            failed = subprocess.run(command + ["--output", str(refs)], capture_output=True, text=True)
            self.assertEqual(failed.returncode, 2)
            for text in ("wrong\tsequence\na\tCA\n", "reference_id\tsequence\na\tCA\na\tAC\n",
                         "reference_id\tsequence\na\t\n", "reference_id\tsequence\na\tCA\textra\n"):
                refs.write_text(text, encoding="utf-8")
                with self.assertRaises(ValueError):
                    read_sequences(refs, "reference_id")

    def test_integration_contract(self):
        index = InvertedIndex({"exact": "CASSLG", "one": "CDSSLG", "extra": "CASWWW"}, 2)
        verifier = lambda q, s, d: oracle_distance(q, s)
        output = search_one(index, "q1", "CASSLG", 1, verifier, annotations={"exact": [{"source_line": 1, "fields": {}}]})
        self.assertEqual([m["reference_id"] for m in output["matches"]], ["exact", "one"])
        self.assertEqual([m["distance"] for m in output["matches"]], [0, 1])
        self.assertEqual(output["matches"][1]["annotations"], [])
        self.assertEqual(search_one(index, "q2", "YYYYYY", 1, verifier), {"query_id": "q2", "matches": []})
        with self.assertRaises(ValueError):
            search_one(index, "q", "CASSLG", 1, lambda q, s, d: True)
        self.assertEqual(search_one(index, "q", "CASSLG", 0, lambda q, s, d: None)["matches"], [])

    def test_length_groups_are_read_only(self):
        index = InvertedIndex({"r": "CA"}, k=2)
        with self.assertRaises(TypeError):
            index.by_length[2] = set()
        with self.assertRaises(AttributeError):
            index.by_length[2].clear()
        self.assertEqual(index.generate_candidates("C", d=1).candidate_ids, {"r"})

    def test_quoted_tsv_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "queries.tsv"
            header = "query_id\tsequence\tnote\n"
            path.write_text(header + 'q1\tCA\t"line one\nline two"\nq2\tAC\tok\n', encoding="utf-8")
            self.assertEqual(read_sequences(path, "query_id"), {"q1": "CA", "q2": "AC"})
            path.write_text(header + 'q1\tCA\t"unclosed\nq2\tAC\tok\n', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, r"queries\.tsv:3: malformed TSV"):
                read_sequences(path, "query_id")

    def test_file_output_and_input_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            refs, queries, output = (root / name for name in ("ref.tsv", "query.tsv", "out.jsonl"))
            refs.write_text("reference_id\tsequence\nr\tCA\n", encoding="utf-8")
            queries.write_text("query_id\tsequence\nq\tCA\n", encoding="utf-8")
            command = [sys.executable, str(Path(__file__).with_name("inverted_index.py")),
                       "--reference", str(refs), "--queries", str(queries), "--d", "0",
                       "--output", str(output)]
            output.write_text("previous result\n", encoding="utf-8")
            subprocess.run(command, check=True, capture_output=True)
            self.assertEqual(json.loads(output.read_text(encoding="utf-8"))["candidate_ids"], ["r"])
            for source in (refs, queries):
                with self.subTest(source=source.name):
                    original = source.read_bytes()
                    output.unlink()
                    try:
                        os.link(source, output)
                    except OSError as error:
                        self.skipTest(f"Hard links unavailable: {error}")
                    run = subprocess.run(command, capture_output=True, text=True)
                    self.assertEqual(run.returncode, 2)
                    self.assertIn("must not overwrite an input file", run.stderr)
                    self.assertEqual(source.read_bytes(), original)

    def test_verified_matches_for_both_seed_strategies(self):
        refs = {"exact": "CASSLG", "insert": "CGASSLG", "delete": "ASSLG", "extra": "CASWWW"}
        index = InvertedIndex(refs, k=2)
        annotations = {"insert": [{"source_line": 7, "fields": {"source": "test"}}]}
        for query, d in (("CASSLG", 1), ("CAS", 1), ("CASS", 1), ("YYYYYY", 0)):
            expected = []
            for rid, sequence in sorted(refs.items()):
                distance = oracle_distance(query, sequence)
                if distance <= d:
                    expected.append({"reference_id": rid, "sequence": sequence, "distance": distance,
                                     "annotations": annotations.get(rid, [])})
            for strategy in ("spread", "prefix"):
                with self.subTest(query=query, d=d, strategy=strategy):
                    result = search_one(index, "q", query, d, lambda a, b, limit: oracle_distance(a, b),
                                        strategy=strategy, annotations=annotations)
                    self.assertEqual(result, {"query_id": "q", "matches": expected})


if __name__ == "__main__":
    unittest.main()
