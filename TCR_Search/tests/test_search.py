"""Exercise the CLI's input selection and complete search output."""

import contextlib
import gzip
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import main


class SearchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.reference = self.root / "reference.tsv"
        self.reference.write_text(
            "reference_id\tsequence\nR1\tCASSF\nR2\tCASAF\n", encoding="utf-8"
        )
        self.queries = self.root / "queries.tsv"
        self.queries.write_text(
            "query_id\tsequence\tcategory\tsource_id\n"
            "Q1\tCASSF\texact_copy\tR2\nQ2\tYYYYF\tcustom\t\n", encoding="utf-8"
        )
        self.records = [
            {"source_line": 12, "fields": {"cdr3": "CASSF", "mhc.a": "A", "note": ""}},
            {"source_line": 13, "fields": {"cdr3": "CASSF", "mhc.a": "B", "note": "中文"}},
        ]
        self.annotation_rows = [
            {"reference_id": "R1", "records": self.records},
            {"reference_id": "R2", "records": [
                {"source_line": 14, "fields": {"cdr3": "CASAF"}}
            ]},
        ]

    def args(self, command="search"):
        return [command, "--reference", str(self.reference), "--queries", str(self.queries)]

    def invoke(self, args):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
            result = main.main(args)
        return result, stdout.getvalue(), stderr.getvalue()

    def write_annotations(self, path):
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "wt", encoding="utf-8") as file:
            for row in self.annotation_rows:
                file.write(json.dumps(row, ensure_ascii=False) + "\n")

    def test_both_methods_return_complete_records_and_unmatched_queries(self):
        self.write_annotations(self.root / "annotations.jsonl.gz")
        results = []
        for method in ("indexed", "exhaustive"):
            code, stdout, _ = self.invoke(self.args() + ["--method", method, "--d", "1", "--k", "2"])
            self.assertEqual(code, 0)
            rows = [json.loads(line) for line in stdout.splitlines()]
            self.assertEqual(rows[0]["query_id"], "Q1")
            self.assertEqual(rows[0]["query_sequence"], "CASSF")
            self.assertEqual(rows[0]["query_metadata"], {"category": "exact_copy", "source_id": "R2"})
            self.assertEqual([(m["reference_id"], m["distance"]) for m in rows[0]["matches"]],
                             [("R1", 0), ("R2", 1)])
            self.assertEqual(rows[0]["matches"][0]["annotations"], self.records)
            self.assertEqual(rows[1]["query_id"], "Q2")
            self.assertEqual(rows[1]["matches"], [])
            results.append(rows)
        self.assertEqual(*results)

    def test_explicit_uncompressed_annotations_and_query_limit(self):
        path = self.root / "custom.jsonl"
        self.write_annotations(path)
        _, stdout, _ = self.invoke(self.args() + ["--annotations", str(path), "--d", "0", "--limit", "1"])
        rows = [json.loads(line) for line in stdout.splitlines()]
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0]["matches"]), 1)
        self.assertEqual(rows[0]["matches"][0]["annotations"], self.records)

    def test_minimal_query_and_absent_annotation_file(self):
        self.queries.write_text("query_id\tsequence\nQ1\tCASSF\n", encoding="utf-8")
        _, stdout, _ = self.invoke(self.args() + ["--d", "0"])
        row = json.loads(stdout)
        self.assertEqual(row["query_metadata"], {})
        self.assertEqual(row["matches"][0]["annotations"], [])

    def test_explicit_missing_annotation_file_fails_before_writing_output(self):
        output = self.root / "result.jsonl"
        with self.assertRaises(SystemExit) as raised:
            self.invoke(self.args() + ["--annotations", str(self.root / "missing.gz"),
                                       "--output", str(output)])
        self.assertEqual(raised.exception.code, 2)
        self.assertFalse(output.exists())

    def test_validation_does_not_load_annotations(self):
        (self.root / "annotations.jsonl.gz").write_bytes(b"invalid annotation file")
        with patch.object(main, "read_annotations", side_effect=AssertionError("not needed")):
            code, _, _ = self.invoke(self.args("validate") + ["--output", str(self.root / "validate.csv")])
        self.assertEqual(code, 0)

    def test_no_arguments_show_help_and_search_requires_input_paths(self):
        code, stdout, _ = self.invoke([])
        self.assertEqual(code, 0)
        self.assertIn("--reference", stdout)
        with self.assertRaises(SystemExit) as raised:
            self.invoke(["search"])
        self.assertEqual(raised.exception.code, 2)


if __name__ == "__main__":
    unittest.main()
