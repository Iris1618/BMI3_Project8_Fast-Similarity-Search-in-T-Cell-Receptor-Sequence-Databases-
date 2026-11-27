"""Check saved-report comparison, including incompatible experiment settings."""

import contextlib
import csv
import io
import json
from pathlib import Path
import tempfile
import unittest

import main
from tcr_search.benchmark import summarize_benchmarks, write_csv


class BenchmarkSummaryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def report(self, percentage, **updates):
        metadata = {
            "reference_count": 1000, "query_count": 2, "query_ids": ["Q1", "Q2"],
            "query_path": str(self.root / "queries.tsv"),
            "reference_path": str(self.root / f"reference_gqe_{percentage}pct.tsv"),
            "d": 2, "k": 3, "strategy": "spread", "repeats": 5,
            "python_version": "3.13", "platform": "test", "verifier": "same DP",
            "timing_scope": "search", "memory_scope": "Python allocations",
            "candidate_count_total": 10, "length_eligible_count_total": 2000,
            "index_build_median_ms": 4,
        }
        metadata.update(updates)
        (self.root / f"gqe_{percentage}pct.metadata.json").write_text(
            json.dumps(metadata), encoding="utf-8")
        rows = [{"method": method, "unit": "ms_per_query", "observations_per_repeat": 2,
                 "median": median, "incremental_peak_python_bytes": memory}
                for method, median, memory in (("exhaustive", 10, 1024), ("indexed", 2, 1048576))]
        path = self.root / f"gqe_{percentage}pct.csv"
        write_csv(path, rows)
        return path

    def test_summary_orders_percentages_and_computes_comparison_metrics(self):
        for percentage in (90, 0, 50):
            self.report(percentage)
        (self.root / "summary.csv").write_text("old summary", encoding="utf-8")
        rows = summarize_benchmarks(self.root)
        self.assertEqual([r["gqe_percentage"] for r in rows], [0, 50, 90])
        row = rows[0]
        self.assertEqual(row["search_speedup"], 5)
        self.assertEqual(row["average_candidates_per_query"], 5)
        self.assertAlmostEqual(row["dp_reduction_percent"], 99.5)
        self.assertEqual(row["indexed_incremental_peak_python_mib"], 1)
        self.assertEqual(Path(row["benchmark_csv"]).name, "gqe_0pct.csv")

    def test_different_settings_and_query_inputs_are_rejected(self):
        for change in ({"d": 1}, {"reference_count": 5000}, {"query_ids": ["Q1", "Q3"]}):
            with self.subTest(change=change), tempfile.TemporaryDirectory() as temp:
                self.root = Path(temp)
                self.report(0)
                self.report(50, **change)
                with self.assertRaisesRegex(ValueError, "differ"):
                    summarize_benchmarks(self.root)

    def test_zero_time_and_no_eligible_references_leave_ratios_empty(self):
        csv_path = self.report(0, candidate_count_total=0, length_eligible_count_total=0)
        with csv_path.open(encoding="utf-8-sig", newline="") as file:
            rows = list(csv.DictReader(file))
        rows[1]["median"] = 0
        with csv_path.open("w", encoding="utf-8-sig", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)
        summary = summarize_benchmarks(self.root)[0]
        self.assertIsNone(summary["search_speedup"])
        self.assertIsNone(summary["dp_reduction_percent"])
        self.assertEqual(summary["average_candidates_per_query"], 0)

    def test_missing_partner_and_invalid_method_rows_are_rejected(self):
        path = self.report(0)
        path.write_text("method\nindexed\nindexed\n", encoding="utf-8")
        with self.assertRaisesRegex(ValueError, "one indexed and one exhaustive"):
            summarize_benchmarks(self.root)
        path.unlink()
        with self.assertRaises(FileNotFoundError):
            summarize_benchmarks(self.root)

    def test_empty_directory_has_clear_error(self):
        with self.assertRaisesRegex(ValueError, "no benchmark metadata"):
            summarize_benchmarks(self.root)

    def test_cli_creates_only_requested_output_and_preserves_sources(self):
        self.report(0)
        self.report(50)
        before = {p.name: p.read_bytes() for p in self.root.iterdir()}
        output = self.root / "summary.csv"
        args = ["summarize", "--input-dir", str(self.root), "--output", str(output)]
        with contextlib.redirect_stdout(io.StringIO()):
            self.assertEqual(main.main(args), 0)
        self.assertEqual(set(p.name for p in self.root.iterdir()), set(before) | {"summary.csv"})
        for name, contents in before.items():
            self.assertEqual((self.root / name).read_bytes(), contents)
        with output.open(encoding="utf-8-sig", newline="") as file:
            self.assertEqual(len(list(csv.DictReader(file))), 2)
        saved = output.read_bytes()
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit) as raised:
            main.main(args)
        self.assertEqual(raised.exception.code, 2)
        self.assertEqual(output.read_bytes(), saved)


class BatchBenchmarkTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.queries = self.root / "queries.tsv"
        self.queries.write_text("query_id\tsequence\nQ1\tCASSF\nQ2\tYYYYF\n", encoding="utf-8")
        self.references = []
        for percentage in (0, 50, 90):
            path = self.root / f"reference_gqe_{percentage}pct.tsv"
            path.write_text("reference_id\tsequence\nR1\tCASSF\nR2\tCASAF\n", encoding="utf-8")
            self.references.append(path)

    def args(self, references=None):
        references = self.references if references is None else references
        return ["benchmark", "--queries", str(self.queries), "--reference",
                *map(str, references), "--d", "1", "--k", "2", "--repeats", "2"]

    def test_batch_from_missing_output_directory_creates_reports_and_summary(self):
        output_dir = self.root / "results" / "comparison"
        sources = {path: path.read_bytes() for path in [self.queries, *self.references]}
        with contextlib.redirect_stdout(io.StringIO()):
            result = main.main(self.args() + ["--output-dir", str(output_dir), "--limit", "1"])
        self.assertEqual(result, 0)
        expected = {f"gqe_{pct}pct{suffix}" for pct in (0, 50, 90)
                    for suffix in (".csv", ".metadata.json")} | {"summary.csv"}
        self.assertEqual({path.name for path in output_dir.iterdir()}, expected)
        with (output_dir / "summary.csv").open(encoding="utf-8-sig", newline="") as file:
            rows = list(csv.DictReader(file))
        self.assertEqual([row["gqe_percentage"] for row in rows], ["0", "50", "90"])
        self.assertTrue(all(row["query_count"] == "1" for row in rows))
        for percentage, row in zip((0, 50, 90), rows):
            with (output_dir / f"gqe_{percentage}pct.csv").open(encoding="utf-8-sig", newline="") as file:
                methods = {item["method"]: item for item in csv.DictReader(file)}
            expected_speedup = float(methods["exhaustive"]["median"]) / float(methods["indexed"]["median"])
            self.assertAlmostEqual(float(row["search_speedup"]), expected_speedup)
        for path, contents in sources.items():
            self.assertEqual(path.read_bytes(), contents)

    def test_single_reference_keeps_csv_and_custom_metadata_output(self):
        output = self.root / "single.csv"
        metadata = self.root / "custom.json"
        with contextlib.redirect_stdout(io.StringIO()):
            result = main.main(self.args(self.references[:1]) +
                               ["--output", str(output), "--metadata", str(metadata)])
        self.assertEqual(result, 0)
        self.assertTrue(output.is_file())
        self.assertTrue(metadata.is_file())
        self.assertFalse((self.root / "summary.csv").exists())

    def test_invalid_batch_options_do_not_write_reports(self):
        output_dir = self.root / "existing"
        output_dir.mkdir()
        existing = output_dir / "keep.txt"
        existing.write_text("preserve", encoding="utf-8")
        cases = [self.args(), self.args() + ["--output-dir", str(output_dir)],
                 self.args() + ["--output-dir", str(self.root / "new"), "--metadata", "custom.json"]]
        for args in cases:
            with self.subTest(args=args), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as raised:
                    main.main(args)
                self.assertEqual(raised.exception.code, 2)
        self.assertEqual(existing.read_text(), "preserve")
        self.assertFalse((self.root / "new").exists())


if __name__ == "__main__":
    unittest.main()
