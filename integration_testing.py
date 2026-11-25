"""Integration test entry point for BMI3 Project 8.

Purpose
-------
This file deliberately focuses on cross-module behavior instead of
white-box unit tests. It wires together Input_Source, Inverted_Index,
both edit-distance verifiers, and Exhaustive_Benchmark, then checks that
the combined pipeline still agrees with an exhaustive baseline. Unit
tests for individual modules remain in their own test_*.py files.

What is covered
---------------
1. Example data across every combination of threshold d, k-mer size k,
   and seed strategy.
2. The two public search_one result builders (indexed vs. exhaustive).
3. The real candidate-generation CLI and its overwrite guard.
4. compare.py validate mode with externally supplied truth files.
5. benchmark and memory measurement contracts.
6. Synthetic data generation followed by a full search.
7. VDJdb-style cleaning on a small handcrafted input.

How to run
----------
From the repository root:

    python integration_testing.py -v
    python -m unittest integration_testing -v

Requirements
------------
Python 3.10+ (verification.py uses ``int | None`` annotations).
Input_Source/examples/reference.tsv and queries.tsv must exist.
No real VDJdb download is required: the suite uses the bundled examples
plus deterministic synthetic data. Every write goes to
tempfile.TemporaryDirectory(), so the repository is never modified.
"""

import csv
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

WORKSPACE = Path(__file__).resolve().parent
#The code directory name is normally stable, but a downloaded or cloned copy may carry a suffix. Try the canonical path first and fall back to a glob match so renaming the folder does not break imports.
PROJECT = (
    WORKSPACE
    / "BMI3_code"
    / "BMI3_Project8_Fast-Similarity-Search-in-T-Cell-Receptor-Sequence-Databases--main"
)
if not PROJECT.is_dir():
    matches = sorted((WORKSPACE / "BMI3_code").glob("BMI3_Project8*"))
    if not matches:
        raise RuntimeError("Cannot locate the BMI3_Project8 code directory")
    PROJECT = matches[0]
INDEX_DIR = PROJECT / "Inverted_Index"
""" Add two entries to sys.path:
   PROJECT   -> lets Exhaustive_Benchmark and Input_Source be imported as
                packages;
   INDEX_DIR -> integration_example.py uses a flat ``import inverted_index``; putting its own directory on sys.path reuses the real production implementation without editing it.
"""
for module_path in (str(PROJECT), str(INDEX_DIR)):
    if module_path not in sys.path:
        sys.path.insert(0, module_path)

# Import style mirrors each module's original layout: the index side is flat, while the benchmark and input-preparation sides are packaged.
# Keeping both styles avoids changing production code just to satisfy the test harness.
from Exhaustive_Benchmark.baseline import exhaustive_matches, search_one as baseline_search_one
from Exhaustive_Benchmark.compare import (
    benchmark,
    load_truth,
    main as compare_main,
    memory_comparison,
)
from Input_Source.prepare_real import prepare as prepare_real
from Input_Source.prepare_synthetic import prepare as prepare_synthetic
from integration_example import search_one as index_search_one
from inverted_index import InvertedIndex, read_sequences
from verification import verify
from verification_basic import verification_basic

EXAMPLES = PROJECT / "Input_Source" / "examples"
REFERENCE_TSV = EXAMPLES / "reference.tsv"
QUERIES_TSV = EXAMPLES / "queries.tsv"
STRATEGIES = ("spread", "prefix")


def basic_verifier(query, reference, max_distance):
    """Adapt the basic DP to the same return contract as verify().

    verification_basic always returns an int. verification.verify returns
    None when the distance exceeds max_distance. Normalizing both to the
    None-on-miss convention lets run_indexed() use either implementation
    and keeps assertions focused on the final hit set and distances
    rather than on each verifier's individual API style.

    Args:
        query: Query amino-acid sequence.
        reference: Reference amino-acid sequence.
        max_distance: Inclusive edit-distance threshold.

    Returns:
        The edit distance when it is <= max_distance, otherwise None.
    """
    distance = verification_basic(query, reference)
    return distance if distance <= max_distance else None


def run_indexed(index, query, d, verifier, strategy="spread"):
    """Reproduce the production retrieval path for one query.

    Steps:
        1. Generate candidates from the inverted index.
        2. Run the supplied DP verifier only on those candidates.
        3. Keep hits whose exact distance is <= d.

    Args:
        index: An InvertedIndex built over the reference database.
        query: Query sequence to search for.
        d: Maximum edit distance.
        verifier: Callable ``(query, reference, d) -> int | None``.
        strategy: Seed layout, either "spread" or "prefix".

    Returns:
        A tuple ``(CandidateResult, {reference_id: distance})``.

    Comparing the returned mapping with an exhaustive baseline checks two
    independent properties at once:
        - the candidate set does not lose a true hit (no false negatives);
        - the verifier reports the correct edit distance for every hit.
    """
    candidate = index.generate_candidates(query, d, strategy)
    matches = {}
    for reference_id in candidate.candidate_ids:
        distance = verifier(query, index.references[reference_id], d)
        if distance is not None and distance <= d:
            matches[reference_id] = distance
    return candidate, matches


def read_rows(path):
    """Read a TSV with a header and keep every original column.

    Used when a test needs metadata that read_sequences() intentionally
    discards, such as category, source_id, reason, or source_distance.

    Args:
        path: Path to a UTF-8 TSV file with a header row.

    Returns:
        A list of dictionaries, one per data row.
    """
    with open(path, encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ExamplePipelineTests(unittest.TestCase):
    """End-to-end checks built on the bundled example dataset.

    The example set is deliberately small (14 references, 8 queries) but
    covers normal queries, short-query fallback, zero-match queries, and
    multi-match queries. That makes it practical to cross-check every
    parameter combination instead of sampling a subset.
    """

    @classmethod
    def setUpClass(cls):
        # Load the read-only inputs once per class. read_sequences() validates IDs, duplicates, alphabet, and TSV shape, so a malformed example fails here before any search assertions run.
        cls.references = read_sequences(REFERENCE_TSV, "reference_id")
        cls.queries = read_sequences(QUERIES_TSV, "query_id")

    def test_all_thresholds_and_strategies_match_exhaustive(self):
        """Core correctness: every d / k / seed-strategy combination must
        agree with the exhaustive baseline.

        The test checks four independent layers:
            1. Both verifiers produce exactly the same
               {reference_id: distance} mapping as exhaustive search.
            2. Every true hit is present in the candidate set. Candidates
               may contain extras, but missing hits would be a false
               negative.
            3. Every candidate passes the length-difference filter.
            4. The chosen mode matches the pigeonhole-principle threshold.

        Scope: 3 thresholds x 2 k-mer sizes x 2 strategies x 8 queries =
        96 search cases, each run through both verifiers.
        """
        for d in (0, 1, 2):
            for k in (2, 3):
                # k changes the index contents; d is only applied at query time, during candidate generation and verification.
                index = InvertedIndex(self.references, k)
                self.assertEqual(index.k, k)
                for strategy in STRATEGIES:
                    for query_id, query in self.queries.items():
                        # subTest pinpoints the exact d/k/strategy/query combination if an assertion fails.
                        with self.subTest(d=d, k=k, strategy=strategy, query=query_id):
                            # Independent oracle: length filter plus full DP over the whole database, never touching the inverted index.
                            expected = exhaustive_matches(query, self.references, d)
                            candidate, got_basic = run_indexed(
                                index, query, d, basic_verifier, strategy
                            )
                            _, got_optimized = run_indexed(index, query, d, verify, strategy)

                            # Both verifier implementations must return the same hit set and the same distances as the exhaustive baseline.
                            self.assertEqual(got_basic, expected)
                            self.assertEqual(got_optimized, expected)
                            # The candidate set may contain extras, but it must never miss a true hit.
                            self.assertLessEqual(set(expected), set(candidate.candidate_ids))
                            # Every candidate must satisfy the safe length filter, even if it came from a common k-mer.
                            self.assertTrue(
                                all(
                                    abs(len(self.references[rid]) - len(query)) <= d
                                    for rid in candidate.candidate_ids
                                )
                            )
                            # The pigeonhole guarantee needs d+1 disjoint seeds of length k. 
                            # If the query is shorter, no such partition exists and the index must fall back to a full scan; otherwise it must stay in indexed mode.
                            expected_mode = (
                                "fallback" if len(query) < k * (d + 1) else "indexed"
                            )
                            self.assertEqual(candidate.mode, expected_mode)

    def test_final_result_helpers_agree(self):
        """Final-output contract: indexed and exhaustive search_one must
        agree field by field.

        This test does not compare latency. It verifies that the two
        public result builders return the same query metadata, hits,
        distances, sequences, annotations, and ordering, so switching
        entry points cannot silently change the output consumed by
        downstream integration.
        """
        index = InvertedIndex(self.references, 2)
        # Only one reference carries annotations. If the indexed path drops or reshapes annotations while filtering candidates, the full-dictionary comparison below fails.
        annotations = {"EDGE_R05": [{"source_line": 5, "fields": {}}]}
        for query_id, query in self.queries.items():
            with self.subTest(query=query_id):
                indexed = index_search_one(
                    index, query_id, query, 1, basic_verifier, annotations=annotations
                )
                baseline = baseline_search_one(
                    query_id, query, self.references, 1, annotations
                )
                # Compare the whole dict, including matches, query metadata, distances, sequences, and annotations.
                self.assertEqual(indexed, baseline)
                # Ordering is part of the output contract. Stable ordering avoids meaningless diffs in downstream comparisons.
                self.assertEqual(
                    indexed["matches"],
                    sorted(indexed["matches"], key=lambda hit: hit["reference_id"]),
                )

    def test_candidate_cli_and_overwrite_guard(self):
        """Real CLI integration and output-safety checks.

        Why a subprocess: the CLI sets up sys.path, parses arguments,
        writes JSONL, and emits a JSON summary on stderr. Calling main()
        in-process would bypass those paths, so the test invokes the
        module exactly as a user would.

        The test runs in two phases:
            1. Generate candidates and independently verify them.
            2. Point --output at the input file and require a fail-closed
               error with the input left byte-for-byte intact.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "candidates.jsonl"
            # Same flags as the documented usage. All output goes to a temporary directory, so the repository is untouched.
            command = [
                sys.executable,
                str(PROJECT / "Inverted_Index" / "inverted_index.py"),
                "--reference", str(REFERENCE_TSV),
                "--queries", str(QUERIES_TSV),
                "--k", "2",
                "--d", "1",
                "--strategy", "spread",
                "--output", str(output),
            ]
            completed = subprocess.run(
                command, capture_output=True, text=True, check=True
            )
            rows = [json.loads(line) for line in output.read_text(encoding="utf-8").splitlines()]
            # Phase 1a: the JSONL stream must contain exactly one record per query.
            self.assertEqual(len(rows), len(self.queries))
            # The CLI writes its summary JSON to stderr; the last line is the machine-readable record.
            summary = json.loads(completed.stderr.strip().splitlines()[-1])
            self.assertEqual(summary["query_count"], len(self.queries))
            self.assertEqual(summary["reference_count"], len(self.references))
            # Every query must take exactly one path, indexed or fallback, so the two counters must add up to the query count.
            self.assertEqual(
                summary["indexed_queries"] + summary["fallback_queries"],
                len(self.queries),
            )

            # Phase 1b: rebuild the same index and re-verify the CLI's own candidates. 
            # This catches a CLI that writes plausible-looking but incomplete candidate lists.
            cli_index = InvertedIndex(self.references, 2)
            for row in rows:
                with self.subTest(query=row["query_id"]):
                    query = self.queries[row["query_id"]]
                    expected = exhaustive_matches(query, self.references, row["d"])
                    # The CLI candidates must cover every true hit.
                    self.assertLessEqual(set(expected), set(row["candidate_ids"]))
                    _, verified = run_indexed(
                        cli_index,
                        query,
                        row["d"],
                        verify,
                        row["strategy"],
                    )
                    self.assertEqual(verified, expected)

            # Phase 2: overwrite guard. Point --output at the same file used by --reference. 
            # The CLI must refuse before truncating anything. Copy the example first so the repository file is never touched by this test.
            guarded_reference = root / "guarded_reference.tsv"
            shutil.copyfile(REFERENCE_TSV, guarded_reference)
            original = guarded_reference.read_bytes()
            failed = subprocess.run(
                [
                    sys.executable,
                    str(PROJECT / "Inverted_Index" / "inverted_index.py"),
                    "--reference", str(guarded_reference),
                    "--queries", str(QUERIES_TSV),
                    "--k", "2",
                    "--d", "1",
                    "--output", str(guarded_reference),
                ],
                capture_output=True,
                text=True,
            )
            # The CLI reserves exit code 2 for input/argument errors.
            self.assertEqual(failed.returncode, 2)
            self.assertIn("must not overwrite an input file", failed.stderr)
            # The important part: the guard must fail closed. One changed byte would mean the input was already truncated.
            self.assertEqual(guarded_reference.read_bytes(), original)

    def test_compare_validate_with_generated_truth(self):
        """compare.py validate mode plus external truth-file loading.

        Truth is generated from the exhaustive baseline rather than
        copied from indexed output, so it is independent of candidate
        generation. Because the baseline and the ``indexed_basic`` method
        share verification_basic, the strongest cross-check is
        ``indexed_optimized`` (banded DP) against this truth. The
        ``exhaustive_vs_supplied_truth`` row is expected to match by
        construction; it exercises file loading, reference-set matching,
        and zero-match handling rather than algorithmic independence.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            truth_path = root / "truth_matches.tsv"
            summary_path = root / "truth_summary.tsv"
            output_path = root / "validation.csv"
            reference_set = "examples"

            # truth_matches.tsv: one row per query-reference hit, with the distance produced by exhaustive DP.
            with open(truth_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["reference_set", "threshold", "query_id", "reference_id", "distance"]
                )
                for query_id, query in self.queries.items():
                    for reference_id, distance in sorted(
                        exhaustive_matches(query, self.references, 1).items()
                    ):
                        writer.writerow(
                            [reference_set, 1, query_id, reference_id, distance]
                        )

            # truth_summary.tsv: one row per query. n_matches and n_reference let the loader verify completeness. 
            #Zero-match queries must still write n_matches=0; omitting them would make an incomplete truth file look valid.
            with open(summary_path, "w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(
                    ["reference_set", "threshold", "query_id", "n_matches", "n_reference"]
                )
                for query_id, query in self.queries.items():
                    writer.writerow(
                        [
                            reference_set,
                            1,
                            query_id,
                            len(exhaustive_matches(query, self.references, 1)),
                            len(self.references),
                        ]
                    )

            # Verify the loader on its own first. It should validate the reference set, threshold, reference count, and hit counts, then return {query_id: {reference_id: distance}}.
            loaded = load_truth(
                truth_path, summary_path, reference_set,
                self.queries, 1, len(self.references),
            )
            for query_id, query in self.queries.items():
                self.assertEqual(
                    loaded[query_id], exhaustive_matches(query, self.references, 1)
                )

            # Then run the full validate CLI: two indexed methods plus a comparison against the supplied truth.
            return_code = compare_main(
                [
                    "validate",
                    "--reference", str(REFERENCE_TSV),
                    "--queries", str(QUERIES_TSV),
                    "--d", "1",
                    "--k", "2",
                    "--strategy", "spread",
                    "--truth", str(truth_path),
                    "--truth-summary", str(summary_path),
                    "--reference-set", reference_set,
                    "--output", str(output_path),
                ]
            )
            self.assertEqual(return_code, 0)

            with open(output_path, encoding="utf-8-sig", newline="") as handle:
                rows = list(csv.DictReader(handle))
            # Each query produces three comparisons: indexed_basic, indexed_optimized, and exhaustive_vs_supplied_truth.
            self.assertEqual(len(rows), 3 * len(self.queries))
            self.assertEqual(
                {row["comparison"] for row in rows},
                {
                    "indexed_basic",
                    "indexed_optimized",
                    "exhaustive_vs_supplied_truth",
                },
            )
            for row in rows:
                with self.subTest(query=row["query_id"], comparison=row["comparison"]):
                    # Python booleans are serialized to CSV as "True"/"False".
                    self.assertEqual(row["match"], "True")
                    # All three difference lists must be empty: no missing hits, no unexpected hits, and no wrong distances.
                    self.assertEqual(row["missing_ids"], "[]")
                    self.assertEqual(row["unexpected_ids"], "[]")
                    self.assertEqual(row["wrong_distances"], "[]")
                    self.assertEqual(
                        row["expected_matches"], row["observed_matches"]
                    )

    def test_benchmark_and_memory_modes(self):
        """Contract checks for timing and memory measurement, not performance.

        Wall-clock numbers depend on machine load, so asserting an
        absolute speedup would make the suite flaky. Instead this test
        verifies that the benchmark produced complete, understandable
        structures: metadata, samples per repeat, method names, units,
        and non-empty observations.
        """
        # d=1, k=2, spread, repeats=2. At least two repeats are needed to report dispersion alongside the mean.
        metadata, rows = benchmark(self.references, self.queries, 1, 2, "spread", 2)
        # Metadata must record input size, query order, index-build samples, and the median, so runs on different datasets remain traceable.
        self.assertEqual(metadata["reference_count"], len(self.references))
        self.assertEqual(metadata["query_count"], len(self.queries))
        self.assertEqual(metadata["query_ids"], list(self.queries))
        self.assertEqual(len(metadata["index_build_ms_samples"]), 2)
        self.assertGreater(metadata["index_build_median_ms"], 0)
        # Exactly five methods: three whole-batch searches plus two verifier comparisons run on the same candidate pairs.
        self.assertEqual(
            {row["method"] for row in rows},
            {
                "exhaustive_basic",
                "indexed_basic",
                "indexed_optimized",
                "candidate_basic",
                "candidate_optimized",
            },
        )
        for row in rows:
            with self.subTest(method=row["method"]):
                # Every repeat must leave one raw sample so dispersion can be audited later.
                self.assertEqual(len(json.loads(row["samples"])), 2)
                if row["method"].startswith("candidate_"):
                    # Candidate-pair methods are timed per pair, so their
                    # unit differs from whole-batch searches.
                    self.assertGreater(row["observations_per_repeat"], 0)
                    self.assertIsNotNone(row["median"])
                else:
                    self.assertEqual(row["unit"], "ms_per_query")
                    self.assertIsNotNone(row["mean"])
                    self.assertIsNotNone(row["median"])

        # Memory comparison uses tracemalloc and includes index-build cost for the indexed methods. 
        #The metric is incremental Python peak memory, not process RSS, and it excludes data loaded before the trace starts.
        memory_rows = memory_comparison(
            self.references, self.queries, 1, 2, "spread"
        )
        self.assertEqual(
            {row["method"] for row in memory_rows},
            {"exhaustive_basic", "indexed_basic", "indexed_optimized"},
        )
        for row in memory_rows:
            self.assertGreater(row["incremental_peak_python_bytes"], 0)
            # The scope string must state the measurement boundary so the number is not mistaken for process memory.
            self.assertIn("not RSS", row["scope"])


class InputPreparationTests(unittest.TestCase):
    """Checks for the synthetic and real-data input-preparation modules."""

    def test_synthetic_generation_and_search(self):
        """Joint validation of the synthetic generator and search pipeline.

        Three libraries are generated with fixed seeds, so both the data
        and the category labels are reproducible. The test deliberately
        does not trust the generator's own claims: after generation it
        re-runs exhaustive search to verify the generated metadata and the
        indexed results.
        """
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "synthetic"
            # Build only the attributes consumed by prepare(), avoiding argparse and the file-publication wrapper used by main().
            args = SimpleNamespace(
                out=output,
                size=40,
                lengths=[12],
                motif_percentages=[0, 50, 90],
                anchors=3,
                max_edits=2,
                probe_count=2,
                seed=20262111,
                probe_seed=20261113,
            )
            prepare_synthetic(args, output)

            # The top-level report must list all three generated datasets; each dataset also has its own per-library report.
            root_report = json.loads(
                (output / "generation_report.json").read_text(encoding="utf-8")
            )
            self.assertEqual(len(root_report["datasets"]), 3)

            for percentage in (0, 50, 90):
                name = f"L12_motif{percentage:03d}"
                library = output / name
                references = read_sequences(library / "reference.tsv", "reference_id")
                queries = read_sequences(library / "queries.tsv", "query_id")
                # rows keeps metadata such as category and source_id for the category-specific checks below.
                rows = {row["query_id"]: row for row in read_rows(library / "queries.tsv")}

                # Generation must produce the requested library size, and query IDs must be unique enough to match metadata rows.
                self.assertEqual(len(references), 40)
                self.assertEqual(len(queries), len(rows))

                # k=3 with d=1 and d=2. Candidate generation plus banded DP must match the exhaustive baseline exactly.
                index = InvertedIndex(references, 3)
                for d in (1, 2):
                    for query_id, query in queries.items():
                        with self.subTest(library=name, d=d, query=query_id):
                            expected = exhaustive_matches(query, references, d)
                            candidate, got = run_indexed(index, query, d, verify)
                            self.assertEqual(got, expected)
                            self.assertLessEqual(set(expected), set(candidate.candidate_ids))

                for query_id, row in rows.items():
                    if row["category"] == "exact_copy":
                        # Exact copy: the query equals its source, the source ID exists, and even d=0 must already produce the source as a candidate.
                        self.assertEqual(row["sequence"], row["source_sequence"])
                        self.assertIn(row["source_id"], references)
                        self.assertEqual(references[row["source_id"]], row["sequence"])
                        candidate = index.generate_candidates(row["sequence"], 0)
                        self.assertIn(row["source_id"], candidate.candidate_ids)
                    if row["category"] == "short_query":
                        # Short queries (length 1-3) cannot supply d+1 disjoint seeds of length 3, so the index must switch to a full-scan fallback.
                        candidate = index.generate_candidates(row["sequence"], 1)
                        self.assertEqual(candidate.mode, "fallback")

            # A motif probe is a constructed query that is close to, but not equal to, any reference sequence.
            probes_path = output / "motif_probes" / "L12" / "queries.tsv"
            self.assertTrue(probes_path.is_file())
            probes = read_rows(probes_path)
            self.assertEqual(len(probes), 2)
            # Collect sequences across all three libraries so the probe check is not limited to a single motif percentage.
            all_reference_sequences = set()
            for percentage in (0, 50, 90):
                all_reference_sequences.update(
                    read_sequences(
                        output / f"L12_motif{percentage:03d}" / "reference.tsv",
                        "reference_id",
                    ).values()
                )
            for probe in probes:
                with self.subTest(probe=probe["query_id"]):
                    # A probe is guaranteed not to be an exact copy of any reference and to be exactly one edit away from its source.
                    #It is not guaranteed to have zero matches at an arbitrary threshold, so no such claim is asserted.
                    self.assertNotIn(probe["sequence"], all_reference_sequences)
                    self.assertEqual(probe["source_distance"], "1")
                    self.assertEqual(
                        verification_basic(probe["source_sequence"], probe["sequence"]), 1
                    )

    def test_real_data_preparation(self):
        """End-to-end check of the real-data cleaner using a fake VDJdb sample.

        The seven input rows cover duplicate sequences, an invalid amino
        acid, an empty CDR3, a wrong chain, and a wrong species. This
        exercises filtering, normalization, deduplication, annotation
        retention, nested subset generation, and query sampling without
        downloading VDJdb.
        """
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "vdjdb.slim.txt"
            # Seven raw rows: two copies of the same sequence, one invalid character, one empty CDR3, one wrong chain, one wrong species.
            source.write_text(
                "gene\tspecies\tcdr3\tweight\n"
                "TRB\tHomoSapiens\tCASS\t1\n"
                "TRB\tHomoSapiens\tCASS\t2\n"
                "TRB\tHomoSapiens\tCATS\t1\n"
                "TRB\tHomoSapiens\tCASSX\t1\n"
                "TRB\tHomoSapiens\t\t1\n"
                "TRA\tHomoSapiens\tCASS\t1\n"
                "TRB\tMusMusculus\tCASS\t1\n",
                encoding="utf-8",
            )
            output = root / "real"
            args = SimpleNamespace(
                input=source,
                out=output,
                sizes=[2],
                queries=2,
                seed=7,
                species="HomoSapiens",
                gene="TRB",
            )
            prepare_real(args, output)

            # The cleaning report must distinguish rows that passed the species/chain filter from rows that passed CDR3 validation.
            report = json.loads(
                (output / "cleaning_report.json").read_text(encoding="utf-8")
            )
            counts = report["counts"]
            self.assertEqual(counts["raw_rows"], 7)
            self.assertEqual(counts["target_rows"], 5)
            self.assertEqual(counts["valid_rows"], 3)
            self.assertEqual(counts["unique_sequences"], 2)
            self.assertEqual(counts["duplicate_rows_merged"], 1)
            self.assertEqual(counts["missing_cdr3"], 1)
            self.assertEqual(counts["invalid_amino_acid"], 1)
            self.assertEqual(counts["excluded_species_or_chain"], 2)
            # The input fingerprint supports experiment traceability; a SHA-256 hex digest is always 64 characters long.
            self.assertEqual(len(report["input_sha256"]), 64)

            # Deduplication must leave only CASS and CATS. When the smallest subset size equals the full database, the subset and the full reference set must be identical.
            references = read_sequences(output / "reference_all.tsv", "reference_id")
            self.assertEqual(sorted(references.values()), ["CASS", "CATS"])
            subset = read_sequences(output / "reference_n2.tsv", "reference_id")
            self.assertEqual(subset, references)

            # Annotations are grouped by reference ID. The duplicated CASS sequence must retain both original records rather than keeping only the first one.
            with gzip.open(output / "annotations.jsonl.gz", "rt", encoding="utf-8") as handle:
                annotations = [json.loads(line) for line in handle]
            self.assertEqual(len(annotations), 2)
            records_by_id = {row["reference_id"]: row["records"] for row in annotations}
            cass_id = next(rid for rid, sequence in references.items() if sequence == "CASS")
            self.assertEqual(len(records_by_id[cass_id]), 2)

            # Queries should be exact copies of reference sequences with
            # distance zero and a traceable source ID.
            query_rows = read_rows(output / "queries.tsv")
            self.assertEqual(len(query_rows), 2)
            for row in query_rows:
                self.assertEqual(row["category"], "exact_copy")
                self.assertIn(row["source_id"], references)
                self.assertEqual(references[row["source_id"]], row["sequence"])
                self.assertEqual(row["source_distance"], "0")

            # Excluded target rows contain only CDR3 problems. Rows removed by the species/chain filter must not leak into this file, or the cleaning report would double-count them.
            excluded = read_rows(output / "excluded_target_rows.tsv")
            self.assertEqual(
                sorted(row["reason"] for row in excluded),
                ["invalid_amino_acid", "missing_cdr3"],
            )
            # Both length-4 sequences fall outside the recommended 8-30 range, so both must be flagged for manual review.
            self.assertEqual(len(read_rows(output / "length_review_flags.tsv")), 2)


if __name__ == "__main__":
    # Running this file directly uses verbose output so every test name
    # and subtest appears in the console.
    unittest.main(verbosity=2)
