"""Compare the exhaustive baseline with the group's indexed search.

Run from the GitHub repository root: python -m Exhaustive_Benchmark.compare
"""

import argparse
import csv
import gc
import json
import math
import platform
import statistics
import sys
import time
import tracemalloc
from collections import defaultdict
from pathlib import Path

from Inverted_Index.inverted_index import InvertedIndex, read_sequences
from verification_basic import verification_basic
from verification import verify

from .baseline import exhaustive_matches


def verifier_basic(query, reference, d):
    """Adapt the full-matrix verifier to the common three-argument interface."""
    distance = verification_basic(query, reference)
    return distance if distance <= d else None


def verify_candidates(query, references, candidate_ids, d, verifier):
    matches = {}
    for rid in candidate_ids:
        distance = verifier(query, references[rid], d)
        if distance is not None:
            if type(distance) is not int or distance < 0:
                raise ValueError("verifier must return a non-negative integer or None")
            if distance <= d:
                matches[rid] = distance
    return matches


def indexed_matches(index, query, d, verifier, strategy="spread"):
    candidates = index.generate_candidates(query, d, strategy)
    return verify_candidates(query, index.references, candidates.candidate_ids, d, verifier)


def match_differences(expected, observed):
    missing = sorted(expected.keys() - observed.keys())
    unexpected = sorted(observed.keys() - expected.keys())
    wrong_distance = sorted(
        (rid, expected[rid], observed[rid])
        for rid in expected.keys() & observed.keys()
        if expected[rid] != observed[rid]
    )
    return missing, unexpected, wrong_distance


def load_truth(path, summary_path, reference_set, queries, d, reference_count):
    """Load truth only after confirming that all query/threshold cases exist."""
    with open(summary_path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference_set", "threshold", "query_id", "n_matches", "n_reference"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("truth summary is missing required columns")
        summaries = {}
        for row in reader:
            if row["reference_set"] == reference_set and int(row["threshold"]) == d:
                qid = row["query_id"]
                if qid in summaries:
                    raise ValueError(f"duplicate truth summary for {qid}")
                summaries[qid] = (int(row["n_matches"]), int(row["n_reference"]))
    missing_queries = set(queries) - summaries.keys()
    if missing_queries:
        raise ValueError(f"truth summary does not cover {len(missing_queries)} query IDs")
    if any(summaries[qid][1] != reference_count for qid in queries):
        raise ValueError("truth summary reference count does not match the input reference file")
    truth = defaultdict(dict)
    with open(path, encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"reference_set", "threshold", "query_id", "reference_id", "distance"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("truth matches is missing required columns")
        for row in reader:
            if row["reference_set"] == reference_set and int(row["threshold"]) == d:
                qid, rid = row["query_id"], row["reference_id"]
                if qid in queries:
                    if rid in truth[qid]:
                        raise ValueError(f"duplicate truth reference ID {rid} for {qid}")
                    truth[qid][rid] = int(row["distance"])
    for qid in queries:
        if len(truth[qid]) != summaries[qid][0]:
            raise ValueError(f"truth match count does not agree with summary for {qid}")
    return truth


def validation_rows(references, queries, index, d, strategy, truth=None):
    """Detailed comparisons; no-match queries produce rows as well."""
    rows = []
    for query_id, query in queries.items():
        baseline = exhaustive_matches(query, references, d)
        candidate = index.generate_candidates(query, d, strategy)
        comparisons = [
            ("indexed_basic", baseline,
             verify_candidates(query, references, candidate.candidate_ids, d, verifier_basic)),
            ("indexed_optimized", baseline,
             verify_candidates(query, references, candidate.candidate_ids, d, verify)),
        ]
        if truth is not None:
            comparisons.append(("exhaustive_vs_supplied_truth", truth[query_id], baseline))
        for name, expected, observed in comparisons:
            missing, unexpected, wrong = match_differences(expected, observed)
            rows.append({
                "query_id": query_id, "query_sequence": query, "d": d,
                "comparison": name, "expected_matches": len(expected),
                "observed_matches": len(observed),
                "match": not (missing or unexpected or wrong),
                "missing_ids": json.dumps(missing),
                "unexpected_ids": json.dumps(unexpected),
                "wrong_distances": json.dumps(wrong),
                "candidate_count": len(candidate.candidate_ids),
                "length_eligible_count": candidate.length_eligible_count,
                "index_mode": candidate.mode,
            })
    return rows


def batch_time(queries, operation):
    start = time.perf_counter()
    for query in queries.values():
        operation(query)
    return (time.perf_counter() - start) * 1000 / len(queries)


def benchmark(references, queries, d, k, strategy, repeats):
    """Run paired whole-search timings; verify correctness before reporting."""
    if not references or not queries or repeats < 2:
        raise ValueError("benchmark needs non-empty inputs and at least two repeats")
    build_ms = []
    for _ in range(repeats):
        start = time.perf_counter()
        index = InvertedIndex(references, k)
        build_ms.append((time.perf_counter() - start) * 1000)

    # Correctness and candidate statistics are outside the timed search loops.
    details = validation_rows(references, queries, index, d, strategy)
    failures = [row for row in details if not row["match"]]
    if failures:
        raise ValueError(f"Correctness check failed for {len(failures)} query/method cases; run validate")

    candidates = {qid: index.generate_candidates(q, d, strategy)
                  for qid, q in queries.items()}
    operations = {
        "exhaustive_basic": lambda q: exhaustive_matches(q, references, d),
        "indexed_basic": lambda q: indexed_matches(index, q, d, verifier_basic, strategy),
        "indexed_optimized": lambda q: indexed_matches(index, q, d, verify, strategy),
    }
    first_query = next(iter(queries.values()))
    for operation in operations.values():
        operation(first_query)
    samples = {name: [] for name in operations}
    names = list(operations)
    for repeat in range(repeats):
        # Rotate the order so the same method is not always measured first.
        for name in names[repeat % len(names):] + names[:repeat % len(names)]:
            samples[name].append(batch_time(queries, operations[name]))

    # Compare verifiers on exactly the same precomputed candidate pairs.
    pairs = [(q, references[rid]) for qid, q in queries.items()
             for rid in sorted(candidates[qid].candidate_ids)]
    verifier_samples = {"candidate_basic": [], "candidate_optimized": []}
    if pairs:
        for repeat in range(repeats):
            variants = [("candidate_basic", verifier_basic), ("candidate_optimized", verify)]
            if repeat % 2:
                variants.reverse()
            for name, verifier in variants:
                start = time.perf_counter()
                for query, sequence in pairs:
                    verifier(query, sequence, d)
                verifier_samples[name].append((time.perf_counter() - start) * 1000 / len(pairs))

    eligible_total = sum(c.length_eligible_count for c in candidates.values())
    candidate_total = sum(len(c.candidate_ids) for c in candidates.values())
    metadata = {
        "reference_count": len(references), "query_count": len(queries),
        "query_ids": list(queries), "d": d, "k": k, "strategy": strategy,
        "repeats": repeats, "fallback_count": sum(c.mode == "fallback" for c in candidates.values()),
        "candidate_count_total": candidate_total, "length_eligible_count_total": eligible_total,
        "dp_reduction_fraction": 1 - candidate_total / eligible_total if eligible_total else None,
        "index_build_ms_samples": build_ms,
        "index_build_median_ms": statistics.median(build_ms),
        "python_version": platform.python_version(), "platform": platform.platform(),
        "timing_scope": "Data loading, validation and output excluded; whole-search includes result construction. Index build reported separately.",
        "memory_scope": "Memory is not measured by this benchmark.",
    }
    rows = []
    for name, values in list(samples.items()) + list(verifier_samples.items()):
        rows.append({
            "method": name,
            "unit": "ms_per_candidate_pair" if name.startswith("candidate_") else "ms_per_query",
            "observations_per_repeat": len(pairs) if name.startswith("candidate_") else len(queries),
            "mean": statistics.mean(values) if values else None,
            "median": statistics.median(values) if values else None,
            "standard_error": statistics.stdev(values) / math.sqrt(len(values)) if len(values) >= 2 else None,
            "samples": json.dumps(values),
        })
    return metadata, rows


def memory_comparison(references, queries, d, k, strategy):
    """Measure incremental peak Python allocations in separate tracing runs.

    The reference and query data are already loaded before tracing begins.
    Indexed runs include one index build plus searching the full query batch.
    These figures are not process RSS and do not include preloaded data.
    """
    if not references or not queries:
        raise ValueError("memory comparison requires non-empty inputs")
    index = InvertedIndex(references, k)
    failures = [row for row in validation_rows(references, queries, index, d, strategy)
                if not row["match"]]
    if failures:
        raise ValueError(f"Correctness check failed for {len(failures)} query/method cases; run validate")
    del index
    rows = []
    for method in ("exhaustive_basic", "indexed_basic", "indexed_optimized"):
        gc.collect()
        if tracemalloc.is_tracing():
            raise RuntimeError("memory comparison requires tracemalloc to be stopped initially")
        tracemalloc.start()
        try:
            index = InvertedIndex(references, k) if method.startswith("indexed_") else None
            for query in queries.values():
                if method == "exhaustive_basic":
                    result = exhaustive_matches(query, references, d)
                else:
                    verifier = verifier_basic if method == "indexed_basic" else verify
                    result = indexed_matches(index, query, d, verifier, strategy)
                del result
            current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        rows.append({"method": method, "reference_count": len(references),
                     "query_count": len(queries), "d": d, "k": k,
                     "strategy": strategy, "incremental_peak_python_bytes": peak,
                     "python_bytes_retained_at_end": current,
                     "scope": "Python allocations after inputs are loaded; indexed methods include build and query batch; not RSS"})
        del index
    return rows


def write_csv(path, rows):
    if not rows:
        raise ValueError("no rows to write")
    with open(path, "x", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def positive_int(value):
    parsed = int(value)
    if parsed < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return parsed


def nonnegative_int(value):
    parsed = int(value)
    if parsed < 0:
        raise argparse.ArgumentTypeError("must be a non-negative integer")
    return parsed


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("validate", "benchmark", "memory"))
    parser.add_argument("--reference", required=True, type=Path)
    parser.add_argument("--queries", required=True, type=Path)
    parser.add_argument("--d", required=True, type=nonnegative_int)
    parser.add_argument("--k", default=3, type=positive_int)
    parser.add_argument("--strategy", choices=("spread", "prefix"), default="spread")
    parser.add_argument("--limit", type=positive_int, help="Use first N queries; omit to use all")
    parser.add_argument("--repeats", type=positive_int, default=5)
    parser.add_argument("--output", required=True, type=Path, help="New CSV file; never overwritten")
    parser.add_argument("--metadata", type=Path, help="New JSON metadata file; required for benchmark")
    parser.add_argument("--truth", type=Path)
    parser.add_argument("--truth-summary", type=Path)
    parser.add_argument("--reference-set")
    args = parser.parse_args(argv)
    if args.mode == "benchmark" and (args.metadata is None or args.repeats < 2):
        parser.error("benchmark requires --metadata and --repeats >= 2")
    if args.mode in ("benchmark", "memory") and any((args.truth, args.truth_summary, args.reference_set)):
        parser.error("truth options are only supported in validate mode")
    if args.mode == "validate" and any((args.truth, args.truth_summary, args.reference_set)) and not all((args.truth, args.truth_summary, args.reference_set)):
        parser.error("truth comparison requires --truth, --truth-summary and --reference-set")
    if args.mode in ("validate", "memory") and args.metadata is not None:
        parser.error("--metadata is only used in benchmark mode")
    if args.output.exists() or (args.metadata and args.metadata.exists()):
        parser.error("output files already exist; choose new paths")
    if args.metadata and args.output.resolve() == args.metadata.resolve():
        parser.error("output and metadata paths must be different")

    references = read_sequences(args.reference, "reference_id")
    queries = read_sequences(args.queries, "query_id")
    if args.limit is not None:
        queries = dict(list(queries.items())[:args.limit])
    if not references or not queries:
        parser.error("reference and query files must not be empty")
    if args.mode == "validate":
        index = InvertedIndex(references, args.k)
        truth = (load_truth(args.truth, args.truth_summary, args.reference_set,
                            queries, args.d, len(references)) if args.truth else None)
        rows = validation_rows(references, queries, index, args.d, args.strategy, truth)
        write_csv(args.output, rows)
        failures = sum(not row["match"] for row in rows)
        print(f"Checked {len(rows)} query comparisons; mismatches: {failures}. Details: {args.output}")
        return 1 if failures else 0
    if args.mode == "memory":
        rows = memory_comparison(references, queries, args.d, args.k, args.strategy)
        write_csv(args.output, rows)
        print(f"Incremental Python allocation peaks: {args.output}; these are not process RSS")
        return 0
    metadata, rows = benchmark(references, queries, args.d, args.k, args.strategy, args.repeats)
    metadata.update(reference_path=str(args.reference), query_path=str(args.queries),
                    reference_set=args.reference_set)
    write_csv(args.output, rows)
    with args.metadata.open("x", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)
    print(f"Benchmark results: {args.output}; metadata: {args.metadata}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
