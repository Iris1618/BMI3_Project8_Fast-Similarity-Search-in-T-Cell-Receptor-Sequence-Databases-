"""Validate and benchmark indexed search against exhaustive search.

Adapted from Exhaustive_Benchmark/compare.py at upstream commit 76316aa.
Both methods use the final verifier in verification.py. The application passes
its indexed-search function in, keeping this module independent of main.py.

References for three packages:
https://docs.python.org/3.12/library/time.html#time.perf_counter
https://docs.python.org/3.12/library/tracemalloc.html
https://docs.python.org/3.12/library/gc.html#gc.collect
"""

import csv
import gc
import json
import math
import platform
import statistics
import time
import tracemalloc

from .exhaustive_search import exhaustive_matches
from .inverted_index import InvertedIndex

#checks whether indexed search gives the same results as exhaustive search.
def validation_rows(references, queries, index, d, strategy, indexed_search):
    """Compare exact IDs and distances, including queries with no matches."""
    rows = []
    for query_id, query in queries.items():
        expected = exhaustive_matches(query, references, d)
        observed = indexed_search(index, query, d, strategy)
        candidates = index.generate_candidates(query, d, strategy)
        # Compare reference IDs first, then compare the distances for shared IDs.
        missing = sorted(expected.keys() - observed.keys())
        unexpected = sorted(observed.keys() - expected.keys())
        wrong = sorted(
            (rid, expected[rid], observed[rid])
            for rid in expected.keys() & observed.keys()
            if expected[rid] != observed[rid]
        )
        rows.append({
            "query_id": query_id, "query_sequence": query, "d": d,
            "comparison": "indexed_vs_exhaustive",
            "expected_matches": len(expected), "observed_matches": len(observed),
            "match": not (missing or unexpected or wrong),
            "missing_ids": json.dumps(missing),
            "unexpected_ids": json.dumps(unexpected),
            "wrong_distances": json.dumps(wrong),
            "candidate_count": len(candidates.candidate_ids),
            "length_eligible_count": candidates.length_eligible_count,
            "index_mode": candidates.mode,
        })
    return rows


#checks that the two search methods agree before the benchmark continues.
def _require_matching_results(references, queries, index, d, strategy, indexed_search):
    rows = validation_rows(references, queries, index, d, strategy, indexed_search)
    failures = sum(not row["match"] for row in rows)
    if failures:
        raise ValueError(
            f"Correctness check failed for {failures} queries; run validate"
        )


#measure the time to process a batch of queries for a method
def _batch_time(queries, operation):
    # perf_counter measures elapsed time; the difference is in seconds.
    start = time.perf_counter()
    for query in queries.values():
        operation(query)
    # Convert seconds to milliseconds and report the average for one query.
    return (time.perf_counter() - start) * 1000 / len(queries)


#the main comparison between exhaustive search and indexed search
def benchmark(references, queries, d, k, strategy, repeats, indexed_search):
    """Validate results, time both methods, then measure memory separately."""
    if not references or not queries or type(repeats) is not int or repeats < 2:
        raise ValueError("benchmark needs non-empty inputs and at least two repeats")

    index = InvertedIndex(references, k)
    _require_matching_results(references, queries, index, d, strategy, indexed_search)

    build_ms = []
    for _ in range(repeats):
        start = time.perf_counter()
        index = InvertedIndex(references, k)
        build_ms.append((time.perf_counter() - start) * 1000)

    # Candidate statistics and validation are outside timed search loops.
    candidates = {
        qid: index.generate_candidates(query, d, strategy)
        for qid, query in queries.items()
    }
    operations = {
        "exhaustive": lambda query: exhaustive_matches(query, references, d),
        "indexed": lambda query: indexed_search(index, query, d, strategy),
    }
    # Both methods run once before measurement to warm up.
    first_query = next(iter(queries.values()))
    for operation in operations.values():
        operation(first_query)
    samples = {name: [] for name in operations}
    names = list(operations)
    for repeat in range(repeats):
        # Alternate execution order while keeping the query batch identical.
        order = names if repeat % 2 == 0 else names[::-1]
        for name in order:
            samples[name].append(_batch_time(queries, operations[name]))

    eligible_total = sum(c.length_eligible_count for c in candidates.values())
    candidate_total = sum(len(c.candidate_ids) for c in candidates.values())
    # Compare the number of DP checks with length-filtered exhaustive search.
    metadata = {
        "reference_count": len(references), "query_count": len(queries),
        "query_ids": list(queries), "d": d, "k": k, "strategy": strategy,
        "repeats": repeats,
        "fallback_count": sum(c.mode == "fallback" for c in candidates.values()),
        "candidate_count_total": candidate_total,
        "length_eligible_count_total": eligible_total,
        "dp_reduction_fraction": 1 - candidate_total / eligible_total if eligible_total else None,
        "index_build_ms_samples": build_ms,
        "index_build_median_ms": statistics.median(build_ms),
        "python_version": platform.python_version(), "platform": platform.platform(),
        "verifier": "verification.verify for both methods",
        "timing_scope": (
            "Data loading, validation and output excluded; whole-search includes "
            "candidate generation and result construction. Index build reported separately."
        ),
        "memory_scope": (
            "Incremental Python allocations after inputs are loaded; indexed includes "
            "one index build and the query batch; not total process memory or RSS. "
            "Measured in separate tracing runs after timing, outside all timing samples."
        ),
    }
    # SE describes variation between repeated timings of the same query batch.
    rows = [{
        "method": name, "unit": "ms_per_query",
        "observations_per_repeat": len(queries),
        "mean": statistics.mean(values), "median": statistics.median(values),
        "standard_error": statistics.stdev(values) / math.sqrt(len(values)),
        "samples": json.dumps(values),
    } for name, values in samples.items()]
    # Memory is measured separately so tracing does not affect the timing samples.
    memory_by_method = {
        row["method"]: row
        for row in memory_comparison(references, queries, d, k, strategy, indexed_search)
    }
    for row in rows:
        memory = memory_by_method[row["method"]]
        row["incremental_peak_python_bytes"] = memory["incremental_peak_python_bytes"]
        row["python_bytes_retained_at_end"] = memory["python_bytes_retained_at_end"]
    return metadata, rows


def memory_comparison(references, queries, d, k, strategy, indexed_search):
    """Measure incremental Python allocations after loading the input data.

    Each method uses a separate tracing run. Indexed search includes one index
    build and the query batch. These peaks exclude preloaded inputs and are not
    total process memory or RSS.
    """
    if not references or not queries:
        raise ValueError("memory comparison requires non-empty inputs")
    if tracemalloc.is_tracing():
        raise RuntimeError("memory comparison requires tracemalloc to be stopped initially")
    index = InvertedIndex(references, k)
    _require_matching_results(references, queries, index, d, strategy, indexed_search)
    del index

    rows = []
    for method in ("exhaustive", "indexed"):
        # Request a collection before each run; this does not guarantee that
        # memory is returned to the operating system or that both runs match.
        gc.collect()
        # Inputs are already loaded, so only later Python allocations are traced.
        tracemalloc.start()
        try:
            index = InvertedIndex(references, k) if method == "indexed" else None
            for query in queries.values():
                if method == "exhaustive":
                    result = exhaustive_matches(query, references, d)
                else:
                    result = indexed_search(index, query, d, strategy)
                del result
            # These values are traced Python allocations, not total process RSS.
            current, peak = tracemalloc.get_traced_memory()
        finally:
            # Stop tracing even if a search raises an error.
            tracemalloc.stop()
        rows.append({
            "method": method, "reference_count": len(references),
            "query_count": len(queries), "d": d, "k": k, "strategy": strategy,
            "incremental_peak_python_bytes": peak,
            "python_bytes_retained_at_end": current,
            "scope": (
                "Python allocations after inputs are loaded; indexed includes "
                "build and query batch; not RSS"
            ),
        })
        del index
    return rows


def write_csv(path, rows):
    """Write a new results CSV without overwriting an existing file."""
    if not rows:
        raise ValueError("no rows to write")
    with open(path, "x", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
