"""Validate and benchmark indexed search against exhaustive search.

Adapted from Exhaustive_Benchmark/compare.py at upstream commit 76316aa.
Both methods use the final verifier in verification.py. The application passes
its indexed-search function in, keeping this module independent of main.py.
"""

import csv
import gc
import json
import math
import platform
import statistics
import time
import tracemalloc
from pathlib import Path

from .exhaustive_search import exhaustive_matches
from .inverted_index import InvertedIndex

COMPARISONS = {
    "gqe": (),
    "size": ("reference_count",),
    "length": ("query_path", "query_ids"),
    "distance": ("d",),
    "seed_length": ("k",),
}


def validation_rows(references, queries, index, d, indexed_search):
    """Compare exact IDs and distances, including queries with no matches."""
    rows = []
    for query_id, query in queries.items():
        expected = exhaustive_matches(query, references, d)
        observed = indexed_search(index, query, d)
        candidates = index.generate_candidates(query, d)
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


def _require_matching_results(references, queries, index, d, indexed_search):
    rows = validation_rows(references, queries, index, d, indexed_search)
    failures = sum(not row["match"] for row in rows)
    if failures:
        raise ValueError(
            f"Correctness check failed for {failures} queries; run validate"
        )


def _batch_time(queries, operation):
    start = time.perf_counter()
    for query in queries.values():
        operation(query)
    return (time.perf_counter() - start) * 1000 / len(queries)


def benchmark(references, queries, d, k, repeats, indexed_search):
    """Validate results, time both methods, then measure memory separately."""
    if not references or not queries or type(repeats) is not int or repeats < 2:
        raise ValueError("benchmark needs non-empty inputs and at least two repeats")

    index = InvertedIndex(references, k)
    _require_matching_results(references, queries, index, d, indexed_search)

    build_ms = []
    for _ in range(repeats):
        start = time.perf_counter()
        index = InvertedIndex(references, k)
        build_ms.append((time.perf_counter() - start) * 1000)

    # Candidate statistics and validation are outside timed search loops.
    candidates = {
        qid: index.generate_candidates(query, d)
        for qid, query in queries.items()
    }
    operations = {
        "exhaustive": lambda query: exhaustive_matches(query, references, d),
        "indexed": lambda query: indexed_search(index, query, d),
    }
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
    metadata = {
        "reference_count": len(references), "query_count": len(queries),
        "average_reference_length": statistics.mean(map(len, references.values())),
        "average_query_length": statistics.mean(map(len, queries.values())),
        "query_ids": list(queries), "d": d, "k": k,
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
    rows = [{
        "method": name, "unit": "ms_per_query",
        "observations_per_repeat": len(queries),
        "mean": statistics.mean(values), "median": statistics.median(values),
        "standard_error": statistics.stdev(values) / math.sqrt(len(values)),
        "samples": json.dumps(values),
    } for name, values in samples.items()]
    memory_by_method = {
        row["method"]: row
        for row in memory_comparison(references, queries, d, k, indexed_search)
    }
    for row in rows:
        memory = memory_by_method[row["method"]]
        row["incremental_peak_python_bytes"] = memory["incremental_peak_python_bytes"]
        row["python_bytes_retained_at_end"] = memory["python_bytes_retained_at_end"]
    return metadata, rows


def memory_comparison(references, queries, d, k, indexed_search):
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
    _require_matching_results(references, queries, index, d, indexed_search)
    del index

    rows = []
    for method in ("exhaustive", "indexed"):
        gc.collect()
        tracemalloc.start()
        try:
            index = InvertedIndex(references, k) if method == "indexed" else None
            for query in queries.values():
                if method == "exhaustive":
                    result = exhaustive_matches(query, references, d)
                else:
                    result = indexed_search(index, query, d)
                del result
            current, peak = tracemalloc.get_traced_memory()
        finally:
            tracemalloc.stop()
        rows.append({
            "method": method, "reference_count": len(references),
            "query_count": len(queries), "d": d, "k": k,
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


def summarize_benchmarks(input_dir, vary="gqe"):
    """Return one summary row for each saved benchmark run."""
    metadata_files = sorted(Path(input_dir).glob("*.metadata.json"))
    if not metadata_files:
        raise ValueError(f"no benchmark metadata files found in {input_dir}")

    summary_rows = []
    first_metadata = None
    settings = ("query_path", "query_ids", "query_count", "reference_count",
                "d", "k", "repeats")
    settings = [setting for setting in settings if setting not in COMPARISONS[vary]]
    if vary in ("distance", "seed_length"):
        settings.append("reference_path")

    for metadata_file in metadata_files:
        # Read the metadata and its matching benchmark CSV.
        with metadata_file.open(encoding="utf-8-sig") as file:
            metadata = json.load(file)
        csv_file = metadata_file.with_name(
            metadata_file.name.removesuffix(".metadata.json") + ".csv")
        with csv_file.open(encoding="utf-8-sig", newline="") as file:
            method_rows = list(csv.DictReader(file))
        methods = {row["method"]: row for row in method_rows}
        if len(method_rows) != 2 or set(methods) != {"exhaustive", "indexed"}:
            raise ValueError(f"{csv_file}: expected one indexed and one exhaustive row")

        # Compare runs that use the same queries and experiment settings.
        if first_metadata is None:
            first_metadata = metadata
        for setting in settings:
            if metadata[setting] != first_metadata[setting]:
                raise ValueError(f"{metadata_file}: benchmark inputs differ in {setting}")

        # Calculate speedup and the average number of candidates per query.
        exhaustive_time = float(methods["exhaustive"]["median"])
        indexed_time = float(methods["indexed"]["median"])
        query_count = metadata["query_count"]
        candidate_count = metadata["candidate_count_total"]
        speedup = exhaustive_time / indexed_time if indexed_time else None
        exhaustive_memory = int(methods["exhaustive"]["incremental_peak_python_bytes"])
        indexed_memory = int(methods["indexed"]["incremental_peak_python_bytes"])
        memory_ratio = indexed_memory / exhaustive_memory if exhaustive_memory else None

        reference_name = Path(metadata["reference_path"]).stem
        summary_row = {"reference": reference_name}
        if vary == "gqe":
            # This label is only relevant to the GQE frequency experiment.
            percentage = reference_name.removeprefix("reference_gqe_").removesuffix("pct")
            summary_row["gqe_percentage"] = int(percentage) if percentage.isdigit() else None

        summary_row.update({
            "comparison": vary,
            "reference_count": metadata["reference_count"], "query_count": query_count,
            "average_reference_length": metadata.get("average_reference_length"),
            "average_query_length": metadata.get("average_query_length"),
            "d": metadata["d"], "k": metadata["k"],
            "fallback_count": metadata.get("fallback_count"),
            "exhaustive_median_ms_per_query": exhaustive_time,
            "indexed_median_ms_per_query": indexed_time,
            "search_speedup": speedup,
            "average_candidates_per_query": candidate_count / query_count,
            "index_build_median_ms": metadata["index_build_median_ms"],
            "exhaustive_incremental_peak_python_mib": exhaustive_memory / 1024**2,
            "indexed_incremental_peak_python_mib": indexed_memory / 1024**2,
            "memory_ratio": memory_ratio,
        })
        summary_rows.append(summary_row)

    # Sort by the condition changed in this experiment.
    sort_column = {"gqe": "gqe_percentage", "size": "reference_count",
                   "length": "average_reference_length", "distance": "d",
                   "seed_length": "k"}[vary]
    summary_rows.sort(key=lambda row: (
        row[sort_column] is None, row[sort_column] or 0, row["reference"]))
    return summary_rows
