"""Run indexed/exhaustive search and compare their performance."""

import argparse
import json
import sys
from pathlib import Path

# Keep the project folders free of automatically generated __pycache__ files.
sys.dont_write_bytecode = True

from tcr_search.benchmark import benchmark, memory_comparison, summarize_benchmarks, validation_rows, write_csv
from tcr_search.data_input import read_annotations, read_records, read_sequences
from tcr_search.exhaustive_search import exhaustive_matches
from tcr_search.inverted_index import InvertedIndex
from tcr_search.verification import verify

PROJECT_ROOT = Path(__file__).resolve().parent
RESULTS = PROJECT_ROOT / "results"


def indexed_matches(index, query, d, strategy="spread"):
    """Generate candidates, then verify each candidate using the final DP."""
    candidates = index.generate_candidates(query, d, strategy)
    matches = {}
    for reference_id in candidates.candidate_ids:
        distance = verify(query, index.references[reference_id], d)
        if distance is not None:
            matches[reference_id] = distance
    return matches


def integer_at_least(minimum):
    def parse(value):
        parsed = int(value)
        if parsed < minimum:
            raise argparse.ArgumentTypeError(f"must be an integer >= {minimum}")
        return parsed
    return parse


def build_parser():
    parser = argparse.ArgumentParser(
        description="Search TCR sequences within an inclusive edit-distance threshold.",
        epilog="Choose a command, then specify --reference and --queries. Use search --help for details.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    summary = commands.add_parser("summarize", help="combine saved benchmark reports into one CSV")
    summary.add_argument("--input-dir", type=Path, required=True, help="directory of benchmark CSV/metadata pairs")
    summary.add_argument("--output", type=Path, required=True, help="new summary CSV; existing files are not overwritten")
    for name in ("search", "validate", "benchmark", "memory"):
        command = commands.add_parser(name)
        command.add_argument("--reference", type=Path, required=True,
                             nargs="+" if name == "benchmark" else None,
                             help="reference TSV; benchmark accepts one or more files")
        command.add_argument("--queries", type=Path, required=True, help="query TSV to read")
        command.add_argument("--d", "--distance", type=integer_at_least(0), default=1)
        command.add_argument("--k", type=integer_at_least(1), default=3)
        command.add_argument("--strategy", choices=("spread", "prefix"), default="spread")
        command.add_argument("--limit", type=integer_at_least(1))
        outputs = command.add_mutually_exclusive_group()
        outputs.add_argument("--output", type=Path, help="new output file; existing files are not overwritten")
        if name == "search":
            command.add_argument("--method", choices=("indexed", "exhaustive"), default="indexed")
            command.add_argument(
                "--annotations", type=Path,
                help="annotation JSONL or JSONL.gz; defaults to annotations.jsonl.gz or "
                     "annotations.jsonl beside the reference TSV when present",
            )
        if name == "benchmark":
            outputs.add_argument("--output-dir", type=Path,
                                 help="empty directory for benchmark reports and an automatic summary.csv")
            command.add_argument("--repeats", type=integer_at_least(2), default=5)
            command.add_argument("--metadata", type=Path)
    return parser


def prepare_outputs(paths):
    resolved = [path.resolve() for path in paths]
    if len(resolved) != len(set(resolved)):
        raise ValueError("output paths must be different")
    for path in paths:
        if path.exists() or path.is_symlink():
            raise ValueError(f"output already exists: {path}; choose a new --output path")
    for path in paths:
        path.parent.mkdir(parents=True, exist_ok=True)


def run_benchmarks(args):
    """Benchmark each reference file, then summarize a batch in one directory."""
    if len(args.reference) > 1 and args.output_dir is None:
        raise ValueError("multiple reference files require --output-dir")
    if args.output_dir is not None and args.metadata is not None:
        raise ValueError("--metadata is only used with a single --output file")
    if args.output_dir is not None and args.output_dir.exists() and any(args.output_dir.iterdir()):
        raise ValueError("output directory must be empty; choose a new --output-dir")

    # Load the same queries once for all reference databases.
    queries = read_sequences(args.queries, "query_id")
    if args.limit is not None:
        queries = dict(list(queries.items())[:args.limit])
    databases = [read_sequences(path, "reference_id") for path in args.reference]
    if not queries or any(not references for references in databases):
        raise ValueError("reference and query files must not be empty")
    if args.output_dir is not None and len({len(references) for references in databases}) > 1:
        raise ValueError("reference files must have the same size for a batch summary")

    # Check every output path before starting the benchmarks.
    if args.output_dir is not None:
        output_files = [args.output_dir / (path.stem.removeprefix("reference_") + ".csv")
                        for path in args.reference]
    else:
        output_files = [args.output if args.output is not None else RESULTS / "benchmark.csv"]
    metadata_files = [args.metadata if args.metadata is not None else
                      path.with_name(path.stem + ".metadata.json") for path in output_files]
    summary_file = args.output_dir / "summary.csv" if args.output_dir is not None else None
    prepare_outputs(output_files + metadata_files + ([summary_file] if summary_file else []))

    for reference_path, references, output, metadata_path in zip(
            args.reference, databases, output_files, metadata_files):
        metadata, rows = benchmark(references, queries, args.d, args.k, args.strategy,
                                   args.repeats, indexed_matches)
        metadata.update(reference_path=str(reference_path.resolve()), query_path=str(args.queries.resolve()))
        write_csv(output, rows)
        with metadata_path.open("x", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        print(f"Benchmark results: {output}; metadata: {metadata_path}")

    if summary_file is not None:
        write_csv(summary_file, summarize_benchmarks(args.output_dir))
        print(f"Summary: {summary_file}")
    return 0


def run_search(args, references, query_records, annotations):
    """Build the index once, search every query, and output matching records."""
    index = InvertedIndex(references, args.k) if args.method == "indexed" else None
    output = sys.stdout
    if args.output is not None:
        prepare_outputs([args.output])
        output = args.output.open("x", encoding="utf-8", newline="\n")
    total_matches = 0
    try:
        for query_id, record in query_records.items():
            query = record["sequence"]
            matches = (indexed_matches(index, query, args.d, args.strategy) if index is not None
                       else exhaustive_matches(query, references, args.d))
            row = {
                "query_id": query_id,
                "query_sequence": query,
                "query_metadata": {key: value for key, value in record.items()
                                   if key not in ("query_id", "sequence")},
                "matches": [
                    {"reference_id": rid, "sequence": references[rid],
                     "distance": matches[rid], "annotations": annotations.get(rid, [])}
                    for rid in sorted(matches)
                ],
            }
            total_matches += len(matches)
            print(json.dumps(row, ensure_ascii=False), file=output)
    finally:
        if output is not sys.stdout:
            output.close()
    print(f"method={args.method}; references={len(references)}; "
          f"queries={len(query_records)}; matches={total_matches}", file=sys.stderr)
    return 0


def main(argv=None):
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        parser.print_help()
        return 0
    args = parser.parse_args(argv)
    try:
        if args.command == "summarize":
            rows = summarize_benchmarks(args.input_dir)
            prepare_outputs([args.output])
            write_csv(args.output, rows)
            print(f"Summarized {len(rows)} benchmark runs: {args.output}")
            return 0
        if args.command == "benchmark":
            return run_benchmarks(args)
        references = read_sequences(args.reference, "reference_id")
        query_records = read_records(args.queries, "query_id")
        if args.limit is not None:
            query_records = dict(list(query_records.items())[:args.limit])
        if not references or not query_records:
            raise ValueError("reference and query files must not be empty")
        if args.command == "search":
            annotations_path = args.annotations
            if annotations_path is None:
                annotations_path = next(
                    (args.reference.parent / name
                     for name in ("annotations.jsonl.gz", "annotations.jsonl")
                     if (args.reference.parent / name).is_file()), None,
                )
            annotations = (read_annotations(annotations_path, references)
                           if annotations_path is not None else {})
            return run_search(args, references, query_records, annotations)

        queries = {qid: record["sequence"] for qid, record in query_records.items()}

        output = args.output if args.output is not None else RESULTS / f"{args.command}.csv"
        if args.command == "validate":
            prepare_outputs([output])
            rows = validation_rows(references, queries, InvertedIndex(references, args.k),
                                   args.d, args.strategy, indexed_matches)
            write_csv(output, rows)
            failures = sum(not row["match"] for row in rows)
            print(f"Checked {len(rows)} comparisons; mismatches: {failures}. Details: {output}")
            return 1 if failures else 0

        if args.command == "memory":
            prepare_outputs([output])
            rows = memory_comparison(references, queries, args.d, args.k, args.strategy, indexed_matches)
            write_csv(output, rows)
            print(f"Memory results: {output} (incremental Python allocations; not process RSS)")
            return 0

    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")
    except KeyboardInterrupt:
        print("\nCancelled.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
