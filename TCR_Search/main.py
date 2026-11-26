"""Run indexed/exhaustive search and compare their performance."""

import argparse
import json
import sys
from pathlib import Path

# Keep the project folders free of automatically generated __pycache__ files.
sys.dont_write_bytecode = True

from tcr_search.benchmark import benchmark, memory_comparison, validation_rows, write_csv
from tcr_search.data_input import read_sequences
from tcr_search.exhaustive_search import exhaustive_matches
from tcr_search.inverted_index import InvertedIndex
from tcr_search.verification import verify

PROJECT_ROOT = Path(__file__).resolve().parent
EXAMPLES = PROJECT_ROOT / "data" / "examples"
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
        epilog="With no arguments, search the examples with d=1 and k=3.",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("search", "validate", "benchmark", "memory"):
        command = commands.add_parser(name)
        command.add_argument("--reference", type=Path, default=EXAMPLES / "reference.tsv")
        command.add_argument("--queries", type=Path, default=EXAMPLES / "queries.tsv")
        command.add_argument("--d", "--distance", type=integer_at_least(0), default=1)
        command.add_argument("--k", type=integer_at_least(1), default=3)
        command.add_argument("--strategy", choices=("spread", "prefix"), default="spread")
        command.add_argument("--limit", type=integer_at_least(1))
        command.add_argument("--output", type=Path, help="new output file; existing files are not overwritten")
        if name == "search":
            command.add_argument("--method", choices=("indexed", "exhaustive"), default="indexed")
        if name == "benchmark":
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


def run_search(args, references, queries):
    """Build the index once, search every query, and output matching records."""
    index = InvertedIndex(references, args.k) if args.method == "indexed" else None
    output = sys.stdout
    if args.output is not None:
        prepare_outputs([args.output])
        output = args.output.open("x", encoding="utf-8", newline="\n")
    total_matches = 0
    try:
        for query_id, query in queries.items():
            matches = (indexed_matches(index, query, args.d, args.strategy) if index is not None
                       else exhaustive_matches(query, references, args.d))
            row = {
                "query_id": query_id,
                "matches": [
                    {"reference_id": rid, "sequence": references[rid],
                     "distance": matches[rid], "annotations": []}
                    for rid in sorted(matches)
                ],
            }
            total_matches += len(matches)
            print(json.dumps(row, ensure_ascii=False), file=output)
    finally:
        if output is not sys.stdout:
            output.close()
    print(f"method={args.method}; references={len(references)}; "
          f"queries={len(queries)}; matches={total_matches}", file=sys.stderr)
    return 0


def main(argv=None):
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(argv or ["search"])
    try:
        references = read_sequences(args.reference, "reference_id")
        queries = read_sequences(args.queries, "query_id")
        if args.limit is not None:
            queries = dict(list(queries.items())[:args.limit])
        if not references or not queries:
            raise ValueError("reference and query files must not be empty")
        if args.command == "search":
            return run_search(args, references, queries)

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

        metadata_path = args.metadata if args.metadata is not None else output.with_name(
            output.stem + ".metadata.json")
        prepare_outputs([output, metadata_path])
        metadata, rows = benchmark(references, queries, args.d, args.k, args.strategy,
                                   args.repeats, indexed_matches)
        metadata.update(reference_path=str(args.reference.resolve()), query_path=str(args.queries.resolve()))
        write_csv(output, rows)
        with metadata_path.open("x", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        print(f"Benchmark results: {output}; metadata: {metadata_path}")
        return 0
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")
    except KeyboardInterrupt:
        print("\nCancelled.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    raise SystemExit(main())
