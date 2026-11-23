"""Generate synthetic references, edited queries and common motif probes."""

import argparse
from collections import Counter
import json
from pathlib import Path
import platform
import random
import tempfile

from verification_basic import verification_basic
from .common import (
    AMINO_ACIDS, QUERY_FIELDS, REFERENCE_FIELDS, query_row,
    sequence_id, write_json, write_tsv,
)

MOTIF = "GQE"


def make_pool(length, size, motif, present, rng):
    sequences = []
    seen = set()
    start = (length - len(motif)) // 2
    for _ in range(max(10000, size * 1000)):
        sequence = "C" + "".join(rng.choice(AMINO_ACIDS) for _ in range(length - 2)) + "F"
        if present:
            sequence = sequence[:start] + motif + sequence[start + len(motif):]
        if (motif in sequence) == present and sequence not in seen:
            seen.add(sequence)
            sequences.append(sequence)
        if len(sequences) == size:
            return sequences
    raise ValueError("Too few distinct sequences for these size, length and motif settings")


def edit_sequence(source, operation, count, rng):
    for _ in range(100):
        sequence = source
        events = []
        if operation == "substitution":
            positions = rng.sample(range(1, len(source) - 1), count)
        else:
            positions = [None] * count
        for position in positions:
            if operation == "insertion":
                position = rng.randrange(1, len(sequence))
                old, new = "", rng.choice(AMINO_ACIDS)
            elif operation == "deletion":
                position = rng.randrange(1, len(sequence) - 1)
                old, new = sequence[position], ""
            else:
                old = sequence[position]
                new = rng.choice(AMINO_ACIDS.replace(old, ""))
            events.append({"operation": operation, "position": position, "from": old, "to": new})
            sequence = sequence[:position] + new + sequence[position + len(old):]
        # Reuse the group's verifier for this source-query pair.
        distance = verification_basic(source, sequence)
        if distance == count:
            return sequence, events, distance
    raise ValueError("Could not produce a query at the requested source distance")


def write_queries(directory, references, sources, args, length, name):
    queries, histories = [], []
    rng = random.Random(args.seed + length * 10 + 3)
    for source in sources:
        examples = [(source, "exact_copy", [], 0)]
        for operation in ("substitution", "insertion", "deletion"):
            for count in range(1, args.max_edits + 1):
                sequence, events, distance = edit_sequence(source, operation, count, rng)
                examples.append((sequence, operation, events, distance))
        for sequence, category, events, distance in examples:
            query_id = f"{name}_Q{len(queries) + 1:04d}"
            queries.append(query_row(query_id, sequence, category,
                                     references[source], source, len(events), distance))
            histories.append({"query_id": query_id, "events": events})
    for n in (1, 2, 3):
        for char in "AW":
            queries.append(query_row(f"{name}_Q{len(queries) + 1:04d}", char * n, "short_query"))
    for _ in range(args.anchors):
        sequence = "".join(rng.choice("WY") for _ in range(length))
        queries.append(query_row(f"{name}_Q{len(queries) + 1:04d}", sequence, "candidate_no_match"))
    write_tsv(directory / "queries.tsv", QUERY_FIELDS, queries)
    with (directory / "edit_history.jsonl").open("w", encoding="utf-8") as handle:
        for history in histories:
            handle.write(json.dumps(history) + "\n")
    return dict(Counter(row["category"] for row in queries))


def make_probes(output, length, libraries, args):
    if args.probe_count == 0:
        return
    motif_start = (length - len(MOTIF)) // 2
    positive = [refs for refs in libraries if any(MOTIF in s for s in refs)]
    if not positive:
        raise ValueError("Motif probes require a library with positive motif percentage")
    shared_sources = set(positive[0])
    for refs in positive[1:]:
        shared_sources.intersection_update(refs)
    sources = sorted(
        s for s in shared_sources
        if s.count(MOTIF) == 1 and s[motif_start:motif_start + len(MOTIF)] == MOTIF
    )
    rng = random.Random(args.probe_seed + length)
    rng.shuffle(sources)
    existing = set().union(*(set(refs) for refs in libraries))
    queries = []
    for source in sources:
        options = []
        for position in range(1, length - 1):
            if motif_start <= position < motif_start + len(MOTIF):
                continue
            for letter in AMINO_ACIDS:
                query = source[:position] + letter + source[position + 1:]
                if letter != source[position] and query not in existing and query.count(MOTIF) == 1:
                    options.append(query)
        if not options:
            continue
        query = rng.choice(options)
        existing.add(query)
        query_id = f"GQE_L{length:02d}_Q{len(queries) + 1:02d}"
        queries.append(query_row(query_id, query, "motif_probe", positive[0][source], source, 1, 1))
        if len(queries) == args.probe_count:
            write_tsv(output / "motif_probes" / f"L{length}" / "queries.tsv", QUERY_FIELDS, queries)
            return
    raise ValueError("Too few shared motif templates. Reduce --probe-count or increase --size")


def prepare(args, output):
    reports = []
    for length in args.lengths:
        absent = make_pool(length, args.size, MOTIF, False, random.Random(args.seed + length * 10))
        present = make_pool(length, args.size, MOTIF, True, random.Random(args.seed + length * 10 + 1))
        libraries = []
        for percentage in args.motif_percentages:
            motif_count = round(args.size * percentage / 100)
            sequences = absent[:args.size - motif_count] + present[:motif_count]
            random.Random(args.seed + length * 10 + 2).shuffle(sequences)
            references = {s: sequence_id("SYN_", s) for s in sequences}
            if len(set(references.values())) != args.size:
                raise ValueError("Duplicate reference ID or sequence")
            libraries.append(references)
            name = f"L{length}_motif{percentage:03d}"
            directory = output / name
            write_tsv(directory / "reference.tsv", REFERENCE_FIELDS,
                      ({"reference_id": rid, "sequence": s} for s, rid in references.items()))
            categories = write_queries(directory, references, absent[:args.anchors], args, length, name)
            report = {"name": name, "reference_count": len(references), "reference_length": length,
                      "motif": MOTIF, "motif_count": sum(MOTIF in s for s in references),
                      "query_categories": categories, "seed": args.seed}
            write_json(directory / "generation_report.json", report)
            reports.append(report)
        make_probes(output, length, libraries, args)
    write_json(output / "generation_report.json", {
        "size": args.size, "lengths": args.lengths, "motif": MOTIF,
        "motif_percentages": args.motif_percentages, "max_edits": args.max_edits,
        "anchors": args.anchors, "seed": args.seed, "probe_seed": args.probe_seed,
        "probe_count_per_length": args.probe_count, "python_version": platform.python_version(),
        "simulation": "Amino acid stress tests with C/F anchors and interior edits",
        "negative_query_status": "candidate_no_match queries await exhaustive validation",
        "source_distance": "Checked using verification_basic for generated edits",
        "datasets": reports,
    })
    print(f"Created {len(reports)} synthetic libraries")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--size", type=int, default=1000)
    parser.add_argument("--lengths", nargs="+", type=int, default=[8, 12, 18])
    parser.add_argument("--motif-percentages", nargs="+", type=int, default=[0, 50, 90])
    parser.add_argument("--anchors", type=int, default=5)
    parser.add_argument("--max-edits", type=int, default=3)
    parser.add_argument("--probe-count", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20262111)
    parser.add_argument("--probe-seed", type=int, default=20261113)
    args = parser.parse_args()
    args.lengths = sorted(set(args.lengths))
    args.motif_percentages = sorted(set(args.motif_percentages))
    if args.size < 1 or not 1 <= args.anchors <= args.size or args.probe_count < 0:
        parser.error("Size and anchors must be positive. Probe count must be nonnegative")
    if min(args.lengths) < 6 or not 0 <= args.max_edits <= min(args.lengths) - 2:
        parser.error("Lengths must fit the motif, anchors and requested edits")
    if args.size > 20 ** (min(args.lengths) - len(MOTIF) - 2):
        parser.error("Size exceeds the possible distinct motif sequences at the shortest length")
    for percentage in args.motif_percentages:
        if not 0 <= percentage <= 100 or args.size - round(args.size * percentage / 100) < args.anchors:
            parser.error("Percentages must be 0-100 and leave enough background source sequences")
    if args.out.exists():
        parser.error("Output directory exists. Choose a new --out directory")
    try:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=args.out.parent) as temporary:
            output = Path(temporary) / "synthetic"
            output.mkdir()
            prepare(args, output)
            output.rename(args.out)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
