"""Generate synthetic references, edited queries and common GQE queries."""

import argparse
import json
from pathlib import Path
import random

from draft.Verification.verification_basic import verification_basic
from .common import (
    AMINO_ACIDS, PYTHON_VERSION, QUERY_FIELDS, REFERENCE_FIELDS,
    query_row, sequence_id, write_json, write_tsv,
)

MOTIF = "GQE"


def make_pool(length, size, motif, present, rng):
    # GENERATE random sequences -> CHECK motif status -> KEEP unique sequences.
    sequences = []
    seen = set()
    motif_start = (length - len(motif)) // 2
    for attempt in range(max(10000, size * 1000)):
        letters = ["C"]
        for position in range(length - 2):
            letters.append(rng.choice(AMINO_ACIDS))
        letters.append("F")
        sequence = "".join(letters)

        if present:
            sequence = sequence[:motif_start] + motif + sequence[motif_start + len(motif):]
        if (motif in sequence) != present:
            continue
        if sequence in seen:
            continue
        seen.add(sequence)
        sequences.append(sequence)
        if len(sequences) == size:
            return sequences
    raise ValueError("Too few distinct sequences for these size, length and motif settings")


def edit_sequence(source, operation, count, rng):
    # APPLY edits -> CHECK distance to the source -> RETRY if needed.
    for attempt in range(100):
        sequence = source
        events = []
        positions = range(count)
        if operation == "substitution":
            positions = rng.sample(range(1, len(source) - 1), count)
        for position in positions:
            if operation == "insertion":
                position = rng.randrange(1, len(sequence))
                old = ""
                new = rng.choice(AMINO_ACIDS)
            elif operation == "deletion":
                position = rng.randrange(1, len(sequence) - 1)
                old = sequence[position]
                new = ""
            else:
                old = sequence[position]
                choices = AMINO_ACIDS.replace(old, "")
                new = rng.choice(choices)
            sequence = sequence[:position] + new + sequence[position + len(old):]
            events.append({"operation": operation, "position": position, "from": old, "to": new})

        # Borrow the group's basic DP; this checks one pair, not the whole library.
        distance = verification_basic(source, sequence)
        if distance == count:
            return sequence, events, distance
    raise ValueError("Could not produce a query at the requested source distance")


def write_queries(directory, references, sources, args, length, name):
    # PLAN edits -> CREATE queries -> ADD short and candidate-negative queries -> SAVE.
    queries = []
    histories = []
    rng = random.Random(args.seed + length * 10 + 3)
    operations = [("exact_copy", 0)]
    for operation in ("substitution", "insertion", "deletion"):
        for count in range(1, args.max_edits + 1):
            operations.append((operation, count))

    for source in sources:
        for operation, count in operations:
            sequence, events, distance = edit_sequence(source, operation, count, rng)
            query_id = f"{name}_Q{len(queries) + 1:04d}"
            row = query_row(query_id, sequence, operation, source_id=references[source],
                            source=source, edits=len(events), distance=distance)
            queries.append(row)
            histories.append({"query_id": query_id, "events": events})

    for length_short in (1, 2, 3):
        for letter in "AW":
            query_id = f"{name}_Q{len(queries) + 1:04d}"
            queries.append(query_row(query_id, letter * length_short, "short_query"))
    for number in range(args.anchors):
        letters = []
        for position in range(length):
            letters.append(rng.choice("WY"))
        query_id = f"{name}_Q{len(queries) + 1:04d}"
        queries.append(query_row(query_id, "".join(letters), "candidate_no_match"))

    write_tsv(directory / "queries.tsv", QUERY_FIELDS, queries)
    with (directory / "edit_history.jsonl").open("w", encoding="utf-8") as handle:
        for history in histories:
            handle.write(json.dumps(history) + "\n")
    categories = {}
    for row in queries:
        category = row["category"]
        if category not in categories:
            categories[category] = 0
        categories[category] += 1
    return categories


def make_probes(output, length, sources, existing, args):
    # CHANGE one residue outside GQE -> REJECT exact copies -> SAVE common queries.
    if args.probe_count == 0:
        return
    motif_start = (length - len(MOTIF)) // 2
    templates = []
    for source in sources:
        if source.count(MOTIF) == 1:
            templates.append(source)
    templates.sort()
    rng = random.Random(args.probe_seed + length)
    rng.shuffle(templates)
    queries = []

    for source in templates:
        options = []
        for position in range(1, length - 1):
            if motif_start <= position < motif_start + len(MOTIF):
                continue
            for letter in AMINO_ACIDS:
                if letter == source[position]:
                    continue
                query = source[:position] + letter + source[position + 1:]
                if query in existing or query.count(MOTIF) != 1:
                    continue
                options.append(query)
        if not options:
            continue
        query = rng.choice(options)
        existing.add(query)
        query_id = f"GQE_L{length:02d}_Q{len(queries) + 1:02d}"
        row = query_row(query_id, query, "motif_probe", source_id=sequence_id("SYN_", source),
                        source=source, edits=1, distance=1)
        queries.append(row)
        if len(queries) == args.probe_count:
            write_tsv(output / "motif_probes" / f"L{length}" / "queries.tsv", QUERY_FIELDS, queries)
            return
    raise ValueError("Too few shared motif templates. Reduce --probe-count or increase --size")


def prepare(args, output):
    # CREATE two pools -> MIX each proportion -> SAVE libraries -> MAKE shared GQE queries.
    reports = []
    motif_counts = []
    for percentage in args.motif_percentages:
        motif_counts.append(round(args.size * percentage / 100))

    for length in args.lengths:
        absent_rng = random.Random(args.seed + length * 10)
        present_rng = random.Random(args.seed + length * 10 + 1)
        absent = make_pool(length, args.size, MOTIF, False, absent_rng)
        present = make_pool(length, args.size, MOTIF, True, present_rng)

        for percentage in args.motif_percentages:
            motif_count = round(args.size * percentage / 100)
            sequences = absent[:args.size - motif_count] + present[:motif_count]
            shuffle_rng = random.Random(args.seed + length * 10 + 2)
            shuffle_rng.shuffle(sequences)
            references = {}
            reference_rows = []
            used_ids = set()
            for sequence in sequences:
                reference_id = sequence_id("SYN_", sequence)
                if reference_id in used_ids:
                    raise ValueError("Duplicate reference ID or sequence")
                used_ids.add(reference_id)
                references[sequence] = reference_id
                reference_rows.append({"reference_id": reference_id, "sequence": sequence})

            name = f"L{length}_motif{percentage:03d}"
            directory = output / name
            write_tsv(directory / "reference.tsv", REFERENCE_FIELDS, reference_rows)
            categories = write_queries(directory, references, absent[:args.anchors], args, length, name)
            actual_motif_count = 0
            for sequence in sequences:
                if MOTIF in sequence:
                    actual_motif_count += 1
            report = {
                "name": name, "reference_count": len(references), "reference_length": length,
                "motif": MOTIF, "motif_count": actual_motif_count,
                "query_categories": categories, "seed": args.seed,
            }
            write_json(directory / "generation_report.json", report)
            reports.append(report)

        # All positive-proportion libraries contain this prefix of the same pool.
        positive_counts = []
        for count in motif_counts:
            if count > 0:
                positive_counts.append(count)
        shared_count = 0
        if positive_counts:
            shared_count = min(positive_counts)
        existing = set(absent[:args.size - min(motif_counts)])
        existing.update(present[:max(motif_counts)])
        make_probes(output, length, present[:shared_count], existing, args)

    write_json(output / "generation_report.json", {
        "size": args.size, "lengths": args.lengths, "motif": MOTIF,
        "motif_percentages": args.motif_percentages, "max_edits": args.max_edits,
        "anchors": args.anchors, "seed": args.seed, "probe_seed": args.probe_seed,
        "probe_count_per_length": args.probe_count, "python_version": PYTHON_VERSION,
        "simulation": "Amino acid stress tests with C/F anchors and interior edits",
        "negative_query_status": "candidate_no_match queries await exhaustive validation",
        "source_distance": "Checked using verification_basic for generated edits",
        "datasets": reports,
    })
    print(f"Created {len(reports)} synthetic libraries")


def main():
    # READ settings -> CHECK settings -> PREPARE files in a new directory.
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

    if args.size < 1 or not 1 <= args.anchors <= args.size:
        parser.error("Size and anchors must be positive, and anchors must fit in the library")
    if args.probe_count < 0:
        parser.error("Probe count must be nonnegative")
    if min(args.lengths) < 6:
        parser.error("Lengths must be at least 6 to fit GQE and the C/F ends")
    if not 0 <= args.max_edits <= min(args.lengths) - 2:
        parser.error("Requested edits must fit the interior of the shortest sequence")
    # Each free position has 20 choices: do not request more than the pool can hold.
    if args.size > 20 ** (min(args.lengths) - len(MOTIF) - 2):
        parser.error("Size exceeds the possible distinct motif sequences at the shortest length")
    has_positive_library = False
    for percentage in args.motif_percentages:
        if not 0 <= percentage <= 100:
            parser.error("Percentages must be between 0 and 100")
        motif_count = round(args.size * percentage / 100)
        if args.size - motif_count < args.anchors:
            parser.error("Each library must leave enough background templates")
        if motif_count > 0:
            has_positive_library = True
    if args.probe_count > 0 and not has_positive_library:
        parser.error("Motif probes require a library with positive motif percentage")
    if args.out.exists() or args.out.is_symlink():
        parser.error("Output directory exists. Choose a new --out directory")
    try:
        prepare(args, args.out)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
