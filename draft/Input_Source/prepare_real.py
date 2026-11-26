"""Prepare reference sequences, queries and annotations from VDJdb."""

import argparse
import gzip
import json
from pathlib import Path
import random

from .common import (
    AMINO_ACIDS, PYTHON_VERSION, QUERY_FIELDS, REFERENCE_FIELDS,
    file_checksum, query_row, sequence_id, valid_sequence, write_json, write_tsv,
)


def prepare(args, output):
    # READ records -> SELECT species and chain -> CLEAN CDR3 -> GROUP duplicates.
    records = {}
    counts = {
        "raw_rows": 0, "target_rows": 0, "valid_rows": 0,
        "missing_cdr3": 0, "invalid_amino_acid": 0,
        "normalised_rows": 0, "excluded_species_or_chain": 0,
    }
    excluded = []
    with args.input.open(encoding="utf-8-sig", newline="") as source:
        columns = source.readline().rstrip("\r\n").split("\t")
        if len(set(columns)) != len(columns):
            raise ValueError("Input column names must be unique")
        for field in ("gene", "species", "cdr3"):
            if field not in columns:
                raise ValueError("Input requires gene, species and cdr3 columns")

        for line_number, line in enumerate(source, start=2):
            # Split literal tabs: VDJdb annotation fields can contain quotes.
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(columns):
                raise ValueError(f"Malformed TSV row at line {line_number}")
            row = {}
            for i in range(len(columns)):
                row[columns[i]] = values[i]
            counts["raw_rows"] += 1

            if row["species"] != args.species or row["gene"] != args.gene:
                counts["excluded_species_or_chain"] += 1
                continue
            counts["target_rows"] += 1
            sequence = row["cdr3"].strip().upper()
            if not valid_sequence(sequence):
                if sequence == "":
                    reason = "missing_cdr3"
                else:
                    reason = "invalid_amino_acid"
                counts[reason] += 1
                excluded.append({"source_line": line_number, "cdr3": row["cdr3"], "reason": reason})
                continue

            counts["valid_rows"] += 1
            if sequence != row["cdr3"]:
                counts["normalised_rows"] += 1
            if sequence not in records:
                records[sequence] = []
            records[sequence].append({"source_line": line_number, "fields": row})

    # ASSIGN IDs -> SAVE unique references and every original annotation record.
    sequences = sorted(records)
    if max(args.sizes) > len(sequences):
        raise ValueError(f"Requested subset exceeds {len(sequences)} available unique sequences")
    ids = {}
    reference_rows = []
    used_ids = set()
    for sequence in sequences:
        reference_id = sequence_id("VDJ_", sequence)
        if reference_id in used_ids:
            raise ValueError("Sequence ID collision")
        used_ids.add(reference_id)
        ids[sequence] = reference_id
        reference_rows.append({"reference_id": reference_id, "sequence": sequence})
    write_tsv(output / "reference_all.tsv", REFERENCE_FIELDS, reference_rows)

    # gzip keeps the annotation file small and compatible with the data package.
    with gzip.open(output / "annotations.jsonl.gz", "wt", encoding="utf-8") as handle:
        for sequence in sequences:
            row = {"reference_id": ids[sequence], "records": records[sequence]}
            handle.write(json.dumps(row, ensure_ascii=False) + "\n")
    write_tsv(output / "excluded_target_rows.tsv", ["source_line", "cdr3", "reason"], excluded)

    # COUNT lengths and flag unusual lengths without removing those sequences.
    length_counts = {}
    flagged_rows = []
    for sequence in sequences:
        length = len(sequence)
        if length not in length_counts:
            length_counts[length] = 0
        length_counts[length] += 1
        if length < 8 or length > 30:
            flagged_rows.append({"reference_id": ids[sequence], "sequence": sequence, "length": length})
    length_rows = []
    for length in sorted(length_counts):
        length_rows.append({"length": length, "unique_sequences": length_counts[length]})
    write_tsv(output / "length_distribution.tsv", ["length", "unique_sequences"], length_rows)
    write_tsv(output / "length_review_flags.tsv", ["reference_id", "sequence", "length"], flagged_rows)

    # SHUFFLE once -> TAKE nested subsets -> SAMPLE queries from the smallest.
    reference_rng = random.Random(args.seed)
    reference_rng.shuffle(sequences)
    for size in args.sizes:
        subset_rows = []
        for sequence in sequences[:size]:
            subset_rows.append({"reference_id": ids[sequence], "sequence": sequence})
        write_tsv(output / f"reference_n{size}.tsv", REFERENCE_FIELDS, subset_rows)

    query_rng = random.Random(args.seed + 1)
    selected = query_rng.sample(sequences[:min(args.sizes)], args.queries)
    query_rows = []
    for number, sequence in enumerate(selected, start=1):
        row = query_row(f"REAL_Q{number:04d}", sequence, "exact_copy",
                        source_id=ids[sequence], source=sequence, edits=0, distance=0)
        query_rows.append(row)
    write_tsv(output / "queries.tsv", QUERY_FIELDS, query_rows)

    counts["unique_sequences"] = len(records)
    counts["duplicate_rows_merged"] = counts["valid_rows"] - len(records)
    write_json(output / "cleaning_report.json", {
        "input_filename": args.input.name, "input_sha256": file_checksum(args.input),
        "species": args.species, "gene": args.gene, "alphabet": AMINO_ACIDS,
        "reference_sizes": args.sizes, "query_count": args.queries,
        "reference_seed": args.seed, "query_seed": args.seed + 1,
        "python_version": PYTHON_VERSION, "counts": counts,
        "length_review": "Lengths outside 8-30 aa are flagged and retained",
    })
    print(f"Kept {counts['valid_rows']} records and {len(records)} unique sequences")


def main():
    # READ settings -> CHECK settings -> PREPARE files in a new directory.
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Path to vdjdb.slim.txt")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--sizes", nargs="+", type=int, default=[1000, 5000, 10000])
    parser.add_argument("--queries", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--species", default="HomoSapiens")
    parser.add_argument("--gene", default="TRB")
    args = parser.parse_args()
    args.sizes = sorted(set(args.sizes))
    if min(args.sizes) < 1 or not 1 <= args.queries <= min(args.sizes):
        parser.error("Sizes must be positive and queries must fit in the smallest subset")
    if args.out.exists() or args.out.is_symlink():
        parser.error("Output directory exists. Choose a new --out directory")
    try:
        prepare(args, args.out)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
