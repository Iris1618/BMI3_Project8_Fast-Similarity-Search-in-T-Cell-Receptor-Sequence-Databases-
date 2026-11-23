"""Prepare reference sequences, queries and annotations from VDJdb."""

import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import platform
import random
import tempfile

from .common import (
    AMINO_ACIDS, QUERY_FIELDS, REFERENCE_FIELDS, query_row,
    sequence_id, write_json, write_tsv,
)


def prepare(args, output):
    # READ records -> FILTER species/chain -> CLEAN CDR3 -> GROUP duplicates.
    records = {}
    counts = Counter(raw_rows=0, target_rows=0, valid_rows=0, missing_cdr3=0,
                     invalid_amino_acid=0, normalised_rows=0, excluded_species_or_chain=0)
    excluded = []
    digest = hashlib.sha256()
    with args.input.open("rb") as source:
        header = source.readline()
        digest.update(header)
        columns = header.decode("utf-8-sig").rstrip("\r\n").split("\t")
        if len(set(columns)) != len(columns) or not {"gene", "species", "cdr3"} <= set(columns):
            raise ValueError("Input requires unique columns including gene, species and cdr3")
        for line_number, line in enumerate(source, 2):
            digest.update(line)
            # VDJdb has literal TSV fields, including quotes and long annotations.
            values = line.decode("utf-8").rstrip("\r\n").split("\t")
            if len(values) != len(columns):
                raise ValueError(f"Malformed TSV row at line {line_number}")
            row = dict(zip(columns, values))
            counts["raw_rows"] += 1
            if row["species"] != args.species or row["gene"] != args.gene:
                counts["excluded_species_or_chain"] += 1
                continue
            counts["target_rows"] += 1
            sequence = row["cdr3"].strip().upper()
            if not sequence or set(sequence) - set(AMINO_ACIDS):
                reason = "invalid_amino_acid" if sequence else "missing_cdr3"
                counts[reason] += 1
                excluded.append({"source_line": line_number, "cdr3": row["cdr3"], "reason": reason})
                continue
            counts["valid_rows"] += 1
            counts["normalised_rows"] += sequence != row["cdr3"]
            records.setdefault(sequence, []).append({"source_line": line_number, "fields": row})

    # SAVE each unique sequence together with all its annotation records.
    sequences = sorted(records)
    if max(args.sizes) > len(sequences):
        raise ValueError(f"Requested subset exceeds {len(sequences)} available unique sequences")
    ids = {sequence: sequence_id("VDJ_", sequence) for sequence in sequences}
    if len(set(ids.values())) != len(sequences):
        raise ValueError("Sequence ID collision")
    write_tsv(output / "reference_all.tsv", REFERENCE_FIELDS,
              ({"reference_id": ids[s], "sequence": s} for s in sequences))
    with gzip.open(output / "annotations.jsonl.gz", "wt", encoding="utf-8") as annotations:
        for sequence in sequences:
            row = {"reference_id": ids[sequence], "records": records[sequence]}
            annotations.write(json.dumps(row, ensure_ascii=False) + "\n")
    write_tsv(output / "excluded_target_rows.tsv", ["source_line", "cdr3", "reason"], excluded)
    lengths = Counter(map(len, sequences))
    write_tsv(output / "length_distribution.tsv", ["length", "unique_sequences"],
              ({"length": n, "unique_sequences": count} for n, count in sorted(lengths.items())))
    write_tsv(output / "length_review_flags.tsv", ["reference_id", "sequence", "length"],
              ({"reference_id": ids[s], "sequence": s, "length": len(s)}
               for s in sequences if len(s) < 8 or len(s) > 30))

    # SHUFFLE once -> TAKE nested subsets -> SAMPLE queries from the smallest.
    random.Random(args.seed).shuffle(sequences)
    for size in args.sizes:
        write_tsv(output / f"reference_n{size}.tsv", REFERENCE_FIELDS,
                  ({"reference_id": ids[s], "sequence": s} for s in sequences[:size]))
    queries = random.Random(args.seed + 1).sample(sequences[:min(args.sizes)], args.queries)
    write_tsv(output / "queries.tsv", QUERY_FIELDS, (
        query_row(f"REAL_Q{i:04d}", s, "exact_copy", ids[s], s, 0, 0)
        for i, s in enumerate(queries, 1)
    ))
    counts["unique_sequences"] = len(records)
    counts["duplicate_rows_merged"] = counts["valid_rows"] - len(records)
    write_json(output / "cleaning_report.json", {
        "input_filename": args.input.name, "input_sha256": digest.hexdigest(),
        "species": args.species, "gene": args.gene, "alphabet": AMINO_ACIDS,
        "reference_sizes": args.sizes, "query_count": args.queries,
        "reference_seed": args.seed, "query_seed": args.seed + 1,
        "python_version": platform.python_version(), "counts": dict(counts),
        "length_review": "Lengths outside 8-30 aa are flagged and retained",
    })
    print(f"Kept {counts['valid_rows']} records and {len(records)} unique sequences")


def main():
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
    if args.out.exists():
        parser.error("Output directory exists. Choose a new --out directory")
    try:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=args.out.parent) as temporary:
            prepare(args, Path(temporary))
            Path(temporary).rename(args.out)
    except (OSError, ValueError) as error:
        parser.exit(2, f"error: {error}\n")


if __name__ == "__main__":
    main()
