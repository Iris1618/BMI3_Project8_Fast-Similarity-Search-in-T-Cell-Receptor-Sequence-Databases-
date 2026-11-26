"""Small helpers for sequence checks, IDs and output files."""

import csv
import hashlib
import json
import sys

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
REFERENCE_FIELDS = ["reference_id", "sequence"]
QUERY_FIELDS = [
    "query_id", "sequence", "category", "source_id", "source_sequence",
    "applied_edits", "source_distance",
]
PYTHON_VERSION = sys.version.split()[0]


def valid_sequence(sequence):
    if sequence == "":
        return False
    for letter in sequence:
        if letter not in AMINO_ACIDS:
            return False
    return True


def sequence_id(prefix, sequence):
    # Keep the same IDs as the existing data package.
    encoded_sequence = sequence.encode("ascii")
    code = hashlib.sha256(encoded_sequence).hexdigest()
    return prefix + code[:16]


def file_checksum(path):
    # Record which exact source file was used.
    checksum = hashlib.sha256()
    with path.open("rb") as source:
        for line in source:
            checksum.update(line)
    return checksum.hexdigest()


def query_row(query_id, sequence, category, source_id="", source="", edits="", distance=""):
    return {
        "query_id": query_id, "sequence": sequence, "category": category,
        "source_id": source_id, "source_sequence": source,
        "applied_edits": edits, "source_distance": distance,
    }


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
