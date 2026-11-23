"""Small file helpers shared by the two input preparation scripts."""

import csv
import hashlib
import json

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
REFERENCE_FIELDS = ["reference_id", "sequence"]
QUERY_FIELDS = [
    "query_id", "sequence", "category", "source_id", "source_sequence",
    "applied_edits", "source_distance",
]


def sequence_id(prefix, sequence):
    return prefix + hashlib.sha256(sequence.encode("ascii")).hexdigest()[:16]


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
        writer.writerows(rows)


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
