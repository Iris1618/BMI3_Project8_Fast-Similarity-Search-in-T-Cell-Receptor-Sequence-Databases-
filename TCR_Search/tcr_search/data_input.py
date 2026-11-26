"""Read reference/query TSV records and their original JSONL annotations."""

import csv
import gzip
import json

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
STRATEGIES = ("spread", "prefix")


def require_integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def check_id(identifier):
    if not isinstance(identifier, str) or not identifier or identifier.strip() != identifier:
        raise ValueError("ID must be a nonempty string without outer whitespace")


def check_sequence(sequence, allow_empty=False):
    if not isinstance(sequence, str) or (not sequence and not allow_empty):
        raise ValueError("expected an amino acid sequence")
    for letter in sequence:
        if letter not in AMINO_ACIDS:
            raise ValueError(f"invalid amino acid: {letter!r}; use uppercase letters")


def read_records(path, id_column):
    """Read validated UTF-8 TSV rows as {ID: row}, retaining every column."""
    records = {}
    with open(path, encoding="utf-8-sig", newline="") as file:
        reader = csv.reader(file, delimiter="\t", strict=True)
        try:
            header = next(reader, [])
            if id_column not in header or "sequence" not in header:
                raise ValueError(f"requires {id_column} and sequence columns")
            if len(header) != len(set(header)):
                raise ValueError("duplicate column names")
            id_col, seq_col = header.index(id_column), header.index("sequence")
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError("malformed TSV row")
                rid, sequence = row[id_col], row[seq_col]
                check_id(rid)
                if rid in records:
                    raise ValueError(f"duplicate ID {rid}")
                check_sequence(sequence)
                records[rid] = dict(zip(header, row))
        except csv.Error as error:
            raise ValueError(f"{path}:{reader.line_num}: malformed TSV: {error}") from error
        except ValueError as error:
            raise ValueError(f"{path}:{reader.line_num}: {error}") from error
    return records


def read_sequences(path, id_column):
    """Read validated TSV records as {ID: sequence} for search algorithms."""
    return {rid: row["sequence"] for rid, row in read_records(path, id_column).items()}


def read_annotations(path, references):
    """Read all source records for the selected {reference_id: sequence} map.

    A full annotation file can serve any reference subset. Only selected IDs
    are retained and their records checked; the complete file is scanned to
    detect duplicate selected IDs. Gzip files are decompressed while reading.
    """
    annotations = {}
    opener = gzip.open if str(path).lower().endswith(".gz") else open
    line_num = 0
    try:
        with opener(path, "rt", encoding="utf-8-sig") as file:
            for line_num, line in enumerate(file, 1):
                if not line.strip():
                    continue
                item = json.loads(line)
                if not isinstance(item, dict):
                    raise ValueError("annotation row must be a JSON object")
                rid = item.get("reference_id")
                check_id(rid)
                if rid not in references:
                    continue
                if rid in annotations:
                    raise ValueError(f"duplicate annotations for reference ID {rid}")
                records = item.get("records")
                if not isinstance(records, list) or not records:
                    raise ValueError(f"reference {rid}: records must be a nonempty list")
                for record in records:
                    if not isinstance(record, dict):
                        raise ValueError(f"reference {rid}: each record must be a JSON object")
                    source_line = record.get("source_line")
                    if type(source_line) is not int or source_line < 1:
                        raise ValueError(f"reference {rid}: source_line must be a positive integer")
                    fields = record.get("fields")
                    if not isinstance(fields, dict) or not all(
                        isinstance(key, str) and isinstance(value, str)
                        for key, value in fields.items()
                    ):
                        raise ValueError(f"reference {rid}: fields must be a string-to-string object")
                    if fields.get("cdr3") != references[rid]:
                        raise ValueError(f"reference {rid}: annotation cdr3 does not match reference sequence")
                annotations[rid] = records
    except (OSError, EOFError, ValueError) as error:
        raise ValueError(f"{path}:{line_num}: {error}") from error
    missing = references.keys() - annotations.keys()
    if missing:
        examples = ", ".join(sorted(missing)[:3])
        raise ValueError(f"{path}:EOF: missing annotations for {len(missing)} reference ID(s): {examples}")
    return annotations
