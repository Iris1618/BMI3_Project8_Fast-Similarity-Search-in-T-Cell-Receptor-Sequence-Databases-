"""Read the group's reference/query TSV format.

The reader is moved from Inverted_Index/inverted_index.py. Input preparation
and downloading real datasets are deliberately outside this example project.
"""

import csv

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


def read_sequences(path, id_column):
    """Read UTF-8 TSV records as {ID: sequence}, validating every input row."""
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
                records[rid] = sequence
        except csv.Error as error:
            raise ValueError(f"{path}:{reader.line_num}: malformed TSV: {error}") from error
        except ValueError as error:
            raise ValueError(f"{path}:{reader.line_num}: {error}") from error
    return records
