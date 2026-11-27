"""Read prepared reference/query TSV files and check their format."""

import csv

# The 20 allowed uppercase amino-acid letters.
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
STRATEGIES = ("spread", "prefix")


def require_integer(value, name, minimum):
    # Use type() so that True and False are not accepted as integers.
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")


def check_id(identifier):
    if not isinstance(identifier, str):
        raise ValueError("ID must be a nonempty string without outer whitespace")
    if identifier == "" or identifier.strip() != identifier:
        raise ValueError("ID must be a nonempty string without outer whitespace")


def check_sequence(sequence, allow_empty=False):
    if not isinstance(sequence, str):
        raise ValueError("expected an amino acid sequence")
    # The index may allow empty sequences; file input does not allow them.
    if sequence == "" and not allow_empty:
        raise ValueError("expected an amino acid sequence")
    for letter in sequence:
        if letter not in AMINO_ACIDS:
            raise ValueError(f"invalid amino acid: {letter!r}; use uppercase letters")


def read_sequences(path, id_column):
    """Return {ID: sequence}; id_column is reference_id or query_id."""
    records = {}
    # utf-8-sig also handles an optional UTF-8 marker at the start of the file.
    with open(path, encoding="utf-8-sig", newline="") as file:
        # TSV columns are separated by tabs.
        reader = csv.reader(file, delimiter="\t", strict=True)
        try:
            # 1. Read the header and locate the two required columns.
            header = next(reader, [])
            if id_column not in header or "sequence" not in header:
                raise ValueError(f"requires {id_column} and sequence columns")
            if len(header) != len(set(header)):
                raise ValueError("duplicate column names")
            id_col = header.index(id_column)
            seq_col = header.index("sequence")

            # 2. Skip blank lines, then check and store each data row.
            for row in reader:
                if not row:
                    continue
                if len(row) != len(header):
                    raise ValueError("malformed TSV row")
                identifier = row[id_col]
                sequence = row[seq_col]
                check_id(identifier)
                # IDs must be unique; different IDs may have the same sequence.
                if identifier in records:
                    raise ValueError(f"duplicate ID {identifier}")
                check_sequence(sequence)
                records[identifier] = sequence
        except csv.Error as error:
            # Add the file path and line number, keeping the original cause.
            raise ValueError(f"{path}:{reader.line_num}: malformed TSV: {error}") from error
        except ValueError as error:
            raise ValueError(f"{path}:{reader.line_num}: {error}") from error
    return records
