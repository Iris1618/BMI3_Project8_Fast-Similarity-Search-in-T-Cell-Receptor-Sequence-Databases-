"""Length-filtered exhaustive search using the final verification function."""

from .verification import verify


def exhaustive_matches(query, references, d):
    """Check every length-eligible reference and return {ID: exact_distance}."""
    if type(d) is not int or d < 0:
        raise ValueError("d must be a non-negative integer")
    matches = {}
    for reference_id, sequence in references.items():
        if abs(len(query) - len(sequence)) > d:
            continue
        distance = verify(query, sequence, d)
        if distance is not None:
            matches[reference_id] = distance
    return matches
