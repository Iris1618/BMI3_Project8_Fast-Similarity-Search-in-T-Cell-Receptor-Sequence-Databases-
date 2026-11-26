"""Length-filtered exhaustive search using the group's basic verifier."""

from verification_basic import verification_basic


def exhaustive_matches(query, references, d):
    """Return {reference_id: exact_distance} for all matches within d.

    references maps the group's reference IDs to their sequences. The length
    filter is safe because insertions/deletions are required to change length.
    """
    if type(d) is not int or d < 0:
        raise ValueError("d must be a non-negative integer")
    matches = {}
    for reference_id, sequence in references.items():
        if abs(len(query) - len(sequence)) > d:
            continue
        distance = verification_basic(query, sequence)
        if distance <= d:
            matches[reference_id] = distance
    return matches


def search_one(query_id, query, references, d, annotations=None):
    """Return the group's STEP2_IO query/matches output format."""
    annotations = {} if annotations is None else annotations
    matches = exhaustive_matches(query, references, d)
    return {
        "query_id": query_id,
        "matches": [
            {"reference_id": rid, "sequence": references[rid],
             "distance": matches[rid], "annotations": annotations.get(rid, [])}
            for rid in sorted(matches)
        ],
    }
