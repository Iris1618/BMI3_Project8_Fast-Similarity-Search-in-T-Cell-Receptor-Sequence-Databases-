from inverted_index import InvertedIndex


def search_one(index: InvertedIndex, query_id: str, query: str, d: int,
               verifier, strategy: str = "spread", annotations=None) -> dict:
    """Produce STEP2_IO's final result using a teammate-provided verifier.

    verifier(query, reference_sequence, d) must return:
      * the exact nonnegative integer distance for every distance <= d;
      * None OR an integer > d for rejected pairs.
    A Boolean match indicator is not sufficient. The adapter cannot establish
    that a supplied verifier is correct; the group must validate final results.
    annotations maps reference_id to the original records list, if available.
    """
    if annotations is None:
        annotations = {}
    result = index.generate_candidates(query, d, strategy)
    matches = []
    for reference_id in sorted(result.candidate_ids):
        sequence = index.references[reference_id]
        distance = verifier(query, sequence, d)
        if distance is None:
            continue
        if type(distance) is not int or distance < 0:
            raise ValueError("verifier must return None or a nonnegative Python int")
        if distance > d:
            continue
        matches.append({
            "reference_id": reference_id,
            "sequence": sequence,
            "distance": distance,
            "annotations": annotations.get(reference_id, []),
        })
    return {"query_id": query_id, "matches": matches}
