"""Return the edit distance if it is <= max_distance.
   Return None if the distance is larger.
"""
def verify(query: str, reference: str, max_distance: int) -> int | None:
    
    # Return None if the lengths differ by more than max_distance.
    if abs(len(query) - len(reference)) > max_distance:
        return None

    # Keep reference as the shorter sequence, which reduces the memory needed.
    if len(reference) > len(query):
        query, reference = reference, query

    reference_length = len(reference)
    outside_band = max_distance + 1

    # Rolling-row DP: store only two rows
    previous_row = [outside_band] * (reference_length + 1)
    current_row = [outside_band] * (reference_length + 1)

    for j in range(min(reference_length, max_distance) + 1):
        previous_row[j] = j

    for i, query_char in enumerate(query, start=1):
        # Banded DP: only calculate cells near the main diagonal
        first_column = max(0, i - max_distance)
        last_column = min(reference_length, i + max_distance)

        if first_column == 0:
            current_row[0] = i
        else:
            current_row[first_column - 1] = outside_band
        if last_column < reference_length:
            current_row[last_column + 1] = outside_band

        # Dynamic programming algorithm
        for j in range(max(1, first_column), last_column + 1):
            cost = 0 if query_char == reference[j - 1] else 1
            deletion = previous_row[j] + 1
            insertion = current_row[j - 1] + 1
            match_or_substitution = previous_row[j - 1] + cost
            current_row[j] = min(deletion, insertion, match_or_substitution)

        previous_row, current_row = current_row, previous_row

    distance = previous_row[-1]
    return distance if distance <= max_distance else None
