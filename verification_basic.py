def verification_basic(query: str, reference: str) -> int:
    m = len(query)
    n = len(reference)
    distance = [[0] * (n + 1) for _ in range(m + 1)]

    for i in range(m + 1):
        distance[i][0] = i

    for j in range(n + 1):
        distance[0][j] = j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            c = 0 if query[i - 1] == reference[j - 1] else 1
            deletion = distance[i - 1][j] + 1
            insertion = distance[i][j - 1] + 1
            match_or_substitution = distance[i - 1][j - 1] + c
            distance[i][j] = min(deletion, insertion, match_or_substitution)

    return distance[m][n]
