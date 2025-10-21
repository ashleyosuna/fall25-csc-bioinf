import numpy as np

def fitting(u, v, match_score = 3, mismatch_score = -3, gap_score = -2):
    u = u.upper()
    v = v.upper()
    # make it so |v| < |u|
    if len(u) < len(v): u, v = v, u
    n, m = len(u) + 1, len(v) + 1

    network = np.full((n, m), 0.0, dtype=float)
    # ignore gaps at the left of v -> first column stays as 0s
    # initialize first row as usual
    for i in range(m): network[0][i] = gap_score * i

    # ignore gaps at the right of v -> keep track of highest score
    max_score = float('-inf')
    max_score_position = (0, 0)

    # fill matrix row by row
    for i in range(1, n):
        for j in range(1, m):
            network[i][j] = max(
                network[i-1][j] + gap_score,
                network[i][j-1] + gap_score,
                network[i-1][j-1] + (match_score if u[i-1] == v[j-1] else mismatch_score),
            )

            if j == m - 1 and network[i][j] > max_score:
                max_score = network[i][j]
                max_score_position = (i, j)

    # backtrack from highest score
    alignment = [[char for char in u], ['-' for _ in u]]
    i, j = max_score_position
    k = max(i - 1, j- 1)

    while j > 0:
        if (j > 0 and i <= 0) and network[i][j] == network[i][j-1] + gap_score:
            j -= 1
            alignment[1][k] = v[j]
        elif (i > 0 and j <= 0) or network[i][j] == network[i-1][j] + gap_score:
            i -= 1
        else:
            i -= 1
            j -= 1
            alignment[1][k] = v[j]
        k -= 1

    return alignment