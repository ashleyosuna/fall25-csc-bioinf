import numpy as np

def fitting(u, v, match_score = 3, mismatch_score = -3, gap_score = -2):
    # assuming |v| < |u|
    n, m = len(u) + 1, len(v) + 1

    network = np.array([[0] * m] * n)
    # ignore gaps at the left of v -> first column stays as 0s
    # initialize first row as usual
    for i in range(m): network[0][i] = gap_score * i

    # ignore gaps at the right of v -> keep track of highest score
    max_score = 0
    max_score_position = (0, 0)

    # fill matrix row by row
    for i in range(1, n):
        for j in range(1, m):
            network[i][j] = max(
                network[i-1][j] + gap_score,
                network[i][j-1] + gap_score,
                network[i-1][j-1] + (match_score if u[i-1] == v[j-1] else mismatch_score),
                0
            )

            if network[i][j] > max_score:
                max_score = network[i][j]
                max_score_position = (i, j)
    
    print(network)
    
    # backtrack from highest score
    alignment = [[_ for _ in u], ["-"] * (n - 1)]
    i, j = max_score_position
    k = j + 1

    while i > 0 or j > 0:
        if network[i][j] == network[i-1][j-1] + mismatch_score or \
            u[i-1] == v[j-1]:
            i -= 1
            j -= 1
            alignment[1][k] = v[j]
        elif network[i][j] == network[i-1][j] + gap_score:
            i -= 1
        elif network[i][j] == network[i][j-1] + gap_score:
            j -= 1
            alignment[1][k] = v[j]
        # current position is where local alignment starts
        else:
            break
        k -= 1

    # if v does not fit entirely into u (i.e., alignment start pos + |v| > len|u|)
    # then add gaps at the end of u
    if max_score_position[1] < m - 1:
        alignment[0] += ["-"] * (m - 1 - max_score_position[1])
        alignment[1] += [_ for _ in v[max_score_position[1]:]]

    return alignment