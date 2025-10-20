import numpy as np

MIN = float('-inf')

def affine_alignment(u, v, match_score = 3, mismatch_score = -3, gap_init_score = -5, gap_extension_score = -3):
    n, m = len(u) + 1, len(v) + 1

    lower = np.ndarray((n, m))
    middle = np.ndarray((n, m))
    upper = np.ndarray((n, m))

    lower[0][0] = MIN
    middle[0][0] = 0
    upper[0][0] = MIN
    # initialize first column
    for i in range(1, n): lower[i][0] = (gap_init_score if i == 1 else lower[i-1][0] + gap_extension_score)

    # initialize first row
    for j in range(1, m): upper[0][j] = (gap_init_score if j == 1 else upper[0][j-1] + gap_extension_score)

    # construct matrices
    for i in range(1, n):
        for j in range(1, m):
            # lower
            lower[i][j] = max(
                lower[i-1][j] + gap_extension_score,
                middle[i-1][j] + gap_init_score
            )

            upper[i][j] = max(
                upper[i][j-1] + gap_extension_score,
                middle[i][j-1] + gap_init_score
            )

            middle[i][j] = max(
                lower[i][j],
                middle[i-1][j-1] + (match_score if u[i-1] == v[j-1] else mismatch_score),
                upper[i][j]
            )

    # backtrace
    i, j = n-1, m-1
    alignment = [[], []]

    while i > 0 or j > 0:
        if middle[i][j] == middle[i-1][j-1] + mismatch_score \
            or u[i-1] == v[j-1]:
            i -= 1
            j -= 1
            alignment[0].append(u[i])
            alignment[1].append(v[j])
        # vertical gap
        elif middle[i][j] == lower[i][j]:
            i -= 1
            alignment[0].append(u[i])
            alignment[1].append("-")
        else:
            j -= 1
            alignment[0].append("-")
            alignment[1].append(v[j])
    
    alignment[0].reverse()
    alignment[1].reverse()
    return alignment