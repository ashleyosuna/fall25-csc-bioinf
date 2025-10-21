import numpy as np

MIN = float('-inf')

def affine_alignment(u, v, match_score = 3, mismatch_score = -3, gap_init_score = -5, gap_extension_score = -1):
    u = u.upper()
    v = v.upper()
    n, m = len(u) + 1, len(v) + 1

    lower = np.full((n, m), MIN, dtype=float)
    middle = np.full((n, m), 0.0, dtype=float)
    upper = np.full((n, m), MIN, dtype=float)

    lower[0][0] = 0.0
    middle[0][0] = 0.0
    upper[0][0] = 0.0
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

    # print(lower, middle, upper, sep="\n")

    while i > 0 or j > 0:
        # vertical gap
        if (i > 0 and j <= 0) or middle[i][j] == lower[i][j]:
            i -= 1
            alignment[0].append(u[i])
            alignment[1].append("-")
        # horizontal gap
        elif (j > 0 and i <= 0) or middle[i][j] == upper[i][j]:
            j -= 1
            alignment[0].append("-")
            alignment[1].append(v[j])
        # match/mismatch
        else:
            i -= 1
            j -= 1
            alignment[0].append(u[i])
            alignment[1].append(v[j])
    
    alignment[0].reverse()
    alignment[1].reverse()
    return alignment

# al = affine_alignment("ACGTGC", "ACGT")
# print(al)