import numpy as np

def local_alignment(u, v, match_score = 3, mismatch_score = -3, gap_score = -2):
    # CONSTRUCT MATRIX
    n, m = len(u) + 1, len(v) + 1
    network = np.ndarray((n, m))

    # keep track of where best local alignment ends so we can backtrack later
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

    # BACKTRACK TO CONSTRUCT ALIGNMENT
    alignment = [[], []]
    i, j = max_score_position

    while i > 0 or j > 0:
        if network[i][j] == network[i-1][j-1] + mismatch_score or \
            u[i-1] == v[j-1]:
            i -= 1
            j -= 1
            alignment[0].append(u[i])
            alignment[1].append(v[j])
        elif network[i][j] == network[i-1][j] + gap_score:
            i -= 1
            alignment[0].append(u[i])
            alignment[1].append('-')
        elif network[i][j] == network[i][j-1] + gap_score:
            j -= 1
            alignment[0].append("-")
            alignment[1].append(v[j])
        # current position is where local alignment starts
        else:
            break
        
    alignment[0].reverse()
    alignment[1].reverse()
    return alignment

al = local_alignment("AAGTG", "GGCTG")
print(al[0], al[1], sep='\n')