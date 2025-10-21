import numpy as np

def global_alignment(u, v, match_score = 3, mismatch_score = -3, gap_score = -2):
    u = u.upper()
    v = v.upper()
    # CONSTRUCT MATRIX
    n, m = len(u) + 1, len(v) + 1
    network = np.full((n, m), 0, dtype=float)

    # fill first row
    for i in range(n): network[i][0] = gap_score * i

    # fill first column
    for j in range(m): network[0][j] = gap_score * j

    # fill matrix row by row
    for i in range(1, n):
        for j in range(1, m):
            network[i][j] = max(
                network[i-1][j] + gap_score,
                network[i][j-1] + gap_score,
                network[i-1][j-1] + (match_score if u[i-1] == v[j-1] else mismatch_score)
            )

    # BACKTRACK TO CONSTRUCT ALIGNMENT
    alignment = [[], []]
    i, j = n - 1, m - 1

    while i > 0 or j > 0:
        # vertical gap, 'move' along u
        if network[i][j] == network[i-1][j] + gap_score:
            i -= 1
            alignment[0].append(u[i])
            alignment[1].append('-')
        # horizontal gap, 'move' along v
        elif network[i][j] == network[i][j-1] + gap_score:
            j -= 1
            alignment[0].append("-")
            alignment[1].append(v[j])
        # if match or mismatch, 'move' along both sequences
        else:
            i -= 1
            j -= 1
            alignment[0].append(u[i])
            alignment[1].append(v[j])
        
    alignment[0].reverse()
    alignment[1].reverse()
    return alignment

# al = global_alignment("ACGTGC", "ACGT")
# print(al)