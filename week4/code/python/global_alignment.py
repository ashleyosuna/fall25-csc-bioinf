import numpy as np

def global_alignment(u, v):
    # CONSTRUCT MATRIX
    n, m = len(u) + 1, len(v) + 1
    network = np.ndarray((n, m))

    # fill first row
    for i in range(n): network[i][0] = 0

    # fill first column
    for j in range(m): network[0][j] = 0

    # fill matrix row by row
    for i in range(1, n):
        for j in range(1, m):
            network[i][j] = max(
                network[i-1][j],
                network[i][j-1],
                network[i-1][j-1] + (1 if u[i-1] == v[j-1] else 0)
            )

    # BACKTRACK TO CONSTRUCT ALIGNMENT
    i, j = n - 1, m - 1
    k = max(n, m)
    alignment = np.array([[""] * k, [""] * k])
    k -= 1

    while i > 0 or j > 0:
        if network[i][j] == network[i-1][j-1] or \
            (u[i-1] == v[j-1] and network[i][j] == network[i-1][j-1] + 1):
            i -= 1
            j -= 1
            alignment[0][k] = u[i]
            alignment[1][k] = v[j]
        elif network[i][j] == network[i-1][j]: 
            i -= 1
            alignment[0][k] = u[i]
            alignment[1][k] = "-"
        else: 
            j -= 1
            alignment[0][k] = "-"
            alignment[1][k] = v[i]
        k -= 1

    # construct alignment
    return alignment

al = global_alignment("ATGTTATA", "ATCGTCC")
global_alignment("ATGTTA", "ATCGT")
    