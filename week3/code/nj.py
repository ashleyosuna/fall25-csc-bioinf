from tree import Tree, TreeNode
import numpy as np

MAX_FLOAT = np.finfo(float).max

def neighbor_joining(distances: np.ndarray):
    i = j = k = u = 0
    i_min = j_min = 0
    dist: float = 0
    dist_sum: float = 0
    dist_min = MAX_FLOAT
    node_dist_i = node_dist_j = node_dist_k = 0

    if distances.shape[0] != distances.shape[1] \
        or not np.allclose(distances.T, distances):
            raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if distances.shape[0] < 4:
        raise ValueError("At least 4 nodes are required")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")
    
    # Keep track on clustered indices
    nodes = np.array(
        [TreeNode(index=i) for i in range(distances.shape[0])]
    )

    is_clustered_v = np.full(
        distances.shape[0], False, dtype=np.uint8
    )

    n_rem_nodes = \
        len(distances) - np.count_nonzero(np.asarray(is_clustered_v))
    
    divergence_v = np.zeros(
        distances.shape[0], dtype=float
    )

    corr_distances_v = np.zeros(
        (distances.shape[0],) * 2, dtype=float
    )

    distances_v = distances.astype(float, copy=True)

    while True:

        # Calculate divergence
        for i in range(distances_v.shape[0]):
            if is_clustered_v[i]:
                continue
            dist_sum = 0.0
            for k in range(distances_v.shape[0]):
                if is_clustered_v[k]:
                    continue
                dist_sum += distances_v[i,k]
            divergence_v[i] = dist_sum
        
        # Calculate corrected distance matrix
        for i in range(distances_v.shape[0]):
            if is_clustered_v[i]:
                    continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                corr_distances_v[i,j] = \
                    (n_rem_nodes - 2) * distances_v[i,j] \
                    - divergence_v[i] - divergence_v[j]

        # Find minimum corrected distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(corr_distances_v.shape[0]):
            if is_clustered_v[i]:
                    continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                dist = corr_distances_v[i,j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        # Check if all nodes have been clustered
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            # -> exit loop
            break
        
        # Cluster the nodes with minimum distance
        # replacing the node at position i_min
        # leaving the node at position j_min empty
        # (is_clustered_v -> True)
        node_dist_i = 0.5 * (
            distances_v[i_min,j_min]
            + 1/(n_rem_nodes-2) * (divergence_v[i_min] - divergence_v[j_min])
        )
        node_dist_j = 0.5 * (
            distances_v[i_min,j_min]
            + 1/(n_rem_nodes-2) * (divergence_v[j_min] - divergence_v[i_min])
        )
        if n_rem_nodes > 3:
            # Clustering is not finished
            # -> Create a node with two children
            nodes[i_min] = TreeNode(
                [nodes[i_min], nodes[j_min]],
                [node_dist_i, node_dist_j]
            )
            # Mark position j_min as clustered
            nodes[j_min] = None
            is_clustered_v[j_min] = True
        else:
            # Clustering is finished
            # Combine ast three nodes into root node
            # Find the index of the remaining one of the three nodes
            # (other than i_min and j_min)
            is_clustered_v[i_min] = True
            is_clustered_v[j_min] = True
            # The index of the remaining one
            k = np.where(~np.asarray(is_clustered_v, dtype=bool))[0][0]
            node_dist_k = 0.5 * (
                distances_v[i_min,k] + distances_v[j_min,k]
                - distances_v[i_min,j_min]
            )
            root = TreeNode(
                [nodes[i_min], nodes[j_min], nodes[k]],
                [node_dist_i, node_dist_j, node_dist_k]
            )
            # Clustering is finished -> put into tree and return
            return Tree(root)
        
        # Update distance matrix
        # Calculate distances of new node to all other nodes
        for k in range(distances_v.shape[0]):
            if not is_clustered_v[k] and k != i_min:
                dist = 0.5 * (
                    distances_v[i_min,k] + distances_v[j_min,k]
                    - distances_v[i_min,j_min]
                )
                distances_v[i_min,k] = dist
                distances_v[k,i_min] = dist

        # Update the amount of remaining nodes
        n_rem_nodes = \
            len(distances) - np.count_nonzero(np.asarray(is_clustered_v))