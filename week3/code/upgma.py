from tree import Tree, TreeNode
import numpy as np

MAX_FLOAT = np.finfo(float).max

def upgma(distances: np.ndarray):
    if distances.shape[0] != distances.shape[1] \
    or not np.allclose(distances.T, distances):
        raise ValueError("Distance matrix must be symmetric")
    if np.isnan(distances).any():
        raise ValueError("Distance matrix contains NaN values")
    if (distances >= MAX_FLOAT).any():
        raise ValueError("Distance matrix contains infinity")
    if (distances < 0).any():
        raise ValueError("Distances must be positive")
    
    nodes = np.array([TreeNode(index=i) for i in range(distances.shape[0])])
    is_clustered_v = np.full(distances.shape[0], False, dtype=np.uint8)
    cluster_size_v = np.ones(distances.shape[0], dtype=float)
    node_heights = np.zeros(distances.shape[0], dtype=float)
    # distances_v = distances.astype(np.float32, copy=True)
    distances_v = distances.astype(float, copy=True)

    while True:
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        for i in range(distances_v.shape[0]):
            if is_clustered_v[i]:
                continue
            for j in range(i):
                if is_clustered_v[j]:
                    continue
                dist = distances_v[i,j]
                if dist < dist_min:
                    dist_min = dist
                    i_min = i
                    j_min = j
        
        if i_min == -1 or j_min == -1:
            # No distance found -> all leaf nodes are clustered
            # -> exit loop
            break

        height = dist_min/2
        nodes[i_min] = TreeNode(
            [nodes[i_min], nodes[j_min]],
            [height-node_heights[i_min], height-node_heights[j_min]]
        )
        node_heights[i_min] = height

        nodes[j_min] = None
        is_clustered_v[j_min] = True
        # Calculate arithmetic mean distances of child nodes
        # as distances for new node and update matrix
        for k in range(distances_v.shape[0]):
            if not is_clustered_v[k] and k != i_min:
                mean = (
                    (
                          distances_v[i_min,k] * cluster_size_v[i_min]
                        + distances_v[j_min,k] * cluster_size_v[j_min]
                    ) / (cluster_size_v[i_min] + cluster_size_v[j_min])
                )
                distances_v[i_min,k] = mean
                distances_v[k,i_min] = mean
        # Updating cluster size of new node
        cluster_size_v[i_min] = cluster_size_v[i_min] + cluster_size_v[j_min]
    

    # As each higher level node is always created on position i_min
    # and i is always higher than j in minimum distance calculation,
    # the root node must be at the last index
    return Tree(nodes[len(nodes)-1])