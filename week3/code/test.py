from tree import Tree
import numpy as np
from upgma import upgma

def test_distances(tree):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        assert leaf.distance_to(tree.root) == dist
    # Example topological distances
    assert tree.get_distance(0, 19, True) == 9
    assert tree.get_distance(4, 2, True) == 10

def test_neighbor_joining():
    return

def test_upgma():
    return

distances = []

with open("week3/data/distances.txt") as f:
    for row in f:
        row = row.strip()
        _row = []
        cols = row.split()
        for col in cols:
            _row.append(float(col))
        distances.append(_row)


tree = upgma(np.array(distances))
test_distances(tree)