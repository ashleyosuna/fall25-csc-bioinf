from tree import Tree, TreeNode
import numpy as np
from upgma import upgma
from nj import neighbor_joining

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
    """
    Compare the results of `neighbor_join()` with a known tree.
    """
    dist = np.array([
        [ 0,  5,  4,  7,  6,  8],
        [ 5,  0,  7, 10,  9, 11],
        [ 4,  7,  0,  7,  6,  8],
        [ 7, 10,  7,  0,  5,  9],
        [ 6,  9,  6,  5,  0,  8],
        [ 8, 11,  8,  9,  8,  0],
    ])  # fmt: skip

    ref_tree = Tree(
        TreeNode(
            [
                TreeNode(
                    [
                        TreeNode(
                            [
                                TreeNode(index=0),
                                TreeNode(index=1),
                            ],
                            [1.0, 4.0],
                        ),
                        TreeNode(index=2),
                    ],
                    [1.0, 2.0],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3.0, 2.0],
                ),
                TreeNode(index=5),
            ],
            [1.0, 1.0, 5.0],
        )
    )

    test_tree = neighbor_joining(dist)
    assert test_tree == ref_tree

def test_upgma(tree, upgma_newick):
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            # assert tree.get_distance(i, j) == pytest.approx(
            #     ref_tree.get_distance(i, j), abs=1e-3
            # )
            assert np.isclose(tree.get_distance(i, j), ref_tree.get_distance(i, j), atol=1e-3)
            assert tree.get_distance(i, j, topological=True) == ref_tree.get_distance(
                i, j, topological=True
            )

distances = []

with open("week3/data/distances.txt") as f:
    for row in f:
        row = row.strip()
        _row = []
        cols = row.split()
        for col in cols:
            _row.append(float(col))
        distances.append(_row)

upgma_newick = None

with open("week3/data/newick_upgma.txt") as f:
    upgma_newick = f.read().strip()

tree = upgma(np.array(distances))
test_distances(tree)
test_upgma(tree, upgma_newick)
test_neighbor_joining()