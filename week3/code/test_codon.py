from typing import Tuple, NoneType
from tree import Tree, TreeNode
import numpy as np
from upgma import upgma
from nj import neighbor_joining
import time
import unittest

# TEST UTILITIES

tests: list[tuple[str, function[Tuple, NoneType]]] = []

def test(name: str):
    def decorator(fn: function[Tuple, NoneType]) -> function[Tuple, NoneType]:
        tests.append((name, fn))
        return fn
    return decorator

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
TEST = unittest.TestCase()

# END OF TEST UTILITIES

# TEST CASES
@test("test_distances")
def test_distances():
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
    dist = tree.root.distance_to(tree.leaves[0])
    for leaf in tree.leaves:
        TEST.assertEqual(leaf.distance_to(tree.root), dist)
    # Example topological distances
    TEST.assertEqual(tree.get_distance(0, 19, True), 9)
    TEST.assertEqual(tree.get_distance(4, 2, True), 10)

@test("test_neighbor_joining")
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
                            [1, 4],
                        ),
                        TreeNode(index=2),
                    ],
                    [1, 2],
                ),
                TreeNode(
                    [
                        TreeNode(index=3),
                        TreeNode(index=4),
                    ],
                    [3, 2],
                ),
                TreeNode(index=5),
            ],
            [1, 1, 5],
        )
    )

    test_tree = neighbor_joining(dist)
    TEST.assertEqual(test_tree, ref_tree)

@test("test_upgma")
def test_upgma():
    """
    Compare the results of `upgma()` with DendroUPGMA.
    """
    ref_tree = Tree.from_newick(upgma_newick)
    # Cannot apply direct tree equality assertion because the distance
    # might not be exactly equal due to floating point rounding errors
    for i in range(len(tree)):
        for j in range(len(tree)):
            # Check for equal distances and equal topologies
            TEST.assertAlmostEqual(tree.get_distance(i, j), ref_tree.get_distance(i, j), delta=1e-3)
            TEST.assertEqual(tree.get_distance(i, j, topological=True), ref_tree.get_distance(i, j, topological=True))
# END TEST CASES

def testRunner():
    passed = total = 0
    start = time.perf_counter()
    for name, test in tests:
        try:
            test()
            # print(f"{name} ... ok")
            passed += 1
        except AssertionError:
            # print(f"{name} ... not ok")
            pass
        except:
            # print(f"{name} ... ERROR")
            pass
        total += 1
    end = time.perf_counter()
    elapsed_time = (end - start) * 1000
    print(f"{elapsed_time} ms")

testRunner()