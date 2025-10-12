import numpy as np
import unittest
import biotite.sequence.phylo as phylo
import os

distances = np.loadtxt("week3/data/distances.txt", dtype=int)
upgma_newick = None
with open("week3/data/newick_upgma.txt") as f:
    upgma_newick = f.read().strip()

tree = phylo.upgma(distances)

class Test(unittest.TestCase):
    def test_distances(self):
    # Tree is created via UPGMA
    # -> The distances to root should be equal for all leaf nodes
        dist = tree.root.distance_to(tree.leaves[0])
        for leaf in tree.leaves:
            self.assertEqual(leaf.distance_to(tree.root), dist)
        # Example topological distances
        self.assertEqual(tree.get_distance(0, 19, True), 9)
        self.assertEqual(tree.get_distance(4, 2, True), 10)

    def test_neighbor_joining(self):
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

        ref_tree = phylo.Tree(
            phylo.TreeNode(
                [
                    phylo.TreeNode(
                        [
                            phylo.TreeNode(
                                [
                                    phylo.TreeNode(index=0),
                                    phylo.TreeNode(index=1),
                                ],
                                [1.0, 4.0],
                            ),
                            phylo.TreeNode(index=2),
                        ],
                        [1.0, 2.0],
                    ),
                    phylo.TreeNode(
                        [
                            phylo.TreeNode(index=3),
                            phylo.TreeNode(index=4),
                        ],
                        [3.0, 2.0],
                    ),
                    phylo.TreeNode(index=5),
                ],
                [1.0, 1.0, 5.0],
            )
        )

        test_tree = phylo.neighbor_joining(dist)
        self.assertEqual(test_tree, ref_tree)

    def test_upgma(self):
        """
        Compare the results of `upgma()` with DendroUPGMA.
        """
        print('in test upgma')
        ref_tree = phylo.Tree.from_newick(upgma_newick)
        # Cannot apply direct tree equality assertion because the distance
        # might not be exactly equal due to floating point rounding errors
        for i in range(len(tree)):
            for j in range(len(tree)):
                # Check for equal distances and equal topologies
                self.assertAlmostEqual(tree.get_distance(i, j), ref_tree.get_distance(i, j), delta=1e-3)
                self.assertEqual(tree.get_distance(i, j, topological=True), ref_tree.get_distance(i, j, topological=True))

tests = Test()
# runner = unittest.TextTestRunner(verbosity=0)
# unittest.main(testRunner=runner)
tests.test_distances()
tests.test_neighbor_joining()
tests.test_upgma()