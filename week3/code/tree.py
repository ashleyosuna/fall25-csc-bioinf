from typing import Optional, List
import numpy as np

class TreeNode:
    _is_root: bool
    _distance: float
    _parent: Optional[TreeNode]
    _index: int
    _children: Optional[List[TreeNode]]

    def __init__(self, children: Optional[List[TreeNode]]=None, distances: Optional[List[float]]=None, index=None):
        self._is_root = False
        self._distance = 0
        self._parent = None

        if index is None:
            if children is None or distances is None:
                raise TypeError("Either reference index (for terminal node) or "
                    "child nodes including the distance "
                    "(for intermediate node) must be set")
            if len(children) == 0:
                raise ValueError("Intermediate nodes must at least contain one child node")
            if len(children) != len(distances):
                raise ValueError(
                    "The number of children must equal the number of distances"
                )
            for i in range(len(children)):
                for j in range(len(children)):
                    if i != j and children[i] is children[j]:
                        raise ValueError(
                            "Two child nodes cannot be the same object"
                        )
            self._index = -1
            self._children = children
            for child, distance in zip(children, distances):
                child._set_parent(parent=self, distance=distance)
        elif index < 0:
            raise ValueError("Index cannot be negative")
        else:
            if children is not None or distances is not None:
                raise TypeError(
                    "Reference index and child nodes are mutually exclusive"
                )
            self._index = index
            self._children = None

    @property
    def index(self):
        return None if self._index == -1 else self._index
    
    @property
    def children(self):
        return self._children

    @property
    def parent(self):
        return self._parent
    
    @property
    def distance(self):
        return None if self._parent is None else self._distance

    def _set_parent(self, parent: TreeNode, distance: float):
        if self._parent is not None or self._is_root:
            raise ValueError("Node already has a parent")
        self._parent = parent
        self._distance = distance

    def __eq__(self, item: TreeNode):
        if self._distance != item._distance: return False
        elif self._index != -1 and self._index != item._index: return False
        elif frozenset(self._children) != frozenset(item._children): return False
        return True
    
    def is_leaf(self):
        return self._index != -1
    
    def is_root(self):
        return bool(self._is_root)
    
    def as_root(self):
        if self._parent is not None:
            raise ValueError("Node has parent, cannot be a root node")
        self._is_root = True
    
    def to_newick(self, labels: Optional[List[str]]=None, include_distance = True, round_distance: Optional[int] = None):
        if self.is_leaf():
            if labels is not None:
                for label in labels:
                    label = labels[self._index]
                    illegal_chars = [",",":",";","(",")"]
                    for char in illegal_chars:
                        if char in label:
                            raise ValueError(
                                f"Label '{label}' contains "
                                f"illegal character '{char}'"
                            )
            else:
                label = str(self._index)
            if include_distance:
                if round_distance is None:
                    return f"{label}:{self._distance:.1f}"
                else:
                    # TODO ?
                    rounded = round(self._distance, round_distance)
                    return f"{label}:{rounded}"
            else:
                return f"{label}"
        else:
            # Build string in a recursive way
            child_strings = [child.to_newick(
                labels, include_distance, round_distance
            ) for child in self._children]
            if include_distance:
                if round_distance is None:
                    return f"({','.join(child_strings)}):{self._distance:.1f}"
                else:
                    # TODO: fix not specifying number of decimal places in f string?
                    rounded = round(self._distance, round_distance)
                    return (
                        f"({','.join(child_strings)}):"
                        f"{rounded}"
                    )
            else:
                return f"({','.join(child_strings)})"
    
    def __str__(self):
        return self.to_newick(labels=None, include_distance=True, round_distance=None)
    
    def _get_leaves(node: TreeNode, leaf_list: List[TreeNode]):
        if node._index == -1:
            for child in node._children:
                TreeNode._get_leaves(child, leaf_list)
        else: leaf_list.append(node)

    def get_leaves(self):
        leaf_list: List[TreeNode] = []
        TreeNode._get_leaves(node=self, leaf_list=leaf_list)
        return leaf_list
    
    def _create_path_to_root(node: TreeNode):
        path: List[TreeNode] = []
        current_node: Optional[TreeNode] = node
        while current_node is not None:
            path.append(current_node)
            current_node = current_node._parent
        return path
    
    def lowest_common_ancestor(self, node: TreeNode):
        lca: Optional[TreeNode] = None
        self_path = TreeNode._create_path_to_root(self)
        other_path = TreeNode._create_path_to_root(node)
        for i in range(-1, -min(len(self_path), len(other_path)) - 1, -1):
            if self_path[i] is other_path[i]:
                lca = self_path[i]
            else:
                break
        return lca
    
    def distance_to(self, node: TreeNode, topological = False):
        distance = 0.0
        current_node: Optional[TreeNode] = None
        lca = self.lowest_common_ancestor(node)

        if lca is None:
            raise ValueError("The nodes do not have a common ancestor")
        current_node = self
        while current_node is not lca:
            if topological: distance += 1.0
            else: distance += current_node._distance

            current_node = current_node._parent
        
        current_node = node
        while current_node is not lca:
            if topological: distance += 1.0
            else: distance += current_node._distance

            current_node = current_node._parent
        
        return distance
    
class Tree:
    _root: TreeNode
    _leaves: List[Optional[TreeNode]]

    def __init__(self, root: TreeNode):
        root.as_root()
        self._root = root

        leaves_unsorted = self._root.get_leaves()
        leaf_count = len(leaves_unsorted)

        indices = np.array([leaf.index for leaf in leaves_unsorted])
        self._leaves = [None] * leaf_count
        for i in range(len(indices)):
            index = indices[i]
            if index >= leaf_count or index < 0:
                raise ValueError("The tree's indices are out of range")
            self._leaves[index] = leaves_unsorted[i]

    @property
    def root(self):
        return self._root
    
    @property
    def leaves(self):
        return self._leaves

    def to_newick(self, labels: Optional[List[str]] = None, include_distance = True, round_distance: Optional[int] = None):
        return self._root.to_newick(labels, include_distance, round_distance) + ";"
    
    def __str__(self):
        return self.to_newick()
    
    def __len__(self):
        return len(self._leaves)
    
    def __eq__(self, item: Tree):
        return self._root == item._root
    
    def get_distance(self, index1, index2, topological = False):
        return self._leaves[index1].distance_to(self._leaves[index2], topological)