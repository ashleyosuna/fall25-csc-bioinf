def reverse_complement(key):
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    key = list(key[::-1])
    for i in range(len(key)):
        key[i] = complements[key[i]]
    return "".join(key)

class Node:
    kmer: str
    _children: set[int]
    _count: int
    visited: bool
    depth: int
    # max_depth_child:
    def __init__(self, kmer: str):
        self.kmer = kmer
        self._children = set()
        self._count = 0
        self.visited = False
        self.depth = 0
        # self.max_depth_child = None
    
    def add_child(self, kmer):
        self._children.add(kmer)
    
    def increase(self):
        self._count += 1

    def reset(self):
        self.visited = False
        self.depth = 0
        # self.max_depth_child = None

class DBG:
    k: int
    kmer2idx: dict[str, int]
    nodes: dict[int, Node]
    kmer_count: int

    def __init__(self, k: int, data_list: List[List[str]]):
        self.k = k
        self.nodes = {}
        self.kmer2idx = {}
        self.kmer_count = 0

        self._check(data_list)
        self._build(data_list)

    def _check(self, data_list):
        assert len(data_list) > 0
        assert self.k <= len(data_list[0][0])

    def _build(self, data_list: List[List[str]]):
        for data in data_list:
            # reads
            for og in data:
                rc = reverse_complement(og)
                # adding edges between consecutive kmers and reverse complememts
                for i in range(len(og) - self.k  - 1):
                    a, b = og[i : i + self.k], og[i + 1 : i + 1 + self.k]
                    self._add_arc(a, b)
                    c, d = rc[i : i + self.k], rc[i + 1 : i + 1 + self.k]
                    self._add_arc(c, d)
                    return
                    
    def _add_node(self, kmer: str):
        if kmer not in self.kmer2idx:
            self.kmer2idx[kmer] = self.kmer_count
            self.nodes[self.kmer_count] = Node(kmer)
            self.kmer_count += 1
        idx = self.kmer2idx[kmer]
        self.nodes[idx].increase()
        return idx
    
    def _add_arc(self, kmer1: str, kmer2: str):
        idx1 = self._add_node(kmer1)
        idx2 = self._add_node(kmer2)
        self.nodes[idx1].add_child(idx2)
    
    def _reset(self):
        for idx in self.nodes.keys():
            self.nodes[idx].reset()
    
    def get_longest_contig(self):
        self._reset()
        return
