import sys
from utils import read_data
from dbg import DBG
# from python import os
# from python import sys as sysp

# sysp.setrecursionlimit(1000000)

argv = sys.argv
short1, short2, long = read_data(argv[1])
k = 25
dbg = DBG(k=k, data_list=[short1, short2, long])

with open(argv[1] + '/contig_codon.fasta', 'w') as f:
    for i in range(20):
        c = dbg.get_longest_contig()
        if c is None:
            break
        print(i, len(c))
        f.write('>contig_' + str(i) + '\n')
        f.write(c + '\n')