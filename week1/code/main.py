import sys
from utils import read_data
from dbg import DBG

argv = sys.argv
# TODO: use os to joing path
short1, short2, long = read_data(argv[1])
k = 25
dbg = DBG(k=k, data_list=[short1, short2, long])
# TODO: use os to join path
with open(argv[1] + '/' + 'contig.fasta', 'w') as f:
    for i in range(20):
        c = dbg.get_longest_contig()
        if c is None:
            break
        print(i, len(c))
        # f.write('>contig_%d\n' % i)
        # f.write(c + '\n')