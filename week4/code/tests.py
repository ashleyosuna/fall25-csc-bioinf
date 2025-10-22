from global_alignment import global_alignment
from local_alignment import local_alignment
from affine_alignment import affine_alignment
from fitting_alignment import fitting
import time
import sys

arg = sys.argv

lang = arg[1]

def read_fasta(file):
    sequences = []
    curr_seq_name = ''
    curr_seq = ''
    with open(file) as f:
        for line in f:
            if line[0] == '>':
                if len(curr_seq): sequences.append((curr_seq_name, curr_seq))
                curr_seq_name = line[1:].strip()
                curr_seq = ''
            else: 
                curr_seq += line.strip()
    if len(curr_seq): sequences.append((curr_seq_name, curr_seq))
    return sequences

# test q1 vs t1, ...., q5 vs q5
sequencesA = read_fasta('week4/code/data/q1.fa') + read_fasta('week4/code/data/MT-human.fa')
sequencesB = read_fasta('week4/code/data/t1.fa') + read_fasta('week4/code/data/MT-orang.fa')

for i in range(len(sequencesA)):
    (nameA, seqA), (nameB, seqB) = sequencesA[i], sequencesB[i]

    start = time.time()
    al = global_alignment(seqA, seqB)

    # print(''.join(al[0]), ''.join(al[1]), sep='\n')
    end = time.time()

    elapsed = (end - start) * 1000
    print(f"{'global-' + nameA:<20} {lang:<15} {elapsed:.4f} ms")

    start = time.time()
    al = local_alignment(seqA, seqB)
    # print(''.join(al[0]), ''.join(al[1]), sep='\n')
    end = time.time()
    elapsed = (end - start) * 1000
    print(f"{'local-' + nameA:<20} {lang:<15} {elapsed:.4f} ms")

    start = time.time()
    al = affine_alignment(seqA, seqB)
    # print(''.join(al[0]), ''.join(al[1]), sep='\n')
    end = time.time()
    elapsed = (end - start) * 1000
    print(f"{'affine-' + nameA:<20} {lang:<15} {elapsed:.4f} ms")

    start = time.time()
    # print(seqA, seqB)
    al = fitting(seqA, seqB)
    # print(al)
    # print(''.join(al[0]), ''.join(al[1]), sep='\n')
    end = time.time()
    elapsed = (end - start) * 1000
    print(f"{'fitting-' + nameA:<20} {lang:<15} {elapsed:.4f} ms")