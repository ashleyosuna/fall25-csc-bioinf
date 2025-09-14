import sys

argv = sys.argv

len_sum = 0
lengths = []

with open(argv[1], 'r') as f:
    lines = f.readlines()
    lines = lines[3:]

    for line in lines:
        _, length = line.strip().split()
        lengths.append(length)
        len_sum += int(length)
    
middle = len_sum / 2
len_sum = 0

for length in lengths:
    len_sum += int(length)

    if len_sum >= middle:
        print(length)
        break