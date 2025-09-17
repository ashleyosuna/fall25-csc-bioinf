#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin

CODE_DIRECTORY="week1/code"
DATA_DIRECTORY="week1/data"
REPOSITORY="week1/resources"
TESTS_DIRECTORY="week1/tests"
data="data1 data2 data3 data4"

ulimit -s 8192000

printf "%-10s %-10s %-10s %-10s\n" "Dataset" "Language" "Runtime" "N50"
echo "---------------------------------------------------------"

for d in $data; do
    #### python ####
    start=$(date +%s)
    
    # redirect output to a file recording lengths of the contigs
    python3 $REPOSITORY/main.py $DATA_DIRECTORY/$d > $TESTS_DIRECTORY/lengths_python_$d.txt
    
    end=$(date +%s)
    
    elapsed=$(( end - start ))
    
    formatted=$(printf "%d:%02d:00" $((elapsed/60)) $((elapsed%60)))
    
    #calculating n50
    n50_python=$(codon run -release $CODE_DIRECTORY/n50.py $TESTS_DIRECTORY/lengths_python_$d.txt)
    printf "%-10s %-10s %-10s %-10d\n" ${d} "python" ${formatted} ${n50_python}

    
    #### codon ####
    start=$(date +%s)

    # redirect output to a file recording lengths of the contigs
    codon run -release $CODE_DIRECTORY/main.py $DATA_DIRECTORY/$d > $TESTS_DIRECTORY/lengths_codon_$d.txt
    end=$(date +%s)

    elapsed=$((end - start))
    formatted=$(printf "%d:%02d:00" $((elapsed/60)) $((elapsed%60)))

    # calculating n50
    n50_codon=$(codon run -release $CODE_DIRECTORY/n50.py $TESTS_DIRECTORY/lengths_codon_$d.txt)
    printf "%-10s %-10s %-10s %-10d\n" ${d} "codon" ${formatted} ${n50_codon}
done