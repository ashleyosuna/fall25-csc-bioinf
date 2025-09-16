#!/bin/bash
data="data1 data2 data3 data4"

printf "%-10s %-10s %-10s %-10s\n" "Dataset" "Language" "Runtime" "N50"
echo "---------------------------------------------------------"

for d in $data; do
    # if [[ "$d" == "data4" ]]; then
    #     echo "hello"
    #     ulimit -s 8192000
    # fi
    # python
    start=$(date +%s)
    
    python3 ./genome-assembly-repo/main.py ./data/$d > ./tests/lengths_python_$d.txt
    
    end=$(date +%s)
    
    elapsed=$(( end - start ))
    
    formatted=$(printf "%d:%02d:00" $((elapsed/3600)) $((elapsed%3600)))
    
    # calculating n50
    n50_python=$(codon run -release ./tests/n50.py ./tests/lengths_python_$d.txt)
    printf "%-10s %-10s %-10s %-10d\n" ${d} "python" ${formatted} ${n50_python}

    
    # codon
    start=$(date +%s)
    codon run -release ./code/main.py ./data/$d > ./tests/lengths_codon_$d.txt
    end=$(date +%s)

    elapsed=$((end - start))
    formatted=$(printf "%d:%02d:00" $((elapsed/3600)) $((elapsed%3600)))

    # calculating n50
    n50_codon=$(codon run -release ./tests/n50.py ./tests/lengths_codon_$d.txt)
    printf "%-10s %-10s %-10s %-10d\n" ${d} "codon" ${formatted} ${n50_codon}
done