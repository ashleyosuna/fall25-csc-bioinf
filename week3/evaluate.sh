#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

printf "%-10s %-10s\n" "Language" "Runtime"
echo "----------------------------"

start=$(date +%s%3N)
# echo -e "\nRunning codon tests\n"
tmp=$(codon run -release week3/code/test_codon.py)

end_time=$(date +%s%3N)

duration_ms=$((end_time - start_time))
printf "%-10s %-10f\n" "codon" ${duration_ms}

echo "----------------------------"
echo -e "\nRunning python tests\n"
python3 week3/code/test_python.py