#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

printf "%-10s %-10s\n" "Language" "Runtime"
echo "----------------------------"
time=$(codon run -release week3/code/test_codon.py)
printf "%-10s %-10s\n" "Codon" ${time}

time=$(python3 week3/code/test_python.py)
printf "%-10s %-10s\n" "Python" ${time}