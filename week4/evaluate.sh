#!/bin/bash
set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

printf "%-16s %-16s %-16s\n" "Method" "Language" "Runtime"
echo "----------------------------"
# time=$(codon run -release week3/code/test_codon.py)
# printf "%-10s %-10s %-2s\n" "Codon" ${time}

# time=$(python3 week3/code/test_python.py)
# printf "%-10s %-10s %-2s\n" "Python" ${time}
python3 week4/code/python/tests.py