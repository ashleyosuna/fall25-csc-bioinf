#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

printf "%-10s %-10s\n" "Language" "Runtime"
echo "----------------------------"
echo "codon "
codon run -release week3/code/test_codon.py

echo "----------------------------"
echo "python "
python3 week3/code/test_python.py