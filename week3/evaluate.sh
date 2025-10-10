#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week3"

echo -e "\nRunning codon tests\n"
codon run -release week3/code/test_codon.py

echo "----------------------------"
echo -e "\nRunning python tests\n"
python3 week3/code/test_python.py