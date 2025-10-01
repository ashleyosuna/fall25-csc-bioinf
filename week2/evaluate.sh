#!/bin/bash
# set -euxo pipefail
PATH=${PATH}:${HOME}/.codon/bin
echo "week2"

printf "\nRunning codon tests\n"
codon run -release week2/code/tests_motif.py

printf "\nRunning python tests\n"
python3 week2/code/tests_python.py